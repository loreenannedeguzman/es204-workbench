"""General ES204 workbench: Tk widgets + a JSON calculation worker.

Read numerics.py for the formulas. This file handles user input and display.
"""
import csv
import io
import json
from pathlib import Path
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from tkinter.scrolledtext import ScrolledText

APP_DIR=Path(__file__).resolve().parent
METHODS=['Euler','Predictor-Corrector','RK2.1 (Heun)','RK2.2 (lesson Ralston)','RK3','RK4']


def enable_dpi():
    # Match Tk window coordinates to the display's pixels on scaled Windows desktops.
    if sys.platform=='win32':
        import ctypes
        try:ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except (OSError,AttributeError):pass


class ScrollForm(ttk.Frame):
    """Scrollable inputs so all fields remain reachable on smaller screens."""
    def __init__(self,parent):
        super().__init__(parent)
        canvas=tk.Canvas(self,highlightthickness=0,bg='#edf2f7',width=390)
        scrollbar=ttk.Scrollbar(self,command=canvas.yview)
        scrollbar.pack(side='right',fill='y')
        canvas.pack(side='left',fill='both',expand=True)
        canvas.configure(yscrollcommand=scrollbar.set)
        self.body=ttk.Frame(canvas,padding=12)
        window=canvas.create_window((0,0),window=self.body,anchor='nw')
        self.body.bind('<Configure>',lambda e:canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>',lambda e:canvas.itemconfigure(window,width=e.width))


class App:
    def __init__(self,root):
        self.root=root
        self.process=None
        self.events=queue.Queue()
        self.generation=0
        self.datasets=[]
        self.result=None
        self.history=[]
        self.fields={}
        self.texts={}
        root.title('ES204 • Numerical Workbench')
        root.geometry(f'{min(1280,root.winfo_screenwidth()-80)}x{min(850,root.winfo_screenheight()-100)}')
        root.minsize(880,580)
        style=ttk.Style();style.theme_use('clam')
        style.configure('.',font=('Segoe UI',10))
        style.configure('TFrame',background='#edf2f7')
        style.configure('TLabel',background='#edf2f7',foreground='#142d47')
        style.configure('Accent.TButton',background='#155e75',foreground='white',padding=7)
        header=ttk.Frame(root,padding=12);header.pack(fill='x')
        ttk.Label(header,text='ES204 / Numerical Workbench',font=('Segoe UI',18,'bold')).pack(side='left')
        self.run_button=ttk.Button(header,text='Solve ODE',style='Accent.TButton',command=self.run)
        self.run_button.pack(side='right',padx=6)
        ttk.Button(header,text='Cancel',command=self.cancel).pack(side='right')
        tools=ttk.Frame(root,padding=(12,0,12,8));tools.pack(fill='x')
        for label,action in [('Save inputs…',self.save_project),('Load inputs…',self.load_project),('Save report…',self.save_report),('Help',self.help)]:
            ttk.Button(tools,text=label,command=action).pack(side='left',padx=(0,6))
        self.status=tk.StringVar(value='Enter inputs, choose methods, then Solve ODE.')
        ttk.Label(root,textvariable=self.status,wraplength=1000,padding=(12,0,12,6)).pack(fill='x')
        panes=ttk.Panedwindow(root,orient='horizontal');panes.pack(fill='both',expand=True,padx=12,pady=(0,12))
        self.tabs=ttk.Notebook(panes);panes.add(self.tabs,weight=0)
        for name in ['ODE','Data / analysis','Roots']:
            frame=ScrollForm(self.tabs);self.tabs.add(frame,text=name)
        self.ode_form=self.tabs.nametowidget(self.tabs.tabs()[0]).body
        self.data_form=self.tabs.nametowidget(self.tabs.tabs()[1]).body
        self.root_form=self.tabs.nametowidget(self.tabs.tabs()[2]).body
        self.build_ode();self.build_analysis();self.build_roots()
        self.field_defaults={k:v.get() for k,v in self.fields.items()}
        self.tabs.bind('<<NotebookTabChanged>>',self.tab_changed)
        right=ttk.Frame(panes);panes.add(right,weight=1)
        self.output_tabs=ttk.Notebook(right);self.output_tabs.pack(fill='both',expand=True)
        self.table_frame=ttk.Frame(self.output_tabs);self.output_tabs.add(self.table_frame,text='Tables')
        self.plot_frame=ttk.Frame(self.output_tabs);self.output_tabs.add(self.plot_frame,text='Plots')
        notes_frame=ttk.Frame(self.output_tabs);self.output_tabs.add(notes_frame,text='Report / formulas')
        self.report=ScrolledText(notes_frame,wrap='word',font=('Consolas',10),state='disabled')
        self.report.pack(fill='both',expand=True)
        tablebar=ttk.Frame(self.table_frame);tablebar.pack(fill='x',pady=6)
        self.table_choice=ttk.Combobox(tablebar,state='readonly');self.table_choice.pack(side='left',fill='x',expand=True)
        self.table_choice.bind('<<ComboboxSelected>>',lambda e:self.show_table())
        ttk.Button(tablebar,text='Export CSV…',command=self.export_csv).pack(side='right',padx=4)
        ttk.Button(tablebar,text='Copy',command=self.copy_table).pack(side='right')
        ttk.Button(tablebar,text='Use as data…',command=lambda:self.table_picker(False)).pack(side='right')
        grid=ttk.Frame(self.table_frame);grid.pack(fill='both',expand=True)
        self.tree=ttk.Treeview(grid,show='headings')
        vs=ttk.Scrollbar(grid,orient='vertical',command=self.tree.yview)
        hs=ttk.Scrollbar(grid,orient='horizontal',command=self.tree.xview)
        self.tree.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
        grid.rowconfigure(0,weight=1);grid.columnconfigure(0,weight=1)
        self.tree.grid(row=0,column=0,sticky='nsew');vs.grid(row=0,column=1,sticky='ns');hs.grid(row=1,column=0,sticky='ew')
        self.chart_choice=ttk.Combobox(self.plot_frame,state='readonly');self.chart_choice.pack(fill='x',pady=6)
        self.chart_choice.bind('<<ComboboxSelected>>',lambda e:self.show_chart())
        plot_options=ttk.Frame(self.plot_frame);plot_options.pack(fill='x')
        self.plot_options={}
        for label,choices in [('Error',['Signed','Absolute']),('X scale',['Auto','Linear','Log','Symlog']),('Y scale',['Auto','Linear','Log','Symlog'])]:
            ttk.Label(plot_options,text=label).pack(side='left')
            var=tk.StringVar(value=choices[0]);self.plot_options[label]=var
            combo=ttk.Combobox(plot_options,textvariable=var,values=choices,state='readonly',width=9)
            combo.pack(side='left',padx=3);combo.bind('<<ComboboxSelected>>',lambda e:self.show_chart())
        self.plot_tolerance=tk.BooleanVar(value=True)
        ttk.Checkbutton(plot_options,text='Tolerance lines',variable=self.plot_tolerance,command=self.show_chart).pack(side='left')
        ttk.Button(self.plot_frame,text='Build custom versus plot…',command=lambda:self.table_picker(True)).pack(fill='x')
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
        self.figure=Figure(figsize=(6,4),dpi=100,layout='constrained')
        self.canvas=FigureCanvasTkAgg(self.figure,master=self.plot_frame)
        toolbar=NavigationToolbar2Tk(self.canvas,self.plot_frame,pack_toolbar=False)
        toolbar.pack(side='bottom',fill='x')
        self.canvas.get_tk_widget().pack(fill='both',expand=True)
        # A hidden notebook tab can map after Matplotlib changes its DPI.
        # Fit the backing image to the actual canvas once the tab is visible.
        self.canvas.get_tk_widget().bind('<Map>',lambda event:root.after_idle(self.fit_chart),add='+')
        self.set_report('No built-in problem solutions. Enter your own ODE or data.\n\nPlots and post-processing use the results of your inputs.\nUse Help for the workflow.')
        root.protocol('WM_DELETE_WINDOW',self.close)
        root.after(60,self.poll)

    def hint(self,parent,text):
        ttk.Label(parent,text=text,wraplength=355).pack(fill='x',pady=(2,8))

    def field(self,parent,key,label,default='',choices=None):
        ttk.Label(parent,text=label,wraplength=355).pack(anchor='w')
        var=tk.StringVar(value=default);self.fields[key]=var
        widget=ttk.Combobox(parent,textvariable=var,values=choices,state='readonly') if choices else ttk.Entry(parent,textvariable=var)
        widget.pack(fill='x',pady=(2,8))
        return widget

    def text(self,parent,key,label,height=5):
        ttk.Label(parent,text=label,wraplength=350).pack(anchor='w')
        widget=ScrolledText(parent,height=height,width=36,font=('Consolas',10),wrap='none')
        widget.pack(fill='x',pady=(2,8));self.texts[key]=widget

    def build_ode(self):
        p=self.ode_form
        self.hint(p,"First-order IVP: y′ = f(x,y). Enter the right side only. Use x, y, sin, cos, exp, sqrt, log; powers ^ or **.")
        self.field(p,'formula',"ODE f(x,y)")
        self.field(p,'output_mode','Output mode','Full tables',choices=['Full tables','Final values only'])
        self.field(p,'grid_mode','Specify the grid using','Step sizes h',choices=['Step sizes h','Interval counts n','Refinement h0 / 2^m','Search n (final error)','Powers of ten n=10^p'])
        self.field(p,'first_p','First p (powers of ten mode)','1')
        self.field(p,'last_p','Last p (powers of ten mode; up to 7)','6')
        self.hint(p,'Powers-of-ten mode produces final-value tables automatically. Final-values mode supports up to 10 million steps per grid without storing every node.')
        self.fields['grid_mode'].trace_add('write',lambda *_: self.fields['output_mode'].set('Final values only') if self.fields['grid_mode'].get()=='Powers of ten n=10^p' else None)
        for key,label,default in [('x0','Initial x','0'),('xf','Final x',''),('y0','Initial y',''),('steps','h values / n values / single h0 (per grid mode)',''),('max_m','Last m (halving mode: m=0 through this value)','3'),('iterations','PC passes: one value or list, e.g. 1,2,3,4','1')]:
            self.field(p,key,label,default)
        self.field(p,'max_n','Maximum n for final-error search','100000')
        self.hint(p,'n-search uses the h/n field as its starting n and assumes decreasing final error. It reports a candidate minimum, with trial results.')
        self.field(p,'stop_metric','Stop at first tested grid meeting tolerance','No early stop',choices=['No early stop','Maximum true error %','Final true error %','Maximum successive error %','Final successive error %'])
        self.field(p,'comparison','Stopping comparison','<',choices=['<','<='])
        self.field(p,'final_reference_source','Final-values reference source','Exact formula if provided',choices=['Exact formula if provided','Given final reference table'])
        self.field(p,'stop_tolerance','Error tolerance in percent (for early stop)','1')
        ttk.Label(p,text='Methods (select one or several)').pack(anchor='w')
        self.method_vars={}
        for method in METHODS:
            var=tk.BooleanVar(value=method=='RK4');self.method_vars[method]=var
            ttk.Checkbutton(p,text=method,variable=var).pack(anchor='w')
        self.field(p,'exact','Exact solution formula y(x) — optional')
        self.text(p,'reference','Given exact values — optional. One x, y pair per line.',6)
        ttk.Button(p,text='Use exact-value table as analysis data',command=self.use_references).pack(fill='x')
        self.field(p,'requested','Report at these x-values (optional, comma-separated)')
        self.field(p,'sampling','Off-grid evaluation policy','Grid nodes only',choices=['Grid nodes only','Linear interpolation','Polynomial interpolation'])
        self.hint(p,'A single final reference is one row: xf, y(xf). Leave the formula blank if only tabulated exact values are provided. Missing grid nodes are flagged, not interpolated.')

    def build_analysis(self):
        p=self.data_form
        self.hint(p,'Use a computed solution or paste your own data. Edit the rows to select the points required by your problem.')
        self.dataset_choice=ttk.Combobox(p,state='readonly');self.dataset_choice.pack(fill='x',pady=4)
        self.field(p,'data_select','Copy only these x nodes (optional, comma-separated)')
        ttk.Button(p,text='Use selected solution as data',command=self.use_dataset).pack(fill='x')
        ttk.Button(p,text='Import two-column CSV…',command=self.import_data).pack(fill='x',pady=4)
        self.text(p,'data','Data points: x, y per line (no header)',8)
        self.field(p,'data_source','Input source','Data table',choices=['Data table','Sample formula'])
        self.field(p,'data_formula','Function to sample (Sample formula mode)')
        for key,label,value in [('data_min','Sampling lower x',''),('data_max','Sampling upper x',''),('data_n','Sampling number of intervals','10')]:self.field(p,key,label,value)
        self.field(p,'operation','Analysis','Polynomial fit',choices=['Polynomial fit','Stationary point (secant)','Finite differences','Trapezoidal rule',"Simpson's 1/3 rule",'Linearized least squares'])
        self.field(p,'fit','Polynomial method','Exact Lagrange',choices=['Exact Lagrange','Exact Vandermonde','Least squares'])
        self.field(p,'degree','Least-squares degree','2')
        self.field(p,'transformed_model','Linearized least-squares model','Exponential: a*exp(b*x)',choices=['Exponential: a*exp(b*x)','Power: a*x^b','Reciprocal: 1/(a+b*x)','Inverse x: a+b/x'])
        self.field(p,'scheme','Finite-difference scheme','Central + higher-order endpoints',choices=['Central + higher-order endpoints','Central only','Forward','Backward','Forward higher-order','Backward higher-order'])
        self.field(p,'analysis_exact','Exact function y(x) for derivative/fit comparisons')
        self.field(p,'exact_integral','Given exact integral — optional')
        self.field(p,'analysis_min','Polynomial plot/integral lower x (blank = data bound)')
        self.field(p,'analysis_max','Polynomial plot/integral upper x (blank = data bound)')
        self.field(p,'query','Evaluate polynomial at x-values (optional, comma-separated)')
        self.hint(p,'Polynomial fitting uses all rows currently in the table. Exact interpolation has degree at most n−1. Polynomial integration uses the entire data domain.')
        self.field(p,'a','Secant first guess')
        self.field(p,'b','Secant second guess')
        self.field(p,'tolerance','Secant residual tolerance |f′(x)|','1e-5')
        self.field(p,'limit','Maximum secant updates','100')

    def build_roots(self):
        p=self.root_form
        self.hint(p,'Solve g(x)=0 for a formula you enter. For a fitted-polynomial minimum, use Data / analysis → Stationary point instead.')
        self.field(p,'root_function','g(x)')
        self.field(p,'root_method','Root method','Secant',choices=['Secant','Bisection','Newton-Raphson'])
        self.field(p,'root_criterion','Stopping criterion','Residual',choices=['Residual','Relative step change %'])
        for key,label,default in [('a','First guess / lower bracket',''),('b','Second guess / upper bracket (unused by Newton)',''),('tolerance','Residual tolerance |g(x)|','1e-5'),('limit','Maximum iterations','100'),('plot_min','Plot lower x',''),('plot_max','Plot upper x','')]:
            self.field(p,'root_'+key,label,default)

    def tab_changed(self,event=None):
        self.run_button.configure(text=['Solve ODE','Analyze data','Find root'][self.tabs.index('current')])

    def inputs(self):
        return {**{k:v.get() for k,v in self.fields.items()},**{k:v.get('1.0','end-1c') for k,v in self.texts.items()},'methods':[k for k,v in self.method_vars.items() if v.get()]}

    def run(self):
        if self.process is not None: return
        raw=self.inputs();tab=self.tabs.index('current')
        action=['solve','analyze','roots'][tab]
        data=raw if tab<2 else {k[5:]:v for k,v in raw.items() if k.startswith('root_')}
        if tab==2 and data['method']=='Newton-Raphson': data['b']=data['a']
        executable=Path(sys.executable).with_name('python.exe')
        try:
            process=subprocess.Popen([str(executable),str(APP_DIR/'numerics.py')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,encoding='utf-8',creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        except OSError as error:
            messagebox.showerror('Could not start',str(error));return
        self.process=process;generation=self.generation
        self.run_button.configure(state='disabled');self.status.set('Calculating… Cancel remains available.')
        request={'action':action,'inputs':data}
        def work():
            try:
                stdout,stderr=process.communicate(json.dumps(request))
                response=json.loads(stdout) if stdout.strip() else {'ok':False,'error':stderr or 'Worker ended without a result.'}
            except Exception as error: response={'ok':False,'error':str(error)}
            self.events.put((generation,response,action))
        threading.Thread(target=work,daemon=True).start()

    def poll(self):
        try:
            generation,response,action=self.events.get_nowait()
        except queue.Empty: pass
        else:
            if generation==self.generation:
                self.process=None;self.run_button.configure(state='normal')
                if not response['ok']:
                    self.status.set('Input or calculation error — previous results retained.')
                    messagebox.showerror('Check inputs',response['error'])
                else:
                    self.render(response['result'])
                    if action=='solve':
                        self.datasets=self.result['datasets']
                        self.dataset_choice.configure(values=[d['label'] for d in self.datasets])
                        if self.datasets: self.dataset_choice.current(0)
                    self.status.set('Complete. Tables and plots reflect the inputs at the last successful run.')
                    if self.result.get('warnings'):
                        self.status.set('Complete with reference-precision warning. See Report / formulas.')
                        messagebox.showwarning('Reference precision', '\n\n'.join(self.result['warnings']))
        self.root.after(60,self.poll)

    @staticmethod
    def display(value):
        if value is None:return 'N/A'
        if isinstance(value,float):return f'{value:.12g}'
        return str(value)

    def render(self,result):
        self.result=result
        self.history.extend(result['tables'])
        self.table_choice.configure(values=[t['title'] for t in result['tables']])
        if result['tables']:self.table_choice.current(0)
        self.show_table()
        self.chart_choice.configure(values=[c['title'] for c in result['charts']])
        if result['charts']:self.chart_choice.current(0)
        self.show_chart()
        self.set_report(self.report_text())

    def show_table(self):
        self.tree.delete(*self.tree.get_children())
        if not self.result or self.table_choice.current()<0:return
        table=self.result['tables'][self.table_choice.current()]
        names=[str(i) for i in range(len(table['columns']))];self.tree.configure(columns=names)
        for name,label in zip(names,table['columns']):
            self.tree.heading(name,text=label);self.tree.column(name,width=max(115,min(240,len(label)*9)),stretch=False)
        # Keep the interface responsive on very large grids; exports retain all rows.
        for row in table['rows'][:5000]:self.tree.insert('', 'end',values=[self.display(v) for v in row])
        if len(table['rows'])>5000:self.tree.insert('','end',values=['First 5000 rows shown; export CSV for all rows.'])

    def show_chart(self):
        self.figure.clear()
        if self.result and self.chart_choice.current()>=0:
            chart=self.result['charts'][self.chart_choice.current()];axis=self.figure.add_subplot(111)
            from plotting import draw_chart
            draw_chart(axis,chart,self.plot_options['Error'].get(),self.plot_options['X scale'].get(),self.plot_options['Y scale'].get(),self.plot_tolerance.get())
        self.canvas.draw_idle()

    def fit_chart(self):
        from types import SimpleNamespace
        widget=self.canvas.get_tk_widget()
        if widget.winfo_width()>10 and widget.winfo_height()>10:
            self.canvas.resize(SimpleNamespace(width=widget.winfo_width(),height=widget.winfo_height()))

    def report_text(self):
        if not self.result:return ''
        parts=['\n\n'.join(self.result['notes'])]
        for table in self.result['tables']:
            parts += ['\n'+table['title'],'\t'.join(table['columns'])]
            parts.extend('\t'.join(self.display(v) for v in row) for row in table['rows'])
        return '\n'.join(parts)

    def set_report(self,text):
        self.report.configure(state='normal');self.report.delete('1.0','end');self.report.insert('1.0',text);self.report.configure(state='disabled')

    def use_dataset(self):
        index=self.dataset_choice.current()
        if index<0:return
        data=self.datasets[index]
        pairs=list(zip(data['x'],data['y']))
        if self.fields['data_select'].get().strip():
            try:
                from numerics import number
                import math
                targets=[number(v.strip()) for v in self.fields['data_select'].get().split(',')]
                selected=[]
                for target in targets:
                    pair=next((p for p in pairs if math.isclose(p[0],target,rel_tol=1e-9,abs_tol=1e-10)),None)
                    if pair is None:raise ValueError(f'x={target:g} is not a node of this solution. Use ODE requested-x evaluation to interpolate explicitly.')
                    selected.append(pair)
                pairs=selected
            except (ValueError,TypeError) as error:messagebox.showerror('Node selection',str(error));return
        self.texts['data'].delete('1.0','end')
        self.texts['data'].insert('1.0','\n'.join(f'{x:.17g}, {y:.17g}' for x,y in pairs))
        self.fields['data_source'].set('Data table')
        self.status.set('Copied '+data['label']+' from the last ODE run. Edit rows to choose a subset.')

    def use_references(self):
        self.texts['data'].delete('1.0','end')
        self.texts['data'].insert('1.0',self.texts['reference'].get('1.0','end-1c'))
        self.fields['data_source'].set('Data table');self.tabs.select(1)
        self.status.set('Exact-value table copied as analysis data.')

    def table_picker(self,plot):
        if not self.history:return
        dialog=tk.Toplevel(self.root);dialog.title('Custom versus plot' if plot else 'Choose result columns for analysis')
        dialog.geometry('720x620' if plot else '720x320')
        container=ScrollForm(dialog);container.pack(fill='both',expand=True);p=container.body
        ttk.Label(p,text='Choose any table from this session, then its numeric columns.').pack(anchor='w')
        choice=ttk.Combobox(p,state='readonly',values=[f'{i+1}. {t["title"]}' for i,t in enumerate(self.history)])
        choice.pack(fill='x',pady=6)
        xchoice=ttk.Combobox(p,state='readonly');ychoice=ttk.Combobox(p,state='readonly')
        ttk.Label(p,text='Horizontal axis / x column').pack(anchor='w');xchoice.pack(fill='x')
        ttk.Label(p,text='Vertical axis / y column').pack(anchor='w');ychoice.pack(fill='x')
        def columns(event=None):
            cols=self.history[choice.current()]['columns']
            names=[f'{i}: {c}' for i,c in enumerate(cols)]
            xchoice.configure(values=names);ychoice.configure(values=names)
            xchoice.current(0);ychoice.current(min(1,len(names)-1))
        choice.bind('<<ComboboxSelected>>',columns);choice.current(len(self.history)-1);columns()
        def pairs():
            table=self.history[choice.current()];ix=xchoice.current();iy=ychoice.current()
            values=[(r[ix],r[iy]) for r in table['rows'] if isinstance(r[ix],(int,float)) and not isinstance(r[ix],bool) and isinstance(r[iy],(int,float)) and not isinstance(r[iy],bool)]
            if not values:raise ValueError('No numeric pairs in these columns. N/A rows are omitted.')
            return table,values
        if not plot:
            def use():
                try: table,values=pairs()
                except ValueError as error:messagebox.showerror('Columns',str(error));return
                self.texts['data'].delete('1.0','end');self.texts['data'].insert('1.0','\n'.join(f'{x:.17g}, {y:.17g}' for x,y in values))
                self.fields['data_source'].set('Data table');self.tabs.select(1);dialog.destroy()
            ttk.Button(p,text='Use these columns as data',command=use).pack(pady=12)
            return
        ttk.Label(p,text='Series label').pack(anchor='w');label=tk.StringVar();ttk.Entry(p,textvariable=label).pack(fill='x')
        curves=[];listing=tk.Listbox(p,height=5);listing.pack(fill='x',pady=6)
        def add():
            try: table,values=pairs()
            except ValueError as error:messagebox.showerror('Columns',str(error));return
            name=label.get() or table['title']+' / '+table['columns'][ychoice.current()]
            curves.append({'label':name,'x':[v[0] for v in values],'y':[v[1] for v in values],'style':'nodes'})
            listing.insert('end',name)
        ttk.Button(p,text='Add series (repeat for an overlay)',command=add).pack(fill='x')
        settings={}
        for key,value in [('Title','Custom comparison'),('X label','x'),('Y label','y')]:
            row=ttk.Frame(p);row.pack(fill='x',pady=3);ttk.Label(row,text=key,width=10).pack(side='left')
            settings[key]=tk.StringVar(value=value);ttk.Entry(row,textvariable=settings[key]).pack(side='left',fill='x',expand=True)
        xs=tk.BooleanVar();ys=tk.BooleanVar()
        ttk.Checkbutton(p,text='Logarithmic X axis (positive values only)',variable=xs).pack(anchor='w')
        ttk.Checkbutton(p,text='Logarithmic Y axis (positive values only)',variable=ys).pack(anchor='w')
        def draw():
            if not curves:add()
            if not curves:return
            if (xs.get() and any(v<=0 for c in curves for v in c['x'])) or (ys.get() and any(v<=0 for c in curves for v in c['y'])):
                messagebox.showerror('Log axis','All plotted values on a logarithmic axis must be positive.');return
            chart={'title':settings['Title'].get(),'xlabel':settings['X label'].get(),'ylabel':settings['Y label'].get(),
                   'curves':curves,'xscale':'log' if xs.get() else 'linear','yscale':'log' if ys.get() else 'linear'}
            if self.result is None:self.result={'tables':[],'charts':[],'notes':[]}
            self.result['charts'].append(chart);self.chart_choice.configure(values=[c['title'] for c in self.result['charts']]);self.chart_choice.current(len(self.result['charts'])-1)
            self.show_chart();self.output_tabs.select(self.plot_frame);dialog.destroy()
        ttk.Button(p,text='Plot these series',command=draw).pack(fill='x',pady=8)

    def import_data(self):
        path=filedialog.askopenfilename(filetypes=[('CSV or text','*.csv *.txt'),('All files','*.*')])
        if path:
            try: value=Path(path).read_text(encoding='utf-8-sig')
            except OSError as error:messagebox.showerror('Import failed',str(error));return
            self.texts['data'].delete('1.0','end');self.texts['data'].insert('1.0',value)

    def selected_csv(self,delimiter=','):
        if not self.result or self.table_choice.current()<0:return ''
        table=self.result['tables'][self.table_choice.current()]
        stream=io.StringIO();writer=csv.writer(stream,delimiter=delimiter)
        writer.writerow(table['columns']);writer.writerows(table['rows']);return stream.getvalue()

    def copy_table(self):
        self.root.clipboard_clear();self.root.clipboard_append(self.selected_csv('\t'))

    def export_csv(self):
        if not self.result:return
        path=filedialog.asksaveasfilename(defaultextension='.csv',filetypes=[('CSV','*.csv')])
        if path:Path(path).write_text(self.selected_csv(),encoding='utf-8',newline='')

    def save_report(self):
        if not self.result:return
        path=filedialog.asksaveasfilename(defaultextension='.txt',filetypes=[('Text report','*.txt')])
        if path:Path(path).write_text(self.report_text(),encoding='utf-8')

    def save_project(self):
        path=filedialog.asksaveasfilename(defaultextension='.json',filetypes=[('Workbench inputs','*.json')])
        if path:Path(path).write_text(json.dumps({'version':2,'inputs':self.inputs()},indent=2),encoding='utf-8')

    def load_project(self):
        if self.process is not None:return
        path=filedialog.askopenfilename(filetypes=[('Workbench inputs','*.json')])
        if not path:return
        try:
            content=json.loads(Path(path).read_text(encoding='utf-8'))
            if content.get('version')!=2:raise ValueError('Unsupported input-file version.')
            data=content['inputs']
            if not isinstance(data,dict):raise ValueError('Invalid inputs file.')
        except (OSError,ValueError,KeyError) as error:messagebox.showerror('Load failed',str(error));return
        for k,v in self.fields.items():v.set(data.get(k,self.field_defaults.get(k,'')))
        for k,w in self.texts.items():w.delete('1.0','end');w.insert('1.0',data.get(k,''))
        for k,v in self.method_vars.items():v.set(k in data.get('methods',[]))
        self.datasets=[];self.dataset_choice.set('');self.dataset_choice.configure(values=[])
        self.result=None;self.history=[];self.tree.delete(*self.tree.get_children());self.table_choice.set('');self.table_choice.configure(values=[])
        self.chart_choice.set('');self.chart_choice.configure(values=[]);self.figure.clear();self.canvas.draw_idle()
        self.set_report('Inputs loaded. Run the calculation to generate results.')
        self.status.set('Inputs loaded; old results cleared. Choose the appropriate tab and run.')

    def cancel(self):
        self.generation+=1
        if self.process and self.process.poll() is None:self.process.terminate()
        self.process=None;self.run_button.configure(state='normal');self.status.set('Calculation cancelled. Previous results retained.')

    def close(self):
        self.cancel();self.root.destroy()

    def help(self):
        messagebox.showinfo('Workflow',
            'ODE: enter a first-order ODE, initial condition, interval and one or more step sizes. Select methods and Solve.\n\n'
            'For large runs: choose Final values only. Powers of ten mode tests n=10^p from First p through Last p and stops separately per method. Select Final true error %, a tolerance in percent, and < or <=. Signed and absolute errors are both shown.\n\n'
            'Final-values reference source lets you explicitly choose your formula or table. A precision warning never changes your reference. Final-only runs keep no node datasets for post-processing.\n\n'
            'Exact references: paste x,y pairs in the exact-values table. A formula is optional. Tables include errors at matching nodes. Multiple grids also produce successive-grid differences.\n\n'
            'Data / analysis: copy a saved ODE solution or paste data. Keep only the required points. Choose interpolation, least squares, finite differences, integration, or a stationary point of the fitted polynomial.\n\n'
            'Roots: enter g(x), a method and starting values. Convergence uses the residual tolerance.\n\n'
            'Plots are embedded. Use the plot toolbar to zoom or save an image. Use CSV for numerical tables. Save inputs stores your setup, not computed results.\n\n'
            'Scroll the input panel to reach all fields. Changing inputs does not recompute results until Run. This app handles scalar first-order IVPs, not systems, boundary-value problems or automatic interpretation of problem text.')


if __name__=='__main__':
    enable_dpi()
    root=tk.Tk()
    try:App(root);root.mainloop()
    except Exception as error:messagebox.showerror('App error',str(error));root.destroy()
