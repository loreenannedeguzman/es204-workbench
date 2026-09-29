"""Forms for general systems, FDM initial-value problems and comparisons."""
import math
import tkinter as tk
from tkinter import ttk, messagebox


def build(app, ScrollForm, methods):
    for title,prefix in [('Systems','systems'),('FDM IVP','fdm_ivp'),('Compare','comparison')]:
        form=ScrollForm(app.tabs);app.tabs.add(form,text=title);p=form.body
        ttk.Button(p,text='Show equations',command=lambda:app.run(equations_only=True)).pack(fill='x',pady=4)
        if prefix=='comparison':
            app.hint(p,'Compare one or more numerical series against supplied reference values at common x nodes. The first column is x, the second is the reference; subsequent columns are numerical series.')
            app.text(p,'comparison_data','Header required: x, reference, series1, series2, ...',12)
            app.hint(p,'Start with just x, reference and their rows, then add numerical columns from previous results below. Or paste a complete table with headers.')
            ttk.Button(p,text='Add column from result history…',command=lambda:history_column(app)).pack(fill='x')
            continue
        if prefix=='systems':
            app.hint(p,'One formula per dependent variable. First order: RHS means its first derivative. Second order: RHS means its second derivative; dNAME denotes the first derivative (e.g. du).')
            app.field(p,'systems_order','System order','First order',choices=['First order','Second order'])
            app.field(p,'systems_names','Dependent-variable names, in order (e.g. u, v)')
            app.text(p,'systems_rhs','Right-hand sides, one per variable (no equals sign)',5)
            app.field(p,'systems_initial','Initial values: first order u,v; second order u,du,v,dv')
            app.field(p,'systems_method','Step method','RK4',choices=methods)
            app.field(p,'systems_passes','PC corrector passes','1')
        else:
            app.hint(p,"Linear second-order IVP: a(x)y'' + b(x)y' + c(x)y = d(x). Both initial conditions are at the left endpoint. Coefficients depend only on x.")
            for key,label,value in [('a',"Coefficient a(x) of y''",'1'),('b',"Coefficient b(x) of y'",'0'),('c','Coefficient c(x) of y','0'),('d','Right side d(x)',''),('y0','Initial y',''),('v0',"Initial y'",'')]:
                app.field(p,prefix+'_'+key,label,value)
            app.field(p,'fdm_ivp_start','Starting formula','Central ghost / Taylor (second order)',choices=['Central ghost / Taylor (second order)','Forward derivative (first order)','Forward derivative (second order)'])
            app.hint(p,'Choose the starting formula required by your lesson. The report shows it explicitly. A first-order start can reduce overall accuracy even with central interior differences.')
        for key,label,value in [('left','Initial x',''),('right','Final x','')]:app.field(p,prefix+'_'+key,label,value)
        app.field(p,prefix+'_grid','Grid definition','Step h',choices=['Step h','Intervals'])
        app.field(p,prefix+'_size','Grid value (h or number of intervals)','')
        if prefix=='systems':
            app.field(p,'systems_plot_x','Optional versus plot: horizontal state name or x')
            app.field(p,'systems_plot_y','Optional versus plot: vertical state name')
            app.hint(p,'Solution and derivative overlays are automatic. Use the optional fields for a phase plot (vertical versus horizontal). Every stage is available as a table.')
        else:
            app.text(p,'fdm_ivp_reference','Optional reference table: x, y per line; no header',5)


def merge_column(text, pairs, label):
    """Preserve existing reference strings and append matched history values."""
    from numerics import number
    lines=[s.strip() for s in text.splitlines() if s.strip()]
    if not lines:raise ValueError('First enter x, reference and their data rows in Compare.')
    def split(s):return [v.strip() for v in (s.split(',') if ',' in s else s.split())]
    header=split(lines[0]);rows=[split(s) for s in lines[1:]]
    if len(header)<2 or not rows or any(len(r)!=len(header) for r in rows):raise ValueError('Enter a rectangular table beginning with x and reference.')
    label=label.strip()
    if not label or ',' in label or '\n' in label or label in header:raise ValueError('Use a distinct column label without commas or line breaks.')
    if len(set(x for x,y in pairs))!=len(pairs):raise ValueError('Selected history column contains repeated x values; choose a solution table.')
    new=[]
    for row in rows:
        x=number(row[0]);matches=[v for at,v in pairs if math.isclose(at,x,rel_tol=1e-10,abs_tol=1e-12)]
        if len(matches)!=1:raise ValueError(f'No unique matching result node at x={x:g}. No interpolation is performed.')
        new.append(', '.join(row+[format(matches[0],'.17g')]))
    return ', '.join(header+[label])+'\n'+'\n'.join(new)


def history_column(app):
    if not app.history:
        messagebox.showinfo('Results','Solve a problem first, or paste a complete comparison table.');return
    dialog=tk.Toplevel(app.root);dialog.title('Add comparison column');dialog.geometry('620x290')
    p=ttk.Frame(dialog,padding=12);p.pack(fill='both',expand=True)
    choice=ttk.Combobox(p,state='readonly',values=[f'{i+1}. {t["title"]}' for i,t in enumerate(app.history)]);choice.pack(fill='x')
    selectors=[]
    for title in ('x column','Numerical value column'):
        ttk.Label(p,text=title).pack(anchor='w');box=ttk.Combobox(p,state='readonly');box.pack(fill='x');selectors.append(box)
    label=tk.StringVar();ttk.Label(p,text='New series label').pack(anchor='w');ttk.Entry(p,textvariable=label).pack(fill='x')
    def columns(event=None):
        cols=app.history[choice.current()]['columns']
        for box in selectors:box.configure(values=[f'{i}: {c}' for i,c in enumerate(cols)])
        selectors[0].current(cols.index('x') if 'x' in cols else 0)
        selectors[1].current(cols.index('y') if 'y' in cols else min(1,len(cols)-1))
    choice.bind('<<ComboboxSelected>>',columns);choice.current(len(app.history)-1);columns()
    def add():
        try:
            data=app.history[choice.current()];ix,iy=[s.current() for s in selectors]
            pairs=[(float(r[ix]),float(r[iy])) for r in data['rows'] if isinstance(r[ix],(int,float)) and not isinstance(r[ix],bool) and isinstance(r[iy],(int,float)) and not isinstance(r[iy],bool)]
            if not pairs:raise ValueError('Select numeric columns.')
            widget=app.texts['comparison_data'];value=merge_column(widget.get('1.0','end-1c'),pairs,label.get())
        except (ValueError,IndexError) as error:messagebox.showerror('Comparison column',str(error));return
        widget.delete('1.0','end');widget.insert('1.0',value);dialog.destroy()
    ttk.Button(p,text='Append matched values',command=add).pack(pady=12)
