"""BVP/matrix forms; kept separate so students can study the GUI wiring."""
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import json


def build(app, ScrollForm, methods):
    for title,prefix in [('BVP','bvp'),('Matrix','matrix')]:
        form=ScrollForm(app.tabs);app.tabs.add(form,text=title);p=form.body
        app.hint(p,'Show equations previews the entered problem. Solve also includes equations and tables in Report / formulas.')
        ttk.Button(p,text='Show equations',command=lambda:app.run(equations_only=True)).pack(fill='x',pady=4)
        if prefix=='bvp':
            app.hint(p,"Linear ODE: a(x)y'' + b(x)y' + c(x)y = d(x). Enter coefficients, not an equals sign. z denotes y'.")
            for k,label,value in [('a',"a(x), coefficient of y''",'1'),('b',"b(x), coefficient of y'",'0'),('c','c(x), coefficient of y','0'),('d','d(x), right side','0'),('left','Left endpoint','0'),('right','Right endpoint','1')]:app.field(p,prefix+'_'+k,label,value)
            app.field(p,'bvp_grid','Grid definition','Points n (lesson)',choices=['Points n (lesson)','Intervals','Step h'])
            app.field(p,'bvp_size','Grid value (n, intervals, or h)','5')
            app.hint(p,"Boundary: alpha*y + beta*y' = gamma. Dirichlet: beta=0; Neumann: alpha=0; Robin: both nonzero. Use increasing-x derivatives at both ends.")
            for side,name in [('l','Left'),('r','Right')]:
                for k,value in [('alpha','1'),('beta','0'),('gamma','0')]:app.field(p,'bvp_'+side+k,name+' '+k,value)
            app.field(p,'bvp_method','BVP method','Finite differences',choices=['Finite differences','Shooting'])
            app.field(p,'bvp_exact','Exact y(x), optional comparison')
            app.text(p,'bvp_reference','Supplied reference table (optional): x, y per line',5)
            app.hint(p,'Formula and supplied table comparisons are reported separately. Neither replaces the other.')
            app.hint(p,'Shooting settings (used only for Shooting):')
            app.field(p,'bvp_ivp','IVP step method','RK4',choices=methods)
            for k,label,value in [('passes','PC corrector passes','1'),('guess0','First shooting parameter guess','0'),('guess1','Second shooting parameter guess','1'),('shoot_tol','Shooting tolerance (not percent)','1e-8'),('shoot_limit','Maximum secant updates','20')]:app.field(p,'bvp_'+k,label,value)
            app.field(p,'bvp_shoot_stop','Shooting criterion','Boundary residual',choices=['Boundary residual','Absolute parameter change'])
            app.hint(p,"Parameter is y'(left) for a Dirichlet left boundary; otherwise y(left). The other initial value is derived from the left boundary.")
            app.hint(p,'Finite-difference matrix solver settings:')
        else:
            app.text(p,'matrix_matrix','Matrix A: one row per line, spaces or commas',7)
            app.text(p,'matrix_rhs','Right-side vector b: spaces, commas or newlines',3)
        app.field(p,prefix+'_solver','Linear-system solver','Direct',choices=['Direct','Jacobi','Gauss-Seidel','SOR'])
        app.field(p,prefix+'_initial','Initial vector, blank = zeros')
        app.field(p,prefix+'_omega','SOR omega (0 < omega < 2)','1.1')
        app.field(p,prefix+'_criterion','Matrix stopping criterion','Residual',choices=['Residual','Absolute update','Relative update %'])
        app.field(p,prefix+'_tolerance','Matrix tolerance (% only for relative update)','1e-6')
        app.field(p,prefix+'_comparison','Stopping comparison','<',choices=['<','<='])
        app.field(p,prefix+'_limit','Maximum matrix iterations','1000')


def export_equations(app):
    equations=(app.result or {}).get('equations',[])
    if not equations:
        messagebox.showinfo('Equations','Use Show equations or Solve in BVP, Matrix, Systems, FDM IVP or Compare first.');return
    path=filedialog.asksaveasfilename(defaultextension='.tex',filetypes=[('LaTeX equations','*.tex')])
    if not path:return
    def escape(s):return s.replace('_',r'\_').replace('&',r'\&').replace('%',r'\%')
    body=[r'\documentclass{article}',r'\usepackage[margin=20mm]{geometry}',r'\usepackage{amsmath,amssymb,graphicx}',
          r'\newsavebox{\equationbox}',r'\newcommand{\fitequation}[1]{\sbox{\equationbox}{$\displaystyle #1$}\par\noindent\makebox[\textwidth][c]{\ifdim\wd\equationbox>\textwidth\resizebox{\textwidth}{!}{\usebox{\equationbox}}\else\usebox{\equationbox}\fi}\par}',r'\begin{document}']
    for eq in equations:body.extend([r'\paragraph{'+escape(eq['title'])+'}',r'\fitequation{'+eq['latex']+'}'])
    body.append(r'\end{document}')
    Path(path).write_text('\n'.join(body),encoding='utf-8')


def save_session(app):
    if app.process is not None:
        messagebox.showinfo('Session','Wait for the calculation or cancel it before saving.');return
    path=filedialog.asksaveasfilename(defaultextension='.json',filetypes=[('Workbench session','*.json')])
    if path:
        data=dict(session_version=1,inputs=app.inputs(),result=app.result,history=app.history,datasets=app.datasets,tab=app.tabs.index('current'))
        Path(path).write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')


def load_session(app):
    if app.process is not None:return
    path=filedialog.askopenfilename(filetypes=[('Workbench session','*.json')])
    if not path:return
    try:
        data=json.loads(Path(path).read_text(encoding='utf-8'))
        if data.get('session_version')!=1 or not isinstance(data['inputs'],dict):raise ValueError('Not a workbench session file.')
        for key in ('history','datasets'):
            if not isinstance(data[key],list):raise ValueError('Invalid session data.')
        result=data['result']
        if result is not None and not all(isinstance(result.get(k),list) for k in ('tables','charts','notes')):raise ValueError('Invalid result data.')
        tab=int(data.get('tab',0))
        if not 0<=tab<len(app.tabs.tabs()):raise ValueError('Unknown tab.')
    except (OSError,ValueError,KeyError,TypeError) as e:messagebox.showerror('Load session',str(e));return
    for k,v in app.fields.items():v.set(data['inputs'].get(k,app.field_defaults.get(k,'')))
    for k,w in app.texts.items():w.delete('1.0','end');w.insert('1.0',data['inputs'].get(k,''))
    for k,v in app.method_vars.items():v.set(k in data['inputs'].get('methods',[]))
    app.history=[]
    if result is not None:app.render(result)
    else:
        app.result=None;app.tree.delete(*app.tree.get_children());app.table_choice.set('');app.table_choice.configure(values=[])
        app.chart_choice.set('');app.chart_choice.configure(values=[]);app.figure.clear();app.canvas.draw_idle();app.set_report('Session has no calculated results.')
    app.history=data['history'];app.datasets=data['datasets']
    app.dataset_choice.configure(values=[d['label'] for d in app.datasets]);app.dataset_choice.set(app.datasets[0]['label'] if app.datasets else '')
    app.tabs.select(tab);app.status.set('Session restored: inputs, results, table history and datasets.')
