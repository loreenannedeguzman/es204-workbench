"""Weak-form input and offline, visually typeset mathematics in Tk."""
import io
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties


def formula_image(latex):
    stream=io.BytesIO()
    math_to_image('$'+latex+'$',stream,prop=FontProperties(family='DejaVu Serif',size=16),dpi=125,format='png',color='#142d47')
    stream.seek(0)
    return Image.open(stream).convert('RGBA')


def build(app,ScrollForm):
    form=ScrollForm(app.tabs);app.tabs.add(form,text='Weak form');p=form.body
    app.hint(p,"Derive a scalar linear second-order weak form in 1D. Standard: a*y''+b*y'+c*y=f. Divergence: -(a*y')'+b*y'+c*y=f; a is the flux coefficient p. No grid is needed.")
    app.field(p,'weak_form','Equation form','Standard',choices=['Standard','Divergence'])
    for key,label,default in [('a','Leading / flux coefficient a(x)',''),('b',"Coefficient b(x) of y'",'0'),('c','Coefficient c(x) of y','0'),('f','Right side f(x)',''),('left','Left endpoint',''),('right','Right endpoint','')]:
        app.field(p,'weak_'+key,label,default)
    app.hint(p,"At each endpoint enter alpha*y + beta*y' = gamma. Derivatives are in increasing x. Dirichlet: beta=0. Neumann: alpha=0. Robin: both nonzero. For an outward derivative, include its sign in beta.")
    for side in ('left','right'):
        for key,default in [('alpha','1'),('beta','0'),('gamma','0')]:app.field(p,'weak_'+side+'_'+key,side.title()+' '+key,default)
    ttk.Button(p,text='Derive and show weak form',command=app.run).pack(fill='x',pady=5)
    app.hint(p,'Assumes regular smooth coefficients. Piecewise/interface, nonlinear, higher-order and multidimensional weak forms are outside this workflow. No built-in example values are supplied.')


def build_view(app):
    frame=ttk.Frame(app.output_tabs);app.output_tabs.add(frame,text='Rendered equations')
    frame.rowconfigure(1,weight=1);frame.columnconfigure(0,weight=1)
    bar=ttk.Frame(frame);bar.grid(row=0,column=0,columnspan=2,sticky='ew')
    app.math_page=0
    app.math_page_label=ttk.Label(bar,text='')
    app.math_page_label.pack(side='left',padx=8)
    def move(delta):
        app.math_page=max(0,app.math_page+delta);show(app,reset=False)
    ttk.Button(bar,text='Previous equations',command=lambda:move(-1)).pack(side='left')
    ttk.Button(bar,text='Next equations',command=lambda:move(1)).pack(side='left')
    text=tk.Text(frame,wrap='none',background='white',font=('Segoe UI',11),padx=15,pady=12)
    vs=ttk.Scrollbar(frame,orient='vertical',command=text.yview);hs=ttk.Scrollbar(frame,orient='horizontal',command=text.xview)
    text.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
    text.grid(row=1,column=0,sticky='nsew');vs.grid(row=1,column=1,sticky='ns');hs.grid(row=2,column=0,sticky='ew')
    text.tag_configure('heading',font=('Segoe UI',12,'bold'),foreground='#155e75',spacing1=12,spacing3=5)
    text.tag_configure('body',font=('Segoe UI',10),spacing3=8)
    app.math_frame=frame;app.math_text=text;app.math_images=[]
    text.insert('end','Derive a weak form to see step-by-step rendered equations here.');text.configure(state='disabled')


def show(app,reset=True):
    import textwrap
    if not hasattr(app,'math_text'):return
    t=app.math_text;t.configure(state='normal');t.delete('1.0','end');app.math_images=[]
    result=app.result or {}
    if 'weak_form' not in result and 'fea' not in result:
        t.insert('end','Use the Weak form tab to derive and display a weak formulation.');t.configure(state='disabled');return
    if reset:app.math_page=0
    equations=result['equations']
    page_size=24 if result.get('fea') else max(1,len(equations))
    app.math_page=min(app.math_page,max(0,(len(equations)-1)//page_size))
    start=app.math_page*page_size
    app.math_page_label.configure(text=f'Equations {start+1}–{min(start+page_size,len(equations))} of {len(equations)}')
    for eq in equations[start:start+page_size]:
        t.insert('end',eq['title']+'\n','heading')
        t.insert('end',textwrap.fill(eq['text'],75)+'\n','body')
        try:
            photo=ImageTk.PhotoImage(equation_image(eq),master=app.root)
            app.math_images.append(photo);t.image_create('end',image=photo);t.insert('end','\n\n')
        except Exception as error:
            t.insert('end','Rendering unavailable for this equation: '+str(error)+'\n'+eq['latex']+'\n\n')
    for note in result.get('notes',[]):
        if note.startswith(('Scope:','WARNING:')):t.insert('end',textwrap.fill(note,75)+'\n','body')
    t.configure(state='disabled');t.yview_moveto(0);t.xview_moveto(0)
    app.output_tabs.select(app.math_frame)


def equation_image(eq):
    if 'matrix' not in eq:
        return formula_image(eq['latex'])
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    a=eq['matrix']; rows=len(a); cols=len(a[0])
    fig=Figure(figsize=(max(3,cols*1.05),max(1,rows*.32+.7)),dpi=120)
    FigureCanvasAgg(fig)
    ax=fig.add_subplot(); ax.axis('off')
    ax.set_title('$'+eq['matrix_label']+'$',fontsize=14)
    tbl=ax.table(cellText=[[f'{v:.7g}' for v in row] for row in a],loc='center',cellLoc='center')
    tbl.auto_set_font_size(False);tbl.set_fontsize(10);tbl.scale(1,1.25)
    stream=io.BytesIO();fig.savefig(stream,format='png',bbox_inches='tight',facecolor='white');stream.seek(0)
    return Image.open(stream).convert('RGBA')
