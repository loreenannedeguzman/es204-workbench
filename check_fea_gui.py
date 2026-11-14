"""Hidden GUI integration; all generated files stay beside this test."""
from pathlib import Path
import importlib.machinery
import json
import time
import tkinter as tk
from unittest.mock import patch
import lesson_ui
from weak_ui import show, equation_image
from test_fea import base

module=importlib.machinery.SourceFileLoader('fea_app',str(Path(__file__).with_name('ES204_App.pyw'))).load_module()
root=tk.Tk();root.withdraw();app=module.App(root)
output=Path(__file__).with_name('fea_test_output');output.mkdir(exist_ok=True)

def run(params):
    app.tabs.select(9)
    for key,value in params.items():
        name='fea_'+key
        if name in app.fields:app.fields[name].set(value)
        elif name in app.texts:app.texts[name].delete('1.0','end');app.texts[name].insert('1.0',value)
    errors=[]
    with patch.object(module.messagebox,'showerror',side_effect=lambda *args:errors.append(args)):
        app.run();deadline=time.monotonic()+60
        while app.process is not None and time.monotonic()<deadline:root.update();time.sleep(.02)
    assert not errors,errors
    assert app.process is None and app.result.get('fea')

try:
    assert len(app.tabs.tabs())==10
    run(base())
    assert 'fea_solution' not in app.result
    run(dict(workflow='Check partition and solve',free_order='2',manual_matrix='4',manual_rhs='1'))
    assert app.result['fea_solution']==[0,.25,0]
    assert app.datasets[0]['label']=='FEA nodal solution'
    result=app.result
    for page in range((len(result['equations'])+23)//24):
        app.math_page=page;show(app,reset=False)
        assert 'Rendering unavailable' not in app.math_text.get('1.0','end')
        assert len(app.math_images)==len(result['equations'][page*24:(page+1)*24])
    for i in range(len(result['charts'])):
        app.chart_choice.current(i);app.show_chart();app.figure.savefig(output/f'chart_{i}.png')
    before=app.inputs();session=output/'session.json';tex=output/'equations.tex'
    with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=str(session)):lesson_ui.save_session(app)
    app.fields['fea_nodes'].set('bad')
    with patch.object(lesson_ui.filedialog,'askopenfilename',return_value=str(session)):lesson_ui.load_session(app)
    assert app.inputs()==before and app.tabs.index('current')==9
    with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=str(tex)):lesson_ui.export_equations(app)
    assert 'Manual partition' in tex.read_text() and r'\begin{bmatrix}' in tex.read_text()
    # Validate compatibility with older saved sessions (no FEA input keys).
    old=json.loads(session.read_text());old['inputs']={k:v for k,v in old['inputs'].items() if not k.startswith('fea_')};old['tab']=0;old['result']=None
    old_path=output/'legacy_session.json';old_path.write_text(json.dumps(old))
    with patch.object(lesson_ui.filedialog,'askopenfilename',return_value=str(old_path)):lesson_ui.load_session(app)
    assert app.fields['fea_nodes'].get()==''
    # Exact rendering pipeline used by Tk, saved for visual inspection.
    from PIL import Image,ImageDraw
    chosen=[result['equations'][i] for i in [0,2,3,5,6,10,12,23,-3,-1] if i<len(result['equations'])]
    images=[equation_image(eq) for eq in chosen]
    sheet=Image.new('RGB',(max(im.width for im in images)+40,sum(im.height+45 for im in images)+20),'white')
    draw=ImageDraw.Draw(sheet);top=10
    for eq,im in zip(chosen,images):
        draw.text((20,top),eq['title'],fill='black');sheet.paste(im,(20,top+25),im);top+=im.height+45
    sheet.save(output/'rendered_working.png')
    print('PASS: FEA worker, manual solve, every equation page rendered, plots, LaTeX export, session round trip, legacy session.')
finally:app.close()
