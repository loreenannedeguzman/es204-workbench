"""Hidden GUI + JSON-worker tests; preserves the user's running application."""
import importlib.machinery
import json
import time
import tkinter as tk
from pathlib import Path
from unittest.mock import patch
import lesson_ui
from test_systems import system, fdm
from systems_ui import merge_column

module=importlib.machinery.SourceFileLoader('workbench_app',str(Path(__file__).with_name('ES204_App.pyw'))).load_module()
root=tk.Tk();root.withdraw();app=module.App(root)
directory=Path(__file__).with_name('gui_test_output');directory.mkdir(exist_ok=True)

def run(index,prefix,params,preview=False):
    app.tabs.select(index)
    for key,value in params.items():
        name=prefix+'_'+key
        if name in app.fields:app.fields[name].set(str(value))
        elif name in app.texts:app.texts[name].delete('1.0','end');app.texts[name].insert('1.0',str(value))
    errors=[]
    with patch.object(module.messagebox,'showerror',side_effect=lambda *args:errors.append(args)):
        app.run(equations_only=preview)
        deadline=time.monotonic()+30
        while app.process is not None and time.monotonic()<deadline:root.update();time.sleep(.02)
        assert not errors,errors
        assert app.process is None and app.result
    return app.result

try:
    assert [app.tabs.tab(t,'text') for t in app.tabs.tabs()]==['ODE','Data / analysis','Roots','BVP','Matrix','Systems','FDM IVP','Compare']
    result=run(5,'systems',system(order='Second order',rhs='-u\n-v',initial='1,0,0,1',plot_x='v',plot_y='u'))
    assert result['tables'][0]['columns']==['i','x','u','du','v','dv']
    assert len(app.datasets)==4
    for i in range(len(result['charts'])):app.chart_choice.current(i);app.show_chart();app.figure.savefig(directory/f'systems_{i}.png')
    text=app.report_text();assert 'Problem-specific K_4' in text and 'Stage evaluations' in text
    result=run(6,'fdm_ivp',fdm(reference='0,1\n0.5,1.64872\n1,2.71828'))
    assert any('reference comparison' in t['title'] for t in result['tables'])
    values=result['datasets'][0]
    merged=merge_column('x, reference\n0,1\n0.5,1.64872\n1,2.71828',list(zip(values['x'],values['y'])),'numerical')
    result=run(7,'comparison',dict(data=merged))
    assert result['tables'][-1]['title']=='Accuracy summary at supplied nodes'
    for i in range(len(result['charts'])):app.chart_choice.current(i);app.show_chart();app.figure.savefig(directory/f'comparison_{i}.png')
    before=app.inputs();history=list(app.history)
    saved=directory/'session.json'
    with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=str(saved)):lesson_ui.save_session(app)
    app.fields['systems_names'].set('changed')
    with patch.object(lesson_ui.filedialog,'askopenfilename',return_value=str(saved)):lesson_ui.load_session(app)
    assert app.inputs()==before and app.history==history and app.tabs.index('current')==7
    for index,prefix,params in [(5,'systems',system()),(6,'fdm_ivp',fdm()),(7,'comparison',dict(data=merged))]:
        out=run(index,prefix,params,True);assert out['equations'] and not out['datasets']
        path=directory/(prefix+'.tex')
        with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=str(path)):lesson_ui.export_equations(app)
        assert r'\begin{document}' in path.read_text()
    # Old sessions do not contain newly introduced fields and keep their tab indices.
    old=json.loads(saved.read_text());old['tab']=3
    old['inputs']={k:v for k,v in old['inputs'].items() if not k.startswith(('systems_','fdm_ivp_','comparison_'))}
    saved.write_text(json.dumps(old))
    with patch.object(lesson_ui.filedialog,'askopenfilename',return_value=str(saved)):lesson_ui.load_session(app)
    assert app.tabs.index('current')==3 and app.fields['systems_method'].get()=='RK4'
    print('New GUI: systems/FDM/comparison workers, plots, equation export, full session round trip and old-session compatibility passed.')
finally:app.close()
