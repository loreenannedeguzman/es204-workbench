"""Exercise the real GUI/worker path and capture only its own window."""
from pathlib import Path
import runpy
import tkinter as tk
import time
from PIL import ImageGrab

module=runpy.run_path(str(Path(__file__).with_name('ES204_App.pyw')),run_name='test_app')
module['enable_dpi']()
root=tk.Tk();root.withdraw();app=module['App'](root)
errors=[]
module['App'].run.__globals__['messagebox'].showerror=lambda title,text:errors.append((title,text))

def wait():
    deadline=time.monotonic()+35
    while app.process is not None and time.monotonic()<deadline:
        root.update();time.sleep(.03)
    assert app.process is None,'Calculation timeout'
    assert not errors,errors

try:
    for k,v in dict(formula='-2*y+sqrt(x)*sin(2*x)',x0='0',xf='2.4',y0='0',steps='.15',iterations='1').items():app.fields[k].set(v)
    for name,var in app.method_vars.items():var.set(name=='Predictor-Corrector')
    app.texts['reference'].insert('1.0','0,0\n0.3,0.0323187\n0.6,0.140207\n0.9,0.270998\n1.2,0.342188\n1.5,0.290355\n1.8,0.104387\n2.1,-0.164493\n2.4,-0.419418')
    app.run();wait()
    assert abs(app.datasets[0]['y'][-1]+.417282787678)<1e-11
    app.use_dataset();app.tabs.select(1);root.update()
    app.fields['operation'].set('Finite differences');app.run();wait()
    assert len(app.result['charts'][0]['curves'])==2
    app.output_tabs.select(app.plot_frame)
    root.geometry('1120x760');root.deiconify();root.update()
    app.canvas.draw();root.update()
    widget=app.canvas.get_tk_widget()
    assert abs(app.figure.bbox.width-widget.winfo_width())<2
    assert abs(app.figure.bbox.height-widget.winfo_height())<2
    bounds=(root.winfo_rootx(),root.winfo_rooty(),root.winfo_rootx()+root.winfo_width(),root.winfo_rooty()+root.winfo_height())
    ImageGrab.grab(bbox=bounds).save(Path(__file__).with_name('gui_preview.png'))
    # Tab changes, table exports and plot rendering must retain the results.
    assert 'Finite-difference derivatives' in app.report_text()
    assert app.selected_csv()
    root.geometry('900x600');root.update()
    assert app.run_button.winfo_viewable()
    assert app.run_button.winfo_rootx()+app.run_button.winfo_width()<=root.winfo_rootx()+root.winfo_width()
    assert app.run_button.winfo_rooty()+app.run_button.winfo_height()<=root.winfo_rooty()+root.winfo_height()
    print('PASS: GUI PC + exact table, saved dataset, finite differences, combined plot; CSV; 900x600 controls.')
finally:
    app.close()
