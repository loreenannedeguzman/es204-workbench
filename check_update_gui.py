from pathlib import Path
import runpy
import time
import tkinter as tk

module=runpy.run_path(str(Path(__file__).with_name('ES204_App.pyw')),run_name='gui_test')
module['enable_dpi']()
root=tk.Tk();root.withdraw();app=module['App'](root)
errors=[];warnings=[]
module['messagebox'].showerror=lambda title,text:errors.append(text)
module['messagebox'].showwarning=lambda title,text:warnings.append(text)
def wait():
    deadline=time.monotonic()+30
    while app.process is not None and time.monotonic()<deadline:
        root.update();time.sleep(.01)
    assert app.process is None,'Timeout'
    assert not errors,errors
try:
    values=dict(formula='20*cos(x)-40*y',x0='0',xf='1',y0='0',
        grid_mode='Powers of ten n=10^p',first_p='1',last_p='6',
        exact='20/1601*(40*cos(x)+sin(x)-40*exp(-40*x))',
        stop_metric='Final true error %',stop_tolerance='1e-5',comparison='<')
    for k,v in values.items():app.fields[k].set(v)
    for k,v in app.method_vars.items():v.set(k in ['Euler','Predictor-Corrector','RK3'])
    assert app.fields['output_mode'].get()=='Final values only'
    app.run();wait()
    assert [r[1] for r in app.result['tables'][0]['rows'] if r[7]]==[1000000,10000,1000]
    assert 'signed true error %' in app.selected_csv()
    assert not warnings
    app.fields['exact'].set('');app.fields['last_p'].set('3')
    app.texts['reference'].insert('1.0','1,0.280494')
    app.run();wait()
    assert warnings and 'unchanged' in warnings[0]
    assert not any(r[7] for r in app.result['tables'][0]['rows'])
    app.fields['grid_mode'].set('Step sizes h');app.fields['output_mode'].set('Full tables')
    app.fields['steps'].set('.05');app.fields['stop_metric'].set('No early stop')
    app.run();wait()
    assert app.datasets and len(app.datasets[0]['x'])==21
    print('PASS: GUI million-step run, separate stopping, CSV, precision warning and full-table compatibility.')
finally:app.close()
