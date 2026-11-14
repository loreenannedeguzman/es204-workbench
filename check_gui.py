"""Hidden-window integration checks, including worker and full-session round trip."""
import importlib.machinery
import tkinter as tk
import time
import tempfile
from pathlib import Path
from unittest.mock import patch
from bvp_methods import bvp
from test_bvp import base
import lesson_ui

module=importlib.machinery.SourceFileLoader('workbench_app',str(Path(__file__).with_name('ES204_App.pyw'))).load_module()
root=tk.Tk();root.withdraw();app=module.App(root)
try:
    assert len(app.tabs.tabs())==10
    app.tabs.select(3)
    for key,val in base().items():
        if 'bvp_'+key in app.fields:app.fields['bvp_'+key].set(val)
    app.run()
    deadline=time.monotonic()+20
    while app.process is not None and time.monotonic()<deadline:root.update();time.sleep(.02)
    assert app.process is None and app.result and app.result['datasets']
    assert len(app.datasets)==1
    app.use_dataset()
    assert '1' in app.texts['data'].get('1.0','end')
    with tempfile.TemporaryDirectory() as td:
        path=str(Path(td)/'session.json')
        before=app.inputs();history=list(app.history)
        with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=path):lesson_ui.save_session(app)
        app.fields['bvp_size'].set('99');app.history=[]
        with patch.object(lesson_ui.filedialog,'askopenfilename',return_value=path):lesson_ui.load_session(app)
        assert app.inputs()==before and app.history==history
        with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=str(Path(td)/'equations.tex')):lesson_ui.export_equations(app)
        assert '\\begin{document}' in (Path(td)/'equations.tex').read_text()
    app.run(equations_only=True)
    deadline=time.monotonic()+20
    while app.process is not None and time.monotonic()<deadline:root.update();time.sleep(.02)
    assert app.result['equations'] and not app.result['datasets']
    assert app.datasets # previously calculated solution retained for post-processing
    print('GUI: ten tabs, actual worker, equation preview, data transfer, session round trip, LaTeX export OK.')
finally:app.close()
