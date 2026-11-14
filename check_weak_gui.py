"""Hidden-window worker, visual-math, export and session checks."""
from pathlib import Path
import importlib.machinery
import tkinter as tk
from unittest.mock import patch
import time
import lesson_ui
from test_weak_forms import base

module=importlib.machinery.SourceFileLoader('weak_app',str(Path(__file__).with_name('ES204_App.pyw'))).load_module()
root=tk.Tk();root.withdraw();app=module.App(root)
output=Path(__file__).with_name('weak_test_output');output.mkdir(exist_ok=True)
try:
    assert len(app.tabs.tabs())==10
    app.tabs.select(8)
    for key,value in base(a='sin(x)',b='-cos(x)',left='1',right='2').items():app.fields['weak_'+key].set(value)
    errors=[]
    with patch.object(module.messagebox,'showerror',side_effect=lambda *args:errors.append(args)):
        app.run();deadline=time.monotonic()+25
        while app.process is not None and time.monotonic()<deadline:root.update();time.sleep(.02)
        assert not errors,errors
        assert app.process is None and 'weak_form' in app.result
    assert len(app.math_images)==len(app.result['equations'])
    assert 'Rendering unavailable' not in app.math_text.get('1.0','end')
    assert app.output_tabs.select()==str(app.math_frame)
    before=app.inputs();saved=output/'session.json';tex=output/'equations.tex'
    with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=str(saved)):lesson_ui.save_session(app)
    app.fields['weak_a'].set('bad')
    with patch.object(lesson_ui.filedialog,'askopenfilename',return_value=str(saved)):lesson_ui.load_session(app)
    assert app.inputs()==before and len(app.math_images)==len(app.result['equations'])
    with patch.object(lesson_ui.filedialog,'asksaveasfilename',return_value=str(tex)):lesson_ui.export_equations(app)
    assert 'Transfer one derivative' in tex.read_text() and '\\int_' in tex.read_text()
    # Render the exact same equation images used in the GUI for visual QA.
    from weak_ui import formula_image
    from PIL import Image,ImageDraw
    eqs=app.result['equations'];imgs=[formula_image(e['latex']) for e in eqs]
    sheet=Image.new('RGB',(max(im.width for im in imgs)+40,sum(im.height+55 for im in imgs)+20),'white')
    draw=ImageDraw.Draw(sheet);top=10
    for e,im in zip(eqs,imgs):
        draw.text((20,top),e['title'],fill='black');sheet.paste(im,(20,top+25),im);top+=im.height+55
    sheet.save(output/'rendered_equations.png')
    print('PASS: ten tabs, actual worker, all equations visually rendered, session restoration, explanatory LaTeX export.')
finally:app.close()
