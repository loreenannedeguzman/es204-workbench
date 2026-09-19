import unittest
import copy
from matplotlib.figure import Figure
from plotting import draw_chart

class PlotTests(unittest.TestCase):
    def test_signed_and_absolute(self):
        chart={'title':'Errors','ylabel':'error','error_plot':True,'xscale':'log',
               'tolerance':1e-5,'curves':[{'label':'Euler','x':[10,100,1000000],
               'y':[1e7,-.01,-1.3e-6],'style':'nodes'}]}
        original=copy.deepcopy(chart)
        ax=draw_chart(Figure().subplots(),chart)
        self.assertEqual(ax.get_yscale(),'symlog')
        self.assertEqual(list(ax.lines[0].get_ydata()),chart['curves'][0]['y'])
        ax=draw_chart(Figure().subplots(),chart,'Absolute')
        self.assertEqual(ax.get_yscale(),'log')
        self.assertTrue(all(v>0 for v in ax.lines[0].get_ydata()))
        self.assertEqual(chart,original)

    def test_log_warning(self):
        chart={'title':'Data','ylabel':'y','curves':[{'label':'Data','x':[1,2],'y':[-1,1]}]}
        ax=draw_chart(Figure().subplots(),chart,yscale='Log')
        self.assertIn('omits',ax.texts[0].get_text())
