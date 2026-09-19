import unittest
import numpy as np
from numpy.polynomial import Polynomial
from polynomial_working import append_working

class WorkingTests(unittest.TestCase):
    def result(self,p,a,b):
        out={'notes':[],'tables':[],'charts':[]}
        xs=np.linspace(a,b,17)
        append_working(out,p,xs,p(xs),a,b,True)
        return {t['title']:t for t in out['tables']}
    def test_crossing(self):
        t=self.result(Polynomial([0,1]),-1,1)
        self.assertAlmostEqual(t['Total geometric area']['rows'][0][2],1)
        self.assertEqual(len(t['Area — split integrals']['rows']),2)
    def test_tangent(self):
        t=self.result(Polynomial([0,0,1]),-1,1)
        self.assertAlmostEqual(t['Total geometric area']['rows'][0][2],2/3)
        self.assertEqual(len(t['Polynomial roots in area interval']['rows']),1)
    def test_no_roots(self):
        t=self.result(Polynomial([-2]),0,3)
        self.assertAlmostEqual(t['Total geometric area']['rows'][0][2],6)
    def test_problem3(self):
        xs=np.linspace(.5,2,17);y=-11.;ys=[y];h=xs[1]-xs[0]
        f=lambda x,y:np.exp(-x)+10-y/x
        for x in xs[:-1]:
            a=f(x,y);b=f(x+h/2,y+h*a/2);c=f(x+h/2,y+h*b/2);d=f(x+h,y+h*c)
            y+=h*(a+2*b+2*c+d)/6;ys.append(y)
        p=Polynomial.fit(xs,ys,8);out={'notes':[],'tables':[],'charts':[]}
        append_working(out,p,xs,np.array(ys),.5,2,True)
        tables={t['title']:t for t in out['tables']}
        self.assertAlmostEqual(tables['Total geometric area']['rows'][0][2],6.0198423911,places=8)
        self.assertAlmostEqual(tables['Polynomial roots in area interval']['rows'][0][0],1.142217718088,places=8)
