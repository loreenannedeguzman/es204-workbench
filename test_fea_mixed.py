"""Regression checks for termwise Auto integration; no file/session writes."""
import unittest
import math
import sympy as sp
from fea_methods import integrate_entry, x


def run(g,mode='Auto'):
    out=dict(equations=[],tables=[],notes=[])
    return integrate_entry(out,g,0,1,mode,'Integral'),out


class MixedIntegration(unittest.TestCase):
    def test_exact_plus_remaining(self):
        g=sp.cos(x)+sp.exp(sp.sqrt(x))
        points=[.5-.5/math.sqrt(3),.5+.5/math.sqrt(3)]
        expected=math.sin(1)+.5*sum(math.exp(math.sqrt(v)) for v in points)
        actual,out=run(g)
        self.assertAlmostEqual(actual,expected,places=13)
        self.assertTrue(any('analytical terms' in e['title'] for e in out['equations']))
        rows=out['tables'][0]['rows']
        for row in rows:self.assertAlmostEqual(row[4],math.exp(math.sqrt(row[2])),places=13)

    def test_product_is_expanded(self):
        g=(1-x)*(sp.cos(x)+sp.exp(sp.sqrt(x)))
        points=[.5-.5/math.sqrt(3),.5+.5/math.sqrt(3)]
        expected=1-math.cos(1)+.5*sum((1-v)*math.exp(math.sqrt(v)) for v in points)
        self.assertAlmostEqual(run(g)[0],expected,places=13)

    def test_explicit_gauss_still_integrates_all(self):
        points=[.5-.5/math.sqrt(3),.5+.5/math.sqrt(3)]
        expected=.5*sum(math.cos(v)+math.exp(math.sqrt(v)) for v in points)
        self.assertAlmostEqual(run(sp.cos(x)+sp.exp(sp.sqrt(x)),'Gauss 2')[0],expected,places=13)

    def test_polynomial_remains_exact(self):
        v,out=run(x**6)
        self.assertAlmostEqual(v,1/7,places=14)
        self.assertEqual(out['tables'],[])

    def test_analytical_only_still_rejects_unsupported(self):
        with self.assertRaises(ValueError):run(sp.cos(x)+sp.exp(sp.sqrt(x)),'Analytical only')


if __name__=='__main__':unittest.main()
