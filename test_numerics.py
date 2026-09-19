import math
import unittest
import numpy as np
from numerics import solve, analyze, roots, METHODS, execute


class NumericalChecks(unittest.TestCase):
    def problem(self,**kw):
        return dict(formula='x-y**2-exp(-x)',x0='0',xf='2',y0='1',steps='.4,.2',iterations='1',methods=['RK4'],exact='',reference='',**kw)

    def test_rk4_known_results(self):
        result=solve(self.problem())
        self.assertAlmostEqual(result['datasets'][0]['y'][-1],1.115384079826,places=11)
        self.assertAlmostEqual(result['datasets'][1]['y'][-1],1.115740577385,places=11)

    def test_all_methods_converge(self):
        p=self.problem();p.update(formula='y',xf='1',steps='.1,.05',methods=METHODS,exact='exp(x)')
        result=solve(p)
        for i in range(0,len(result['datasets']),2):
            old,new=result['datasets'][i:i+2]
            self.assertLess(abs(new['y'][-1]-math.e),abs(old['y'][-1]-math.e))

    def test_exact_table_pc_and_missing(self):
        p=self.problem();p.update(formula='-2*y+sqrt(x)*sin(2*x)',xf='2.4',y0='0',steps='.15',methods=['Predictor-Corrector'],reference='0,0\n0.3,0.0323187\n0.31,0.04\n2.4,-0.419418')
        result=solve(p)
        self.assertAlmostEqual(result['datasets'][0]['y'][-1],-.417282787678,places=11)
        refs=result['tables'][1]['rows']
        self.assertIsNone(refs[0][4])
        self.assertIsNone(refs[2][1])
        self.assertAlmostEqual(refs[-1][4],.5090893385,places=8)

    def data(self,op,**kwargs):
        return dict(data='0,1\n1,2\n2,5\n3,10\n4,17',operation=op,fit='Exact interpolation',degree='2',query='1.5',a='-1',b='.5',tolerance='1e-9',limit='100',**kwargs)

    def test_polynomial_and_minimum(self):
        r=analyze(self.data('Stationary point (secant)'))
        final=r['tables'][-1]['rows'][0]
        self.assertAlmostEqual(final[0],0,places=8)
        self.assertAlmostEqual(final[1],1,places=8)
        self.assertTrue(final[3])
        evals=next(t for t in r['tables'] if t['title']=='Polynomial evaluation')['rows'][0]
        self.assertAlmostEqual(evals[1],3.25,places=10)

    def test_finite_differences(self):
        result=analyze(self.data('Finite differences'))
        for x,y,d1,d2 in next(t for t in result['tables'] if t['title']=='Finite-difference derivatives')['rows']:
            self.assertAlmostEqual(d1,2*x,places=10)
            self.assertAlmostEqual(d2,2,places=10)

    def test_integrals(self):
        r=analyze(self.data("Simpson's 1/3 rule"))
        self.assertAlmostEqual(r['tables'][1]['rows'][0][-1],4+64/3,places=10)
        p=self.data('Trapezoidal rule');p['data']='0,0\n.1,.1\n1,1'
        self.assertAlmostEqual(analyze(p)['tables'][1]['rows'][0][-1],.5)

    def test_root_methods(self):
        for method in ['Secant','Bisection','Newton-Raphson']:
            r=roots(dict(function='x**2-2',method=method,a='1',b='2',tolerance='1e-10',limit='100',plot_min='0',plot_max='2'))
            self.assertAlmostEqual(r['tables'][-1]['rows'][0][0],math.sqrt(2),places=9)
            self.assertTrue(r['tables'][-1]['rows'][0][2])

    def test_bad_inputs(self):
        for key,value in [('steps','.3'),('formula','y_prime+y'),('exact','-.419418')]:
            p=self.problem();p[key]=value
            with self.assertRaises(ValueError):solve(p)
        p=self.data("Simpson's 1/3 rule");p['data']='0,0\n1,1\n3,9'
        with self.assertRaises(ValueError):analyze(p)

    def test_grid_modes_and_sampling(self):
        p=self.problem();p.update(grid_mode='Interval counts n',steps='5,10',requested='.5,1,1.5',sampling='Linear interpolation')
        r=solve(p)
        self.assertEqual(len(r['datasets'][1]['x']),11)
        sampled=next(t for t in r['tables'] if 'requested' in t['title'])
        self.assertTrue(all(row[1] is not None for row in sampled['rows']))
        p.update(grid_mode='Refinement h0 / 2^m',steps='.4',max_m='2')
        self.assertEqual(len(solve(p)['datasets']),3)

    def test_pc_passes_and_n_search(self):
        p=self.problem();p.update(iterations='1,2,4',methods=['Predictor-Corrector'])
        r=execute({'action':'solve','inputs':p})
        self.assertEqual(len(r['datasets']),6)
        p.update(formula='y',xf='1',iterations='1',methods=['Euler'],exact='exp(x)',grid_mode='Search n (final error)',steps='1',max_n='1000',stop_tolerance='2')
        r=execute({'action':'solve','inputs':p})
        n=len(r['datasets'][0]['x'])-1
        self.assertLessEqual(abs((1+1/n)**n-math.e)/math.e*100,2)
        self.assertGreater(abs((1+1/(n-1))**(n-1)-math.e)/math.e*100,2)

    def test_lagrange_vandermonde_and_transformed(self):
        for fit in ['Exact Lagrange','Exact Vandermonde','Least squares']:
            p=self.data('Polynomial fit');p['fit']=fit
            r=analyze(p)
            val=next(t for t in r['tables'] if t['title']=='Polynomial evaluation')['rows'][0][1]
            self.assertAlmostEqual(val,3.25,places=9)
        p=self.data('Linearized least squares');p['data']='\n'.join(f'{x},{2*math.exp(.3*x)}' for x in range(5))
        r=analyze(p)
        self.assertLess(r['tables'][1]['rows'][0][0],1e-20)

    def test_formula_integration_and_derivative_comparison(self):
        p=self.data("Simpson's 1/3 rule");p.update(data_source='Sample formula',data_formula='x**2',data_min='0',data_max='1',data_n='10')
        self.assertAlmostEqual(analyze(p)['tables'][1]['rows'][0][-1],1/3,places=12)
        p=self.data('Finite differences');p['analysis_exact']='x**2+1'
        r=analyze(p)
        self.assertEqual(len(r['charts'][0]['curves']),2)
        for name in ['Forward','Backward','Central only','Forward higher-order','Backward higher-order']:
            p['scheme']=name;r=analyze(p)
            rows=next(t for t in r['tables'] if t['title']=='Finite-difference derivatives')['rows']
            self.assertTrue(any(row[2] is not None for row in rows))

    def test_relative_root_criterion(self):
        for method in ['Secant','Bisection','Newton-Raphson']:
            r=roots(dict(function='x**2-2',method=method,a='1',b='2',tolerance='1e-6',criterion='Relative step change %',limit='100',plot_min='0',plot_max='2'))
            self.assertTrue(r['tables'][-1]['rows'][0][2])
            self.assertLess(r['tables'][0]['rows'][-1][-1],1e-6)


if __name__=='__main__':unittest.main()
