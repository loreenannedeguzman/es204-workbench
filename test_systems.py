import json
import math
import unittest
import numpy as np
from systems_methods import systems, fdm_ivp, comparison
from systems_ui import merge_column
from numerics import execute, METHODS
from bvp_methods import bvp
from test_bvp import base


def system(**kw):
    p=dict(order='First order',names='u,v',rhs='v\n-u',initial='1,0',left='0',right='1',size='0.05',grid='Step h',method='RK4')
    p.update(kw);return p


def fdm(**kw):
    p=dict(a='1',b='0',c='-1',d='0',y0='1',v0='1',left='0',right='1',size='0.05',grid='Step h')
    p.update(kw);return p


class SystemsTests(unittest.TestCase):
    def test_coupled_rotation_and_rk4_order(self):
        errors=[]
        for h in ('0.1','0.05'):
            out=systems(system(size=h,plot_x='v',plot_y='u'))
            last=out['tables'][0]['rows'][-1]
            errors.append(max(abs(last[2]-math.cos(1)),abs(last[3]+math.sin(1))))
            self.assertEqual(len(out['tables'][1]['rows']),4*(len(out['tables'][0]['rows'])-1))
            self.assertEqual(out['charts'][-1]['xlabel'],'v')
            self.assertEqual(out['charts'][-1]['curves'][0]['y'],out['datasets'][0]['y'])
            json.dumps(out,allow_nan=False)
        self.assertGreater(errors[0]/errors[1],15)

    def test_second_order_coupled_polynomials(self):
        out=systems(system(order='Second order',rhs='v+2-x^3\nu+6*x-x^2',initial='0,0,0,0',size='0.01'))
        data=np.array(out['tables'][0]['rows']);x=data[:,1]
        np.testing.assert_allclose(data[:,2:].T,[x*x,2*x,x**3,3*x*x],atol=2e-9)
        self.assertEqual(out['tables'][0]['columns'],['i','x','u','du','v','dv'])
        self.assertEqual(len(out['charts'][1]['curves']),2)

    def test_every_vector_method_and_pc_passes(self):
        for method in METHODS:
            out=systems(system(rhs='1\n2',initial='0,1',method=method,passes='3'))
            np.testing.assert_allclose(out['tables'][0]['rows'][-1][2:],[1,3],atol=1e-12)
            self.assertTrue(out['equations'])
        one=systems(system(method='Predictor-Corrector',size='0.5',passes='1'))
        two=systems(system(method='Predictor-Corrector',size='0.5',passes='2'))
        # Independent trapezoidal fixed-point updates for the first step.
        np.testing.assert_allclose(one['tables'][0]['rows'][1][2:],[.875,-.5])
        np.testing.assert_allclose(two['tables'][0]['rows'][1][2:],[.875,-.46875])

    def test_preview_and_invalid_inputs(self):
        out=execute({'action':'systems','inputs':system(equations_only=True)})
        self.assertTrue(out['equations']);self.assertFalse(out['tables'])
        for kw in (dict(names='u,u'),dict(names='sin,v'),dict(rhs='v\nw'),dict(initial='1'),dict(size='0.3'),dict(rhs='1/(x-x)\nv'),dict(order='Second order',names='u,du',initial='0,0,0,0'),dict(plot_x='bad',plot_y='u')):
            with self.subTest(kw=kw), self.assertRaises((ValueError,TypeError,KeyError,ZeroDivisionError)):
                systems(system(**kw))


class FDMTests(unittest.TestCase):
    def test_polynomial_variable_coefficients(self):
        # y=x²+2x+3, a=1+x, b=x, c=2; exact for both second-order starts.
        for start in ('Central ghost / Taylor (second order)','Forward derivative (second order)'):
            out=fdm_ivp(fdm(a='1+x',b='x',c='2',d='8+8*x+4*x^2',y0='3',v0='2',start=start))
            xs=np.array(out['datasets'][0]['x'])
            np.testing.assert_allclose(out['datasets'][0]['y'],xs**2+2*xs+3,atol=1e-11)
            json.dumps(out,allow_nan=False)

    def test_accuracy_orders(self):
        for start,low,high in [('Central ghost / Taylor (second order)',3.7,4.3),('Forward derivative (second order)',3.7,4.3),('Forward derivative (first order)',1.9,2.2)]:
            errors=[]
            for h in ('.05','.025'):
                out=fdm_ivp(fdm(size=h,start=start));errors.append(abs(out['datasets'][0]['y'][-1]-math.e))
            self.assertTrue(low<errors[0]/errors[1]<high,(start,errors))

    def test_matches_independent_full_matrix(self):
        # Solve all central equations plus the second-order forward initial derivative.
        n=20;h=1/n;A=np.zeros((n+1,n+1));rhs=np.zeros(n+1)
        A[0,0]=1;rhs[0]=1
        A[1,:3]=[-3,4,-1];rhs[1]=2*h
        for i in range(1,n):A[i+1,i-1:i+2]=[1,-2-h*h,1]
        expected=np.linalg.solve(A,rhs)
        out=fdm_ivp(fdm(start='Forward derivative (second order)'))
        np.testing.assert_allclose(out['datasets'][0]['y'],expected,atol=1e-12)

    def test_preview_and_singular_recurrence(self):
        out=fdm_ivp(fdm(equations_only=True));self.assertTrue(out['equations']);self.assertFalse(out['datasets'])
        for kw in (dict(a='0'),dict(b='-40'),dict(size='0.3'),dict(start='unknown')):
            with self.assertRaises(ValueError):fdm_ivp(fdm(**kw))

    def test_supplied_reference_and_off_grid(self):
        out=fdm_ivp(fdm(c='0',d='0',y0='0',v0='1',reference='0,0\n0.05,0.05\n0.123,1'))
        rows=next(t['rows'] for t in out['tables'] if 'reference comparison' in t['title'])
        self.assertIsNone(rows[0][5]);self.assertAlmostEqual(rows[1][3],0)
        self.assertIsNone(rows[2][1]);json.dumps(out,allow_nan=False)


class CompareTests(unittest.TestCase):
    def test_signed_errors_summary_and_zero(self):
        out=comparison(dict(data='x, reference, first, second\n0,0,1,-1\n1,2,3,1'))
        first=out['tables'][1]['rows']
        self.assertEqual(first[1][3:],[-1,1,-50]);self.assertIsNone(first[0][-1])
        self.assertEqual(out['tables'][-1]['rows'][0],['first',2,1,1,50])
        self.assertEqual(len(out['charts'][-1]['curves'][0]['x']),1)
        json.dumps(out,allow_nan=False)

    def test_history_merge_preserves_reference_and_matches_nodes(self):
        original='x, reference\n0.1, 0.123450\n0.2, 0.567890'
        merged=merge_column(original,[(.10000000000000002,4),(.2,5)],'new')
        self.assertIn('0.123450',merged);self.assertIn('0.567890',merged)
        out=comparison(dict(data=merged));self.assertEqual(out['tables'][0]['rows'][0][2],4)
        with self.assertRaises(ValueError):merge_column(original,[(.1,1)],'new')

    def test_bvp_reference_and_combined_derivative_plot(self):
        out=bvp(base(reference='0,0\n0.5,0.5\n1,1'))
        chart=next(c for c in out['charts'] if c['title']=='Solution and first derivative')
        self.assertEqual(len(chart['curves']),2)
        errors=next(t for t in out['tables'] if t['title'].endswith('reference comparison'))
        self.assertEqual(errors['rows'][-1][3],0)


if __name__=='__main__':unittest.main()
