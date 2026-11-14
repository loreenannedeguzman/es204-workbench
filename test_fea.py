import unittest
import numpy as np
import sympy as sp
from fea_methods import calculate, expression, x, integrate_entry


def base(**kw):
    p=dict(form='Self-adjoint',p='1',r='0',f='2',nodes='0,.5,1',integration='Auto',left_alpha='1',left_beta='0',left_gamma='0',right_alpha='1',right_beta='0',right_gamma='0',workflow='Assemble only')
    p.update(kw)
    return p


def solve(p):
    # Test harness independently partitions the exposed assembly.
    result=calculate(p);s=result['fea_system'];k=np.array(s['matrix']);f=np.array(s['rhs']);u=np.zeros(len(f))
    for i,v in s['prescribed'].items():u[int(i)-1]=v
    order=[i-1 for i in s['free_nodes']]
    reduced=k[np.ix_(order,order)];rhs=(f-k@u)[order]
    return calculate(dict(p,workflow='Check partition and solve',free_order=','.join(str(i+1) for i in order),manual_matrix='\n'.join(' '.join(str(v) for v in row) for row in reduced),manual_rhs=' '.join(str(v) for v in rhs)))


class FEA(unittest.TestCase):
    def test_quadratic(self):
        r=solve(base());np.testing.assert_allclose(r['fea_solution'],[0,.25,0],atol=1e-14)
        self.assertEqual(len(r['datasets']),1)

    def test_no_automatic_partition(self):
        r=calculate(base());self.assertNotIn('fea_solution',r)
        np.testing.assert_allclose(r['fea_system']['matrix'],[[2,-2,0],[-2,4,-2],[0,-2,2]])
        np.testing.assert_allclose(r['fea_system']['rhs'],[.5,1,.5])

    def test_nonzero_essential(self):
        np.testing.assert_allclose(solve(base(f='0',left_gamma='1',right_gamma='3'))['fea_solution'],[1,2,3])

    def test_robin_both_ends(self):
        # y=1+x, -y''=0. Left y+y'=2; right y+2y'=4.
        r=solve(base(f='0',left_alpha='1',left_beta='1',left_gamma='2',right_alpha='1',right_beta='2',right_gamma='4'))
        np.testing.assert_allclose(r['fea_solution'],[1,1.5,2])

    def test_neumann(self):
        r=solve(base(f='0',left_gamma='1',right_alpha='0',right_beta='1',right_gamma='2'))
        np.testing.assert_allclose(r['fea_solution'],[1,2,3])

    def test_standard_transformation(self):
        # x*y''+2*y'=2, solution y=x on [1,2]. p=x^2, f=-2*x.
        p=base(form='Standard',a='x',b='2',c='0',f='2',nodes='1,1.5,2',left_gamma='1',right_gamma='2')
        r=solve(p);np.testing.assert_allclose(r['fea_solution'],[1,1.5,2])
        coeff=next(e['latex'] for e in r['equations'] if e['title']=='Region 1: coefficients')
        self.assertIn('x^{2}',coeff)

    def test_piecewise_flux(self):
        r=solve(base(f='0',regions='0; .5; 1; 0; 0\n.5; 1; 2; 0; 0',right_gamma='1'))
        np.testing.assert_allclose(r['fea_solution'],[0,2/3,1])

    def test_region_interface_required(self):
        with self.assertRaisesRegex(ValueError,'interface'):calculate(base(loads='0;.3;1\n.3;1;2'))

    def test_gauss_cubic(self):
        out=dict(tables=[],equations=[],notes=[])
        v=integrate_entry(out,x**3,0,2,'Gauss 2','test')
        self.assertAlmostEqual(v,4)
        self.assertEqual(len(out['tables'][0]['rows']),2)

    def test_nonpolynomial_gauss(self):
        r=calculate(base(r='exp(sqrt(x))',integration='Gauss 2'))
        pts=.25+.25*np.array([-1,1])/np.sqrt(3)
        expected=2+.25*np.sum(np.exp(np.sqrt(pts))*(1-pts/.5)**2)
        self.assertAlmostEqual(r['fea_system']['matrix'][0][0],expected)

    def test_wrong_partition_rejected(self):
        with self.assertRaisesRegex(ValueError,'does not match'):calculate(base(workflow='Check partition and solve',free_order='2',manual_matrix='3',manual_rhs='1'))

    def test_user_rounding_used(self):
        r=calculate(base(workflow='Check partition and solve',free_order='2',manual_matrix='4',manual_rhs='1.0000005'))
        self.assertAlmostEqual(r['fea_solution'][1],1.0000005/4,places=14)

    def test_point_load(self):
        np.testing.assert_allclose(solve(base(f='0',point_loads='2,1'))['fea_solution'],[0,.25,0])

    def test_singular_neumann(self):
        with self.assertRaisesRegex(ValueError,'singular'):solve(base(f='0',left_alpha='0',left_beta='1',right_alpha='0',right_beta='1'))

    def test_singular_coefficient(self):
        with self.assertRaises(ValueError):calculate(base(p='x-.25'))
        with self.assertRaises(ValueError):calculate(base(p='1/(x-.25)'))
        with self.assertRaises(ValueError):calculate(base(form='Standard',a='x',b='0',c='0',f='0'))

    def test_expression_security(self):
        for s in ['__import__("os")','x.__class__','[1][0]','open("file")']:
            with self.assertRaises(ValueError):expression(s)
        self.assertEqual(expression('ln(x)+x^2'),sp.log(x)+x*x)

    def test_reordered_partition(self):
        r=calculate(base(f='0',left_gamma='1',right_alpha='0',right_beta='1',right_gamma='1',workflow='Check partition and solve',free_order='3,2',manual_matrix='2,-2\n-2,4',manual_rhs='1,2'))
        np.testing.assert_allclose(r['fea_solution'],[1,1.5,2])


if __name__=='__main__':unittest.main()
