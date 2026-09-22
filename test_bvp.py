import unittest
import json
import numpy as np
from bvp_methods import bvp,matrix


def base(**kw):
    p=dict(a='1',b='0',c='0',d='0',left='0',right='1',grid='Points n (lesson)',size='5',
           lalpha='1',lbeta='0',lgamma='0',ralpha='1',rbeta='0',rgamma='1',
           method='Finite differences',solver='Direct')
    p.update(kw);return p


def get(out,title):return next(t for t in out['tables'] if t['title']==title)['rows']


class BVPTests(unittest.TestCase):
    def test_dirichlet_line(self):
        out=bvp(base());data=out['datasets'][0]
        np.testing.assert_allclose(data['y'],data['x'],atol=1e-12)
        json.dumps(out,allow_nan=False)

    def test_lecture_example3(self):
        out=bvp(base(a='5+x',b='-2',c='1',d='8*x',right='1.5',size='7',lgamma='1',rgamma='10'))
        table=get(out,'Linear system — augmented matrix')
        # The lesson writes the matrix divided by h^2.
        np.testing.assert_allclose(np.array(table)[:,1:-1]/.25**2,
            [[-167,80,0,0,0],[92,-175,84,0,0],[0,96,-183,88,0],[0,0,100,-191,92],[0,0,0,104,-199]])
        self.assertAlmostEqual(out['datasets'][0]['y'][-1],10)

    def test_lecture_robin_example4(self):
        out=bvp(base(c='-2',d='cos(x)',right='pi/2',grid='Step h',size='pi/8',lgamma='3',ralpha='1',rbeta='1',rgamma='0'))
        ys=out['datasets'][0]['y']
        # Slides solve with rounded five-decimal matrix entries; retain full precision.
        np.testing.assert_allclose(ys,[3,1.60184,.84021,.44676,.25013],atol=2e-5)
        self.assertLess(get(out,'Discrete system checks')[0][1],1e-12)
        M=np.array(get(out,'Linear system — augmented matrix'))[:,1:]
        np.testing.assert_allclose(M[:,:-1],[[-2.30843,1,0],[1,-2.30843,1],[0,.73583,-1.25174]],atol=1e-5)

    def test_boundary_types(self):
        # y=x^2; derivative conditions and Robin on both sides.
        cases=[dict(lalpha='0',lbeta='1',lgamma='0',rgamma='1'),
               dict(lalpha='1',lbeta='1',lgamma='0',ralpha='2',rbeta='1',rgamma='4'),
               dict(lalpha='1',lbeta='0',lgamma='0',ralpha='0',rbeta='1',rgamma='2')]
        for extra in cases:
            out=bvp(base(d='2',**extra));data=out['datasets'][0]
            np.testing.assert_allclose(data['y'],np.array(data['x'])**2,atol=1e-12)

    def test_neumann_unique(self):
        # y=x^2, y''-y=2-x^2 with derivatives at both ends.
        out=bvp(base(c='-1',d='2-x^2',lalpha='0',lbeta='1',lgamma='0',ralpha='0',rbeta='1',rgamma='2'))
        np.testing.assert_allclose(out['datasets'][0]['y'],np.array(out['datasets'][0]['x'])**2,atol=1e-12)

    def test_singular_neumann_rejected(self):
        for method in ('Finite differences','Shooting'):
            with self.assertRaises(ValueError):bvp(base(method=method,lalpha='0',lbeta='1',lgamma='0',ralpha='0',rbeta='1',rgamma='0'))

    def test_shooting_example1_converges_with_refinement(self):
        errors=[]
        for n in (5,9,17):
            out=bvp(base(method='Shooting',b='-2',c='1',d='-8+7*x-x^2',right='2',size=str(n),rgamma='2'))
            xs=np.array(out['datasets'][0]['x']);ys=out['datasets'][0]['y']
            errors.append(max(abs(ys-(3*xs-xs**2))))
        self.assertLess(errors[2],errors[1]/8)
        self.assertLess(errors[1],errors[0]/8)

    def test_shooting_mixed_pc(self):
        out=bvp(base(method='Shooting',b='-2/x',d='-2',left='1',right='3',size='81',lalpha='0',lbeta='1',lgamma='-1',rgamma='3',ivp='Predictor-Corrector'))
        xs=np.array(out['datasets'][0]['x'])
        self.assertLess(max(abs(np.array(out['datasets'][0]['y'])-(-xs**3+xs**2+21))),.01)

    def test_lecture_neumann_shooting(self):
        out=bvp(base(method='Shooting',c='1',d='1',right='pi/2',grid='Step h',size='pi/10',lalpha='0',lbeta='1',lgamma='-3',ralpha='0',rbeta='1',rgamma='-2',guess0='2',guess1='4'))
        trajectory=get(out,'Shooting trajectory 0')
        self.assertAlmostEqual(trajectory[-1][2],-1.999778,places=5)
        self.assertAlmostEqual(trajectory[-1][3],-1.00034,places=5)
        self.assertLess(abs(get(out,'Shooting trials')[-1][4]),1e-8)

    def test_lecture_matrix(self):
        p=dict(matrix='-2.30843 1 0\n1 -2.30843 1\n0 .73583 -1.25174',rhs='-2.85753 .10905 .05901',tolerance='1e-10')
        for method in ('Jacobi','Gauss-Seidel','SOR'):
            out=matrix(dict(p,solver=method,omega='1.15'))
            actual=[r[1] for r in get(out,'Matrix solution')]
            expected=np.linalg.solve([[-2.30843,1,0],[1,-2.30843,1],[0,.73583,-1.25174]],[-2.85753,.10905,.05901])
            np.testing.assert_allclose(actual,expected,atol=1e-9)
            np.testing.assert_allclose(actual,[1.60184,.84021,.44676],atol=2e-5)

    def test_grid_modes(self):
        for grid,size in [('Points n (lesson)','5'),('Intervals','4'),('Step h','.25')]:
            self.assertEqual(len(bvp(base(grid=grid,size=size))['datasets'][0]['x']),5)
        with self.assertRaises(ValueError):bvp(base(grid='Step h',size='.3'))

    def test_equations_only(self):
        for method in ('Finite differences','Shooting'):
            out=bvp(base(method=method,equations_only=True,solver='Jacobi'))
            self.assertFalse(out['datasets']);self.assertTrue(out['equations'])
            self.assertFalse(any(t['title'].endswith(' iterations') for t in out['tables']))

    def test_matrix_methods(self):
        p=dict(matrix='4 -1 0\n-1 4 -1\n0 -1 3',rhs='15 10 10',tolerance='1e-10')
        for method in ('Direct','Jacobi','Gauss-Seidel','SOR'):
            out=matrix(dict(p,solver=method,omega='1.2'))
            np.testing.assert_allclose([r[1] for r in get(out,'Matrix solution')],[5,5,5],atol=1e-9)
            json.dumps(out,allow_nan=False)

    def test_matrix_initial_and_one_iteration(self):
        p=dict(matrix='4 1\n2 3',rhs='1 2',initial='0 0',limit='1',tolerance='1e-20')
        j=matrix(dict(p,solver='Jacobi'));g=matrix(dict(p,solver='Gauss-Seidel'))
        self.assertAlmostEqual(get(j,'Matrix solution')[1][1],2/3)
        self.assertAlmostEqual(get(g,'Matrix solution')[1][1],.5)
        self.assertTrue(j['warnings'])

    def test_matrix_stopping_comparison(self):
        p=dict(matrix='1',rhs='1',solver='Jacobi',criterion='Absolute update',tolerance='1',limit='1')
        self.assertFalse(get(matrix(dict(p,comparison='<')),'Matrix checks')[0][0])
        self.assertTrue(get(matrix(dict(p,comparison='<=')),'Matrix checks')[0][0])

    def test_invalid_matrix(self):
        for p in [dict(matrix='1 1\n1 1',rhs='1 2'),dict(matrix='0 1\n1 2',rhs='1 2',solver='Jacobi'),dict(matrix='1',rhs='1',solver='SOR',omega='2')]:
            with self.assertRaises(ValueError):matrix(p)

if __name__=='__main__':unittest.main()
