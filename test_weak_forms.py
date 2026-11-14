import json
import unittest
import sympy as sp
from weak_forms import derive
from weak_ui import formula_image
from numerics import execute

x=sp.Symbol('x')


def base(**kw):
    p=dict(form='Standard',a='1',b='0',c='0',f='2',left='0',right='1',
           left_alpha='1',left_beta='0',left_gamma='0',right_alpha='1',right_beta='0',right_gamma='1')
    p.update(kw);return p


class WeakFormTests(unittest.TestCase):
    def test_variable_coefficient_product_rule(self):
        out=derive(base(a='sin(x)',b='-cos(x)',left='1',right='2'))
        terms=out['weak_form']['volume_coefficients']
        self.assertEqual(sp.simplify(sp.sympify(terms[1])+2*sp.cos(x)),0)
        self.assertEqual(sp.simplify(sp.sympify(terms[0])+sp.sin(x)),0)
        self.assertEqual(out['weak_form']['boundary_load'],'0')
        self.assertTrue(any('v(2)=0' in e['latex'] for e in out['equations']))
        json.dumps(out,allow_nan=False)

    def test_manufactured_weak_identity_all_boundary_types(self):
        # Substitute a known exact solution into the produced weak identity,
        # integrate independently, and test several admissible variations.
        y=x*x+2*x+1;yp=sp.diff(y,x)
        for form in ('Standard','Divergence'):
            a=1+x;b=2-x;c=x
            forcing=a*sp.diff(y,x,2)+b*yp+c*y if form=='Standard' else -sp.diff(a*yp,x)+b*yp+c*y
            for bc in [((1,0),(1,0)),((1,0),(0,1)),((0,1),(1,0)),((1,2),(3,4)),((0,1),(0,1))]:
                p=base(form=form,a=str(a),b=str(b),c=str(c),f=str(forcing))
                for side,pt,(al,be) in zip(('left','right'),(0,1),bc):
                    p.update({side+'_alpha':str(al),side+'_beta':str(be),side+'_gamma':str((al*y+be*yp).subs(x,pt))})
                out=derive(p);K,B,C=map(sp.sympify,out['weak_form']['volume_coefficients'])
                for trial in (1+x,1+x*x):
                    v=trial*(x if bc[0][1]==0 else 1)*((1-x) if bc[1][1]==0 else 1)
                    residual=sp.integrate(K*yp*sp.diff(v,x)+B*yp*v+C*y*v-forcing*v,(x,0,1))
                    for row in out['tables'][0]['rows']:
                        at=sp.sympify(row[0]);residual+=(sp.sympify(row[2])*y.subs(x,at)-sp.sympify(row[3]))*v.subs(x,at)
                    self.assertEqual(sp.simplify(residual),0,(form,bc))

    def test_divergence_and_expanded_standard_equivalence(self):
        div=derive(base(form='Divergence',a='1+x^2',b='x'))
        standard=derive(base(a='-1-x^2',b='-x'))
        self.assertEqual(div['weak_form']['volume_coefficients'],standard['weak_form']['volume_coefficients'])

    def test_render_every_step(self):
        cases=[base(),base(a='sin(x)',b='-cos(x)',left='1',right='2'),
               base(form='Divergence',a='1+x',b='exp(x)',c='x^2',f='cos(x)',left_alpha='1',left_beta='2',left_gamma='3',right_alpha='0',right_beta='1',right_gamma='4'),
               base(left='pi/2',right='pi',a='exp(x)',b='1/x',f='sqrt(x)')]
        for p in cases:
            out=derive(p)
            for e in out['equations']:
                with self.subTest(e=e['title'],p=p):
                    im=formula_image(e['latex']);self.assertGreater(im.width,10);self.assertGreater(im.height,10)

    def test_validation(self):
        for kw in [dict(left='1',right='0'),dict(a='0'),dict(a='x'),dict(a='1/x'),dict(a='sin(x)',right='4'),dict(a='abs(x)'),dict(b='y'),dict(left_alpha='0',left_beta='0')]:
            with self.subTest(kw=kw),self.assertRaises(ValueError):derive(base(**kw))

    def test_worker_dispatch(self):
        out=execute(dict(action='weak',inputs=base()))
        self.assertGreater(len(out['equations']),10)
        self.assertFalse(out['datasets']);json.dumps(out,allow_nan=False)

    def test_decimal_zero_boundary_coefficients(self):
        out=derive(base(left_alpha='1.0',left_beta='0.0',right_alpha='0.0',right_beta='1.0',right_gamma='2.0'))
        self.assertEqual(out['tables'][0]['rows'][0][1],'Dirichlet')
        self.assertEqual(out['tables'][0]['rows'][1][1],'Neumann')
        self.assertEqual(float(out['tables'][0]['rows'][1][3]),-2)
        for e in out['equations']:formula_image(e['latex'])


if __name__=='__main__':unittest.main()
