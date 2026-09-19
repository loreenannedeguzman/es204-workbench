"""Problem Set 2 regression checks: python validate_problem_set.py."""
import ast
import contextlib
import io
import math
from pathlib import Path
from unittest.mock import patch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

path=Path(__file__).with_name('numerical_methods_program.py')
text=path.read_text(encoding='utf-8')
ns={'__name__':'validation'}
exec(compile(text,str(path),'exec'),ns)
output=[]
checks=[]

def check(name,condition):
    assert condition,name
    checks.append(name)

def run(name,answers=()):
    stream=io.StringIO();answers=iter(map(str,answers))
    with patch('builtins.input',lambda prompt='':next(answers)),contextlib.redirect_stdout(stream),patch.object(plt,'show',lambda:None):
        result=ns[name]()
    output.append(stream.getvalue())
    return result,stream.getvalue()

# Each dependent block must explain its prerequisite without running it.
for name in ['ps2_1b','ps2_1c','ps2_1d','ps2_2b','ps2_2c']:
    result,out=run(name)
    check(name+' prerequisite',result is None and 'first' in out)

p1,_=run('ps2_1a',[''])
check('RK4 grids have 6 and 11 nodes',len(p1['runs'][0]['x'])==6 and len(p1['runs'][1]['x'])==11)
check('RK4 endpoints',all(math.isclose(v['x'][-1],2) for v in p1['runs'].values()))

# Independent high-precision RK4 arithmetic for the two tutorial meshes.
import mpmath as mp
mp.mp.dps=40
for m in [0,1]:
    h=mp.mpf('0.4')/2**m;x=mp.mpf(0);y=mp.mpf(1);ys=[y]
    def rhs(x,y):return x-y*y-mp.exp(-x)
    for i in range(5*2**m):
        a=h*rhs(x,y)
        b=h*rhs(x+h/2,y+a/2)
        c=h*rhs(x+h/2,y+b/2)
        d=h*rhs(x+h,y+c)
        y+=(a+2*b+2*c+d)/6;x+=h;ys.append(y)
    check('High-precision RK4 m='+str(m),np.allclose(p1['runs'][m]['y'],list(map(float,ys)),rtol=1e-12,atol=1e-12))

common,_=run('ps2_1b',['2'])
check('Shared nodes use no interpolation',all(r['coarse_source']==r['fine_source']=='node' for r in common['rows']))
requested,out=run('ps2_1b',['1'])
check('Requested grid explicit assumption','ASSUMPTION' in out and any(r['coarse_source']=='interpolated' for r in requested['rows']))
check('Requested x values',[r['x'] for r in requested['rows']]==[0,.5,1,1.5,2])
for row in common['rows']:
    expected=100*abs((row['fine']-row['coarse'])/row['fine'])
    check('Approximate error at '+str(row['x']),math.isclose(row['error'],expected,abs_tol=1e-14))

fit,_=run('ps2_1c')
check('Problem 1 polynomial degree 10',sp.degree(fit['polynomial'])==10)
check('Polynomial evaluation reproduces RK4 nodes',np.allclose([fit['function'](x) for x in fit['x']],fit['y'],rtol=1e-9,atol=1e-9))
minimum,_=run('ps2_1d',['',''])
check('Secant verifies local minimum',minimum['minimum'] and minimum['residual']<1e-5 and minimum['curvature']>0)
near=[fit['function'](minimum['x']+d) for d in [-.001,.001]]
check('Stationary point is locally lower',all(y>minimum['y'] for y in near))

p2,_=run('ps2_2a')
check('PC h .15 yields 17 points',len(p2['y'])==17 and math.isclose(p2['h'],.15))
# Independent trapezoidal corrector calculation using 40-digit arithmetic.
x=mp.mpf(0);y=mp.mpf(0);h=mp.mpf('.15');ys=[y]
def rhs2(x,y):return mp.sqrt(x)*mp.sin(2*x)-2*y
for i in range(16):
    slope=rhs2(x,y)
    predicted=y+h*slope
    y+=h*(slope+rhs2(x+h,predicted))/2
    x+=h;ys.append(y)
check('High-precision one-pass PC matches',np.allclose(p2['y'],list(map(float,ys)),rtol=1e-12,atol=1e-12))
check('Zero exact reference marked undefined',math.isnan(p2['errors'][0]))
for i in range(1,9):
    check('Tabulated exact error '+str(i),math.isclose(p2['errors'][i],100*abs((p2['exact_y'][i]-p2['y'][2*i])/p2['exact_y'][i]),rel_tol=1e-12))
derivatives,_=run('ps2_2b')
integrals,_=run('ps2_2c')
check('Both derivative arrays have 17 values',len(derivatives['first'])==len(derivatives['second'])==17)
check('Trapezoidal independent interval sum',math.isclose(integrals['trapezoidal'],sum(.15*(p2['y'][i]+p2['y'][i+1])/2 for i in range(16)),abs_tol=1e-14))
check('Simpson independent panel sum',math.isclose(integrals['simpson'],sum(.15*(p2['y'][i]+4*p2['y'][i+1]+p2['y'][i+2])/3 for i in range(0,16,2)),abs_tol=1e-14))

least,_=run('ps2_2d')
check('Least-squares degree six',sp.degree(least['polynomial'])==6)
# A scaled polynomial fit is an independent, well-conditioned numerical reference.
reference=np.polynomial.Polynomial.fit(least['x'],least['y'],6)
check('Least-squares agrees with scaled fit',np.allclose([least['function'](x) for x in least['x']],reference(least['x']),rtol=1e-10,atol=1e-10))
check('Termwise integral agrees with scaled fit',math.isclose(least['integral'],reference.integ()(2.4)-reference.integ()(0),abs_tol=1e-11))
check('Antiderivative differentiation exact',sp.expand(sp.diff(least['antiderivative'],ns['x_symbol'])-least['polynomial'])==0)

# Cubic FD data: central first differences have known h**2 truncation error;
# all the second derivative stencils here must be exact on a cubic.
saved=ns['problem_set_results']['2a']
xs=np.linspace(0,2.4,17); h=.15
ns['problem_set_results']['2a']={'x':list(xs),'y':list(xs**3),'h':h}
cubic,_=run('ps2_2b')
check('FD second derivatives exact on cubic',np.allclose(cubic['second'],6*xs,atol=1e-11))
check('Central first derivative truncation term',np.allclose(cubic['first'][1:-1],3*xs[1:-1]**2+h*h,atol=1e-11))
check('Endpoint first derivative truncation term',np.allclose([cubic['first'][0],cubic['first'][-1]],[3*xs[0]**2-2*h*h,3*xs[-1]**2-2*h*h],atol=1e-11))
ns['problem_set_results']['2a']=saved
ns['problem_set_results']['2b']=derivatives

# Plotting must not call a numerical method, differentiate, or solve a system.
def forbidden(*a,**kw):raise AssertionError('Plot attempted to recalculate')
with patch.dict(ns,{'rk4':forbidden,'predictor_corrector':forbidden}),patch.object(sp,'diff',forbidden),patch.object(np.linalg,'solve',forbidden):
    run('ps2_plots',['all'])
check('Seven requested figures rendered from saved data',len(plt.get_fignums())==7)
figs=[plt.figure(i) for i in plt.get_fignums()]
check('Derivative graph has both curves',any(f.axes[0].get_title().startswith('Problem 2b') and len(f.axes[0].lines)==2 for f in figs))
folder=Path(__file__).parent/'tmp'/'problem_set_plots';folder.mkdir(parents=True,exist_ok=True)
for i,fig in enumerate(figs,1):fig.savefig(folder/f'plot_{i}.png',dpi=120,bbox_inches='tight')
plt.close('all')

# The alternative h0=.5 is only used when explicitly entered by the user.
run('ps2_1a',['.5'])
check('New RK4 grids invalidate dependents',all(k not in ns['problem_set_results'] for k in ['1b','1c','1d']))
alternative,_=run('ps2_1b',['1'])
check('Clarified h0 .5 requested points need no interpolation',all(r['coarse_source']==r['fine_source']=='node' for r in alternative['rows']))
run('ps2_1a',[''])
run('ps2_1c')
run('ps2_1d',['',''])

tree=ast.parse(text)
menu=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='menu')
choices=[int(n.test.comparators[0].value) for n in ast.walk(menu) if isinstance(n,ast.If) and isinstance(n.test,ast.Compare) and isinstance(n.test.left,ast.Name) and n.test.left.id=='choice']
check('One handler per option 0-48',sorted(choices)==list(range(49)))
_,out=run('menu',['','44','y','45','y','46','y','47','y','48','2b','n'])
check('Problem 2 actual menu sequence','Integral of fitted f(x)' in out)
plt.close('all')
print('\n'.join('PASS '+x for x in checks))
print(f'\n{len(checks)} problem-set checks passed.')
print('Problem 1 minimum:',minimum)
print('Problem 2 integrals:',integrals)
print('Problem 2 least-squares integral:',least['integral'])
