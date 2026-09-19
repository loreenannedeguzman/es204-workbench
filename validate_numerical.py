"""Focused lesson regression checks; run with python validate_numerical.py."""
import ast
import builtins
import contextlib
import io
import math
from pathlib import Path
import sys
from unittest.mock import patch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

path=Path(__file__).with_name('numerical_methods_program.py')
text=path.read_text(encoding='utf-8')
ns={'__name__':'study_validation'}
exec(compile(text,str(path),'exec'),ns)
checks=[]

def check(label, condition):
    assert condition,label
    checks.append(label)

def close(a,b,tol=1e-8):
    return math.isclose(a,b,rel_tol=tol,abs_tol=tol)

def call(name, answers=(), *args, **kwargs):
    answers=iter(map(str,answers)); captured={}; stream=io.StringIO()
    target=ns[name].__code__
    def trace(frame,event,arg):
        if event=='return' and frame.f_code is target:
            captured.update(frame.f_locals)
    with patch('builtins.input',lambda prompt='':next(answers)), patch.object(plt,'show',lambda:None), contextlib.redirect_stdout(stream):
        sys.setprofile(trace)
        try: result=ns[name](*args,**kwargs)
        finally: sys.setprofile(None)
    plt.close('all')
    return result,captured,stream.getvalue()

ode=lambda x,y:4*math.exp(.8*x)-.5*y
for name,expected,h in [('euler',3.702440882143687,.1),('predictor_corrector',3.7535337829889754,.1),('rk2_1',3.7535337829889754,.1),('rk4',3.751532463952617,.25)]:
    value=ns[name](ode,0,2,.5,h)[1][-1]
    check('A-1 Example 4: '+name,close(value,expected))
check('A-1 Example 2 Euler unstable h=.2',ns['euler'](lambda x,y:-10*y,0,1,1,.2)[1]==[1,-1,1,-1,1,-1])
check('A-1 Example 2 Euler h=.05',close(ns['euler'](lambda x,y:-10*y,0,1,1,.05)[1][-1],2**-20))
check('A-1 Example 3 PC h=.05',close(ns['predictor_corrector'](lambda x,y:-10*y,0,1,1,.05)[1][-1],.625**20))
for name in ['euler','predictor_corrector','rk2_1','rk2_2','rk3','rk4']:
    for h in [0,-.1,.3,3]:
        try: ns[name](ode,0,2,1,h)
        except ValueError: pass
        else: raise AssertionError((name,h,'invalid mesh accepted'))
    check(name+' rejects invalid meshes',True)

# A nonautonomous nonlinear one-step check distinguishes Ralston from midpoint.
f=lambda x,y:x*x+y*y
h=.1;x=.2;y=.4
k1=f(x,y); k2=f(x+3*h/4,y+3*h*k1/4)
expected=y+h*(k1/3+2*k2/3)
check('A-1 RK2.2 Ralston coefficients',close(ns['rk2_2'](f,x,y,x+h,h)[1][-1],expected,1e-13))
# All six solvers' convergence orders on y'=y.
for name,order in [('euler',1),('predictor_corrector',2),('rk2_1',2),('rk2_2',2),('rk3',3),('rk4',4)]:
    errors=[abs(ns[name](lambda x,y:y,0,1,1,h)[1][-1]-math.e) for h in [.1,.05]]
    check(name+' convergence order',abs(math.log2(errors[0]/errors[1])-order)<.12)

saved={k:ns[k] for k in ['f','exact','x0','xf','y0','h0','exact_at_xf','HAS_EXACT_FUNCTION']}
ns.update(f=lambda x,y:y/x+x*math.sin(x),exact=lambda x:(math.sqrt(2)*x-2*x*math.cos(x))/2,x0=math.pi/4,xf=3*math.pi/2,y0=0,h0=5*math.pi/24,HAS_EXACT_FUNCTION=True)
_,loc,out=call('euler_m_refinement')
check('A-2 Example 2 Euler m=7',loc['answer'][0]==7)
result,loc,out=call('lagrange_hmax',(),ns['euler'],ns['f'],ns['exact'],ns['x0'],ns['y0'],ns['xf'])
check('A-2 Lagrange h_max near .005648',close(result[1],.005648,2e-6))
check('A-2 h_max verified integer mesh',result[2]<=1)
_,loc,out=call('rk4_m_refinement')
check('A-2 Example 2 RK4 m=1',loc['answer'][0]==1)
result,loc,out=call('rk4_successive_m_error',['0','.5'],ns['f'],ns['x0'],ns['y0'],ns['xf'],ns['h0'])
check('A-2 successive RK4 agrees',result[0]==1)
ns.update(saved)

# Exercise the independent finite difference blocks using the lesson's printed data.
xs=[0,.5,1,1.5,2,2.5,3]
ys=[-1,-.569661,-.103466,.5811,1.59478,3.00478,4.85186]
xtext=','.join(map(str,xs));ytext=','.join(map(str,ys));exact='x^2-2*x-3*exp(-x)+2'
_,loc,_=call('finite_difference_approximations',['x^3',1,.1])
check('Function first differences',close(loc['central'],3.01))
_,loc,_=call('finite_difference_tabulated',[xtext,ytext,exact])
check('A-2 Example 4 forward first derivative',close(loc['first_derivatives'][0],.860678))
check('A-2 Example 4 backward first derivative',close(loc['first_derivatives'][-1],3.69416))
check('A-2 Example 4 forward second derivative',close(loc['second_derivatives'][0],.143424))
_,loc,_=call('central_second_derivative_tabulated',[xtext,ytext,exact])
check('A-2 Example 5 central second derivative',close(loc['second_derivatives'][0],.143424))
for endpoint in ['1','2']:
    _,loc,_=call('higher_order_endpoint_differences',[xtext,ytext,endpoint,exact])
    i=0 if endpoint=='1' else -1
    if i==0:
        d1=(-3*ys[0]+4*ys[1]-ys[2]);d2=(2*ys[0]-5*ys[1]+4*ys[2]-ys[3])/.25
    else:
        d1=3*ys[-1]-4*ys[-2]+ys[-3];d2=(2*ys[-1]-5*ys[-2]+4*ys[-3]-ys[-4])/.25
    check('A-2 higher-order endpoint '+endpoint,close(loc['higher_first'],d1) and close(loc['higher_second'],d2))
for name,expected in [('trapezoidal_rule',.25*(ys[0]+ys[-1]+2*sum(ys[1:-1]))),('simpson_one_third_rule',(ys[0]+ys[-1]+4*sum(ys[1:-1:2])+2*sum(ys[2:-1:2]))/6)]:
    _,loc,_=call(name,['2',xtext,ytext,''])
    check('A-2 Example 6 '+name,close(loc['integral'],expected))
    _,loc,_=call(name,['1','x',-1,1,2,'0'])
    check(name+' zero integral error does not crash',math.isnan(loc['error_percent']))
_,loc,_=call('compare_integration_methods',['exp(x)*ln(x^2+1)/cos(x)',0,1,7,6,''])
check('A-2 Exercise 2 trapezoidal',close(loc['trapezoidal_result'],.83613,1e-5))
check('A-2 Exercise 2 Simpson',close(loc['simpson_result'],.81376,1e-5))

points=[str(len(xs))]+[str(v) for pair in zip(xs,ys) for v in pair]
_,vand,_=call('vandermonde_approximation',points+['n'])
_,lag,_=call('lagrangian_polynomial',points+['n'])
vf=vand['polynomial_function']
lf=ns['sp'].lambdify(ns['x_symbol'],lag['polynomial'],'math')
check('A-3 Vandermonde interpolates all nodes',all(close(vf(x),y) for x,y in zip(xs,ys)))
check('A-3 Lagrange and Vandermonde agree',all(close(vf(x),lf(x),1e-7) for x in np.linspace(0,3,31)))
_,loc,_=call('lagrangian_finite_difference',points+['n'])
check('A-3 Lagrange/FD first derivative endpoints',close(loc['fd_first'][0],(-3*ys[0]+4*ys[1]-ys[2])))
data=[4,3,4.234,5,6.980,7,11.509,9,18.975]
_,loc,_=call('least_squares',data+['5','0'])
check('A-3 Example 4 exponential lowest loss',loc['best']['name']=='Exponential')
check('A-3 Example 4 exponential parameters',close(loc['best']['a'],1.999,.001) and close(loc['best']['b'],.25,.001))

_,loc,_=call('bisection_method',['exp(x)-x^2-10',2.5,3,.01,100,''])
check('A-3 Example 6 bisection',close(loc['root'],2.9189453125))
poly='-1+.997393*x-.488922*x^2+.480639*x^3-.106544*x^4+.0149587*x^5-.000990578*x^6'
_,loc,_=call('newton_raphson_method',[poly,1.5,.0001,100,''])
check('A-3 Example 8 Newton',close(loc['root'],1.090277,2e-6))
_,loc,_=call('secant_method',['x^2-2',1,2,.00001,100,''])
check('Secant residual criterion',abs(loc['f'](loc['root']))<.00001)
_,loc,_=call('secant_method',['x^2+1',-1,1,.00001,100])
check('Secant equal function values handled',loc['root'] is None)
_,loc,_=call('newton_raphson_method',['x^2',0,.0001,100,''])
check('Newton exact zero root handled',loc['root']==0)

# A-3 Example 10: independently calculate each fixed-point corrector pass.
for iterations in range(1,5):
    ode=lambda x,y:-x*x*math.cos(y)
    xx=0; yy=2
    for step in range(8):
        guess=yy+.25*ode(xx,yy)
        for k in range(iterations): guess=yy+.125*(ode(xx,yy)+ode(xx+.25,guess))
        yy=guess;xx+=.25
    actual=ns['predictor_corrector'](ode,0,2,2,.25,iterations)[1][-1]
    check('A-3 Example 10 corrector passes '+str(iterations),close(actual,yy,1e-13))

# Menu dispatch smoke test: every displayed entry calls a defined callable.
tree=ast.parse(text)
menu=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='menu')
choices=[n for n in ast.walk(menu) if isinstance(n,ast.If) and isinstance(n.test,ast.Compare) and isinstance(n.test.left,ast.Name) and n.test.left.id=='choice']
check('Menu options 0 through 48 have unique handlers',sorted(int(n.test.comparators[0].value) for n in choices)==list(range(49)))
for node in ast.walk(menu):
    if isinstance(node,ast.Call) and isinstance(node.func,ast.Name):
        check('Menu callable '+node.func.id,callable(ns.get(node.func.id,getattr(builtins,node.func.id,None))))

# Plotting uses saved results, never re-enters input or reruns a solver.
for name in ['plot_vandermonde','plot_lagrangian','plot_finite_difference','plot_least_squares','plot_root_finding','plot_root_convergence']:
    call(name)
    check(name+' uses saved results without prompts',True)
saved_cache=ns['saved_results'].copy()
ns['saved_results'].clear()
for name in ['plot_vandermonde','plot_lagrangian','plot_finite_difference','plot_least_squares','plot_root_finding','plot_root_convergence']:
    _,_,out=call(name)
    check(name+' explains missing results', 'first' in out)
ns['saved_results'].update(saved_cache)
call('plot_error_comparison',['RK2.1,RK4'],ns['f'],ns['exact'],ns['x0'],ns['y0'],ns['xf'],ns['h0'],ns['methods'])
check('All repaired plot routes render with Agg',True)
call('menu',['','0'])
check('Menu startup and exit',True)
# Additional boundary and cache-only checks.
_,loc,_=call('trapezoidal_rule',['2','0,1,3','0,1,3',''])
check('Unequal-spacing trapezoidal rule integrates a line',close(loc['integral'],4.5))
_,_,out=call('simpson_one_third_rule',['2','0','1'])
check('Simpson rejects a one-point table','at least three' in out)
_,_,out=call('finite_difference_approximations',['x',0,0])
check('Finite difference rejects zero h','positive' in out)
_,loc,_=call('bisection_method',['x',0,1])
check('Bisection caches exact endpoint root',ns['saved_results']['root']['root']==0)
_,loc,_=call('bisection_method',['x^2-2',1,2,1e-15,1])
check('Unconverged bisection cache is labelled',not ns['saved_results']['root']['converged'])
with patch.object(np.linalg,'solve',side_effect=AssertionError('Plot recomputed a linear solve')):
    for name in ['plot_vandermonde','plot_lagrangian','plot_finite_difference','plot_least_squares','plot_root_finding','plot_root_convergence']:
        call(name)
check('Cached plots perform no linear solves',True)
_,_,out=call('menu',['','27']+points+['n','y','34','0'])
check('Menu calculation then cached plotting then exit','Exiting program' in out)
Path('validation_results.txt').write_text('\n'.join('PASS '+x for x in checks)+'\n',encoding='utf-8')
print(f'{len(checks)} checks passed (including menu callable checks).')
