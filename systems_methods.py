"""User-defined ODE systems, initial-value finite differences and comparisons.

All examples belong in tests/documentation, never in the solver's inputs.
"""
import math
import re
import numpy as np
import sympy as sp
from numerics import expr, number, integer, finite, table, curve, report, ode_step, METHODS
from bvp_methods import output, equation


def uniform_grid(p):
    left, right = number(p['left']), number(p['right'])
    if right <= left:
        raise ValueError('Final x must be greater than initial x.')
    if p.get('grid', 'Step h') == 'Step h':
        h = number(p['size'])
        if h <= 0:
            raise ValueError('h must be positive.')
        count = (right-left)/h
        if not math.isclose(count, round(count), rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError('h must divide the interval; no endpoint adjustment is made.')
        n = round(count)
    else:
        n = integer(p['size'])
    if not 2 <= n <= 10000:
        raise ValueError('Use 2–10000 intervals for full system/FDM tables.')
    return np.linspace(left, right, n+1)


def grid_equation(out, xs):
    h = xs[1]-xs[0]
    equation(out, 'Uniform grid', f'x_i = {xs[0]:.12g} + i*{h:.12g}, i=0,...,{len(xs)-1}',
             rf'x_i={xs[0]:.12g}+i({h:.12g}),\quad i=0,\ldots,{len(xs)-1}')


def system_spec(p):
    names = [s.strip() for s in p['names'].replace(',', ' ').split()]
    reserved = {'x', 'e', 'E', 'I', 'pi', 'oo', 'nan', 'True', 'False'}
    if not names or len(names)>6 or len(set(names)) != len(names):
        raise ValueError('Enter 1–6 distinct dependent-variable names.')
    if any(not re.fullmatch('[A-Za-z][A-Za-z0-9_]*', s) or s in reserved or hasattr(sp, s) for s in names):
        raise ValueError('Use simple variable names such as u, v, w; do not use x, constants or function names.')
    second = p.get('order', 'First order') == 'Second order'
    states = [v for s in names for v in (s, 'd'+s)] if second else names
    if len(set(states)) != len(states) or any(s in reserved or hasattr(sp,s) for s in states):
        raise ValueError('Derivative names (d + variable name) must be distinct and unreserved.')
    lines = [s.strip() for s in p['rhs'].splitlines() if s.strip()]
    if len(lines) != len(names):
        raise ValueError('Enter exactly one right-hand-side formula per named variable, in the same order.')
    raw = [expr(s, 'x '+' '.join(states)) for s in lines]
    rhs = [v for s, f in zip(names,raw) for v in (sp.Symbol('d'+s),f)] if second else raw
    initial = [number(s) for s in p['initial'].replace(',', ' ').split()]
    if len(initial) != len(states):
        raise ValueError('Initial values must follow this state order: '+', '.join(states))
    return names, states, rhs, np.array(initial, dtype=float), second


def stepping_equations(out, states, rhs, method, passes):
    x, h = sp.symbols('x_i h')
    u = sp.Matrix([sp.Symbol(s+'_i') for s in states])
    symbols = [sp.Symbol(s) for s in states]
    def field(at, values):
        mapping = dict(zip(symbols,values));mapping[sp.Symbol('x')] = at
        return sp.Matrix([r.subs(mapping, simultaneous=True) for r in rhs])
    def slope(k):
        return sp.Matrix([sp.Symbol(f'k{k}_{s}') for s in states])
    k1,k2,k3,k4 = [slope(k) for k in range(1,5)]
    stages = [('K_1',field(x,u))]
    if method in ('Predictor-Corrector','RK2.1 (Heun)'):
        stages.append(('K_2',field(x+h,u+h*k1)))
        update = u+h*(k1+k2)/2
    elif method == 'RK2.2 (lesson Ralston)':
        stages.append(('K_2',field(x+3*h/4,u+3*h*k1/4)))
        update = u+h*(k1/3+2*k2/3)
    elif method in ('RK3','RK4'):
        stages.append(('K_2',field(x+h/2,u+h*k1/2)))
        if method=='RK3':
            stages.append(('K_3',field(x+h,u-h*k1+2*h*k2)))
            update=u+h*(k1+4*k2+k3)/6
        else:
            stages += [('K_3',field(x+h/2,u+h*k2/2)),('K_4',field(x+h,u+h*k3))]
            update=u+h*(k1+2*k2+2*k3+k4)/6
    else:
        update=u+h*k1
    for label,value in stages:
        equation(out,'Problem-specific '+label, label+' = '+str(value),label+'='+sp.latex(value))
    if method=='Predictor-Corrector':
        equation(out,'Euler predictor', 'U^(0) = U_i + h K_1',r'U^{(0)}=U_i+hK_1')
        previous=sp.Matrix([sp.Symbol(s+'^{(j)}') for s in states])
        corrected=u+h*(k1+field(x+h,previous))/2
        equation(out,'Trapezoidal corrector',f'U^(j+1) = {corrected}; j=0,...,{passes-1}; U_(i+1)=U^({passes})',
                 r'U^{(j+1)}='+sp.latex(corrected)+rf',\quad j=0,\ldots,{passes-1},\quad U_{{i+1}}=U^{{({passes})}}')
    else:
        equation(out,method+' update','U_(i+1) = '+str(update),r'U_{i+1}='+sp.latex(update))


def systems(p):
    out=output();xs=uniform_grid(p);h=xs[1]-xs[0]
    names,states,rhs,u,second=system_spec(p)
    method=p.get('method','RK4');passes=integer(p.get('passes','1'))
    if method not in METHODS or passes>100:
        raise ValueError('Choose a supported method and 1–100 PC passes.')
    if len(xs)*len(states)*(passes+1 if method=='Predictor-Corrector' else 4)>500000:
        raise ValueError('Reduce the grid or corrector passes: stage reports are limited to 500000 scalar slopes.')
    grid_equation(out,xs)
    equation(out,'State order','U = ['+', '.join(states)+']',r'U='+sp.latex(sp.Matrix([sp.Symbol(s) for s in states])))
    for s,r,val in zip(states,rhs,u):
        equation(out,'Differential equation for '+s, s+"' = "+str(r),sp.latex(sp.Symbol(s))+r'^\prime='+sp.latex(r))
        equation(out,'Initial '+s,f'{s}({xs[0]:.12g})={val:.12g}',rf'{sp.latex(sp.Symbol(s))}({xs[0]:.12g})={val:.12g}')
    stepping_equations(out,states,rhs,method,passes)
    out['notes'].append('All components advance simultaneously at every stage. k values are slopes, not h times slopes. '+
                        ('dNAME is the first derivative of NAME.' if second else 'First-order state system.'))
    if p.get('equations_only'):
        out['notes'].append('Equation preview only; no numerical solution computed.');return out
    f=sp.lambdify([sp.Symbol('x')]+[sp.Symbol(s) for s in states],rhs,'math')
    values=[u.copy()];stages=[]
    for i,x in enumerate(xs[:-1]):
        counter=[0]
        def evaluate(at,state):
            result=np.array([finite(v) for v in f(float(at),*map(float,state))])
            counter[0]+=1
            label=('K1' if counter[0]==1 else f'corrector {counter[0]-1} slope') if method=='Predictor-Corrector' else f'K{counter[0]}'
            stages.append([i,float(x),label,float(at),*map(float,state),*map(float,result)])
            return result
        u=np.array([finite(v) for v in ode_step(evaluate,float(x),u,h,method,passes)])
        values.append(u.copy())
    values=np.array(values)
    out['tables'].append(report('System solution — '+method,['i','x']+states,[[i,float(x),*map(float,row)] for i,(x,row) in enumerate(zip(xs,values))]))
    out['tables'].append(report('Stage evaluations — '+method,['step','x_i','stage','stage x']+['stage '+s for s in states]+['slope '+s for s in states],stages))
    selected=[states.index(s) for s in names]
    out['charts'].append({'title':'System solution','ylabel':'State value','curves':[curve(states[j],xs,values[:,j],'nodes') for j in selected]})
    if second:
        out['charts'].append({'title':'First derivatives','ylabel':'First derivative','curves':[curve(states[j+1],xs,values[:,j+1],'nodes') for j in selected]})
    horizontal=p.get('plot_x','').strip();vertical=p.get('plot_y','').strip()
    if horizontal or vertical:
        if horizontal not in ['x']+states or vertical not in states:
            raise ValueError('Versus plot: horizontal must be x or a state name; vertical must be a state name.')
        vx=xs if horizontal=='x' else values[:,states.index(horizontal)]
        out['charts'].append({'title':vertical+' versus '+horizontal,'xlabel':horizontal,'ylabel':vertical,'curves':[curve(vertical,vx,values[:,states.index(vertical)],'nodes')]})
    out['datasets']=[{'label':method+' / '+s,'x':xs.tolist(),'y':values[:,j].tolist()} for j,s in enumerate(states)]
    return out


def add_reference(out, xs, ys, reference, label):
    """Match supplied nodes only, never interpolate or replace a reference."""
    if not str(reference).strip():return
    refs=table(reference);rows=[];mx=[];my=[];ry=[];err=[];px=[];pe=[]
    for x,r in refs:
        ix=int(np.argmin(abs(xs-x)))
        if not math.isclose(float(xs[ix]),float(x),rel_tol=1e-10,abs_tol=1e-12):
            rows.append([float(x),None,float(r),None,None,None,'No matching grid node']);continue
        y=float(ys[ix]);e=float(r-y);percent=e/r*100 if r!=0 else None
        rows.append([float(x),y,float(r),e,abs(e),percent,'Zero reference: percent undefined' if r==0 else ''])
        mx.append(float(x));my.append(y);ry.append(float(r));err.append(e)
        if percent is not None:px.append(float(x));pe.append(float(percent))
    out['tables'].append(report(label+' reference comparison',['x','numerical','supplied reference','reference - numerical','absolute error','signed error %','note'],rows))
    out['charts'].append({'title':label+' and supplied reference','ylabel':'y','curves':[curve(label,xs,ys,'nodes'),curve('Supplied reference',refs[:,0],refs[:,1],'points')]})
    if mx:out['charts'].append({'title':label+' signed error','ylabel':'Reference - numerical','curves':[curve(label,mx,err,'nodes')]})
    if px:out['charts'].append({'title':label+' percentage error','ylabel':'Error (%)','error_plot':True,'curves':[curve(label,px,pe,'nodes')]})
    out['notes'].append('Reference values are used exactly as entered, at matching grid nodes only. Percentage error is undefined for a zero reference. Rounded reference data limit the accuracy of error estimates.')


def fdm_ivp(p):
    out=output();xs=uniform_grid(p);h=xs[1]-xs[0]
    coeff=[expr(p[k]) for k in ('a','b','c','d')]
    functions=[sp.lambdify(sp.Symbol('x'),v,'math') for v in coeff]
    cv=np.array([[finite(f(x)) for f in functions] for x in xs])
    if np.any(cv[:,0]==0):raise ValueError("The coefficient of y'' must be nonzero at every node.")
    y0,v0=number(p['y0']),number(p['v0'])
    scheme=p.get('start','Central ghost / Taylor (second order)')
    if scheme not in ('Central ghost / Taylor (second order)','Forward derivative (first order)','Forward derivative (second order)'):
        raise ValueError('Choose a supported starting formula.')
    grid_equation(out,xs)
    x=sp.Symbol('x');y=sp.Function('y')(x)
    equation(out,'Governing equation',f'({coeff[0]}) y\'\' + ({coeff[1]}) y\' + ({coeff[2]}) y = {coeff[3]}',
             sp.latex(coeff[0]*sp.diff(y,x,2)+coeff[1]*sp.diff(y,x)+coeff[2]*y)+'='+sp.latex(coeff[3]))
    equation(out,'Initial conditions',f'y({xs[0]:.12g})={y0:.12g}; y\'({xs[0]:.12g})={v0:.12g}',rf'y({xs[0]:.12g})={y0:.12g},\quad y^\prime({xs[0]:.12g})={v0:.12g}')
    equation(out,'Central finite differences',"y'_i ≈ (y_(i+1)-y_(i-1))/(2h); y''_i ≈ (y_(i+1)-2y_i+y_(i-1))/h²",
             r"y'_i\approx\frac{y_{i+1}-y_{i-1}}{2h},\quad y''_i\approx\frac{y_{i+1}-2y_i+y_{i-1}}{h^2}")
    equation(out,'Marching recurrence','y_(i+1) = P_i*y_i + Q_i*y_(i-1) + R_i; i=1,...,n-1',
             r'y_{i+1}=\frac{(2a_i-h^2c_i)y_i-(a_i-hb_i/2)y_{i-1}+h^2d_i}{a_i+hb_i/2}=P_i y_i+Q_i y_{i-1}+R_i')
    denominators=cv[:,0]+h*cv[:,1]/2
    if any(abs(denominators[i]) <= 1e-14*max(abs(cv[i,0]),abs(h*cv[i,1]/2)) for i in range(1,len(xs)-1)):
        raise ValueError('A recurrence denominator a_i+h*b_i/2 is zero or nearly zero. Change h or the formulation.')
    P=(2*cv[:,0]-h*h*cv[:,2])[1:-1]/denominators[1:-1]
    Q=-(cv[:,0]-h*cv[:,1]/2)[1:-1]/denominators[1:-1]
    R=h*h*cv[1:-1,3]/denominators[1:-1]
    if scheme=='Forward derivative (first order)':
        first=y0+h*v0
        formula=r'y_1=y_0+h v_0'
        explanation='First-order forward derivative: y_1=y_0+h*v_0. This start generally reduces global accuracy to first order.'
    elif scheme=='Forward derivative (second order)':
        if abs(4-P[0]) <= 1e-14*max(4,abs(P[0])):raise ValueError('The second-order forward starting equations are singular. Choose another start or h.')
        first=((3+Q[0])*y0+2*h*v0+R[0])/(4-P[0])
        formula=r'\frac{-3y_0+4y_1-y_2}{2h}=v_0,\quad y_1=\frac{(3+Q_1)y_0+2hv_0+R_1}{4-P_1}'
        explanation='Use the second-order forward derivative and the recurrence at i=1 together to determine y_1.'
    else:
        aa,bb,cc,dd=cv[0];first=y0+h*v0+h*h*(dd-bb*v0-cc*y0)/(2*aa)
        formula=r'y_{-1}=y_1-2hv_0,\quad y_1=y_0+hv_0+\frac{h^2}{2a_0}(d_0-b_0v_0-c_0y_0)'
        explanation='Eliminate the ghost node using the centered initial derivative and the ODE at x_0; equivalent to a second-order Taylor start.'
    equation(out,'Starting formula',explanation,formula)
    equation(out,'Substituted starting value',f'h={h:.12g}, y_0={y0:.12g}, v_0={v0:.12g}; y_1={first:.12g}',rf'h={h:.12g},\quad y_0={y0:.12g},\quad v_0={v0:.12g},\quad y_1={first:.12g}')
    for i,(pp,qq,rr) in enumerate(zip(P,Q,R),1):
        equation(out,f'Recurrence at i={i}',f'y_{i+1}=({pp:.12g})*y_{i}+({qq:.12g})*y_{i-1}+({rr:.12g})',
                 rf'y_{{{i+1}}}=({pp:.12g})y_{{{i}}}+({qq:.12g})y_{{{i-1}}}+({rr:.12g})')
    out['notes'].append('Both initial conditions are imposed at the left endpoint. No right-end value is imposed. '+explanation)
    out['notes'].append('Interior stencils are second order; stability and accumulated error depend on the ODE, grid and starting formula. This method is not an RK or shooting calculation.')
    if p.get('equations_only'):return out
    ys=[y0,finite(first)];steps=[]
    for i,(pp,qq,rr) in enumerate(zip(P,Q,R),1):
        new=finite(pp*ys[i]+qq*ys[i-1]+rr)
        steps.append([i,float(xs[i]),float(pp),float(qq),float(rr),ys[i-1],ys[i],new]);ys.append(new)
    ys=np.array(ys);ds=np.gradient(ys,h,edge_order=2);ds[0]=v0
    out['tables'].append(report('FDM initial-value solution',['i','x','y','y prime (given initially; estimated elsewhere)'],[[i,float(x),float(y),float(v)] for i,(x,y,v) in enumerate(zip(xs,ys,ds))]))
    out['tables'].append(report('FDM recurrence evaluations',['i','x_i','P_i','Q_i','R_i','y_(i-1)','y_i','y_(i+1)'],steps))
    out['charts'].append({'title':'FDM initial-value solution','ylabel':'y','curves':[curve('FDM IVP',xs,ys,'nodes')]})
    out['datasets'].append({'label':'FDM IVP','x':xs.tolist(),'y':ys.tolist()})
    add_reference(out,xs,ys,p.get('reference',''),'FDM IVP')
    return out


def comparison(p):
    out=output()
    lines=[s.strip() for s in p['data'].splitlines() if s.strip()]
    if not lines:raise ValueError('Enter a header and rows: x, reference, estimate1, estimate2, ...')
    def split(s):return [v.strip() for v in (s.split(',') if ',' in s else s.split())]
    headers=split(lines[0])
    if len(headers)<3 or len(headers)>14 or len(set(headers))!=len(headers):raise ValueError('Use distinct headers: x, reference, then 1–12 numerical series.')
    data=[]
    for line in lines[1:]:
        cells=split(line)
        if len(cells)!=len(headers):raise ValueError('Every data row must have the same number of columns as the header.')
        data.append([number(v) for v in cells])
    if not data or len(data)>10000:raise ValueError('Use 1–10000 data rows.')
    data.sort();values=np.array(data);xs=values[:,0];refs=values[:,1]
    if len(set(xs))!=len(xs):raise ValueError('Use distinct x nodes.')
    equation(out,'Error definitions','E = reference - numerical; absolute error = |E|; signed error % = 100*E/reference (undefined at zero reference)',
             r'E_i=r_i-u_i,\quad |E_i|=|r_i-u_i|,\quad \varepsilon_i=100\frac{r_i-u_i}{r_i}\%\quad(r_i\ne0)')
    if p.get('equations_only'):return out
    out['tables'].append(report('Comparison input (reference retained)',headers,data))
    curves=[curve(headers[1],xs,refs,'nodes')];ec=[];pc=[];summary=[]
    for j,label in enumerate(headers[2:],2):
        ys=values[:,j];errors=refs-ys;nonzero=refs!=0
        rows=[[float(x),float(r),float(y),float(e),float(abs(e)),float(100*e/r) if r else None] for x,r,y,e in zip(xs,refs,ys,errors)]
        out['tables'].append(report(label+' errors',['x','reference','numerical','signed error','absolute error','signed error %'],rows))
        summary.append([label,len(xs),float(max(abs(errors))),float(np.sqrt(np.mean(errors**2))),float(max(abs(100*errors[nonzero]/refs[nonzero]))) if nonzero.any() else None])
        curves.append(curve(label,xs,ys,'nodes'));ec.append(curve(label,xs,errors,'nodes'))
        if nonzero.any():pc.append(curve(label,xs[nonzero],100*errors[nonzero]/refs[nonzero],'nodes'))
    out['tables'].append(report('Accuracy summary at supplied nodes',['series','nodes','maximum absolute error','RMSE','maximum absolute error % (nonzero references)'],summary))
    out['charts']=[{'title':'Numerical solutions and reference','ylabel':'Value','curves':curves},
                   {'title':'Signed errors at supplied nodes','ylabel':'Reference - numerical','curves':ec}]
    if pc:out['charts'].append({'title':'Percentage errors at supplied nodes','ylabel':'Error (%)','error_plot':True,'curves':pc})
    out['notes'].append('All series are compared at the supplied common nodes. No interpolation or reference substitution is performed. Zero-reference percentage errors are N/A and omitted from percentage plots. Rounded references limit meaningful error precision. Error statistics alone do not establish the cause of a discrepancy.')
    return out
