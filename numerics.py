"""General numerical methods. No problem-specific data or solutions.

Functions return JSON-compatible tables and curves for the desktop interface.
Run as a worker process to keep calculations away from Tk's event loop.
"""
import json
import math
import sys
import numpy as np
import sympy as sp
from numpy.polynomial import Polynomial

METHODS = ['Euler', 'Predictor-Corrector', 'RK2.1 (Heun)', 'RK2.2 (lesson Ralston)', 'RK3', 'RK4']


def expr(text, variables='x'):
    symbols = {name: sp.Symbol(name) for name in variables.split()}
    value = sp.sympify(str(text).replace('^', '**'), locals={**symbols, 'e': sp.E, 'ln': sp.log})
    if not isinstance(value, sp.Expr) or value.free_symbols - set(symbols.values()):
        raise ValueError('Use only '+variables+' in the expression. Enter the right side, without an equals sign.')
    if value.atoms(sp.core.function.AppliedUndef):
        raise ValueError('Unknown function in formula. Use functions such as sin(x), exp(x), sqrt(x), log(x).')
    return value


def number(text):
    value = float(expr(text, ''))
    if not math.isfinite(value):
        raise ValueError('All numerical inputs must be finite.')
    return value


def integer(text, minimum=1):
    value = number(text)
    if value != int(value) or value < minimum:
        raise ValueError(f'Enter an integer at least {minimum}.')
    return int(value)


def table(text, minimum=1):
    rows = []
    for line in text.strip().splitlines():
        if not line.strip():
            continue
        cells = line.split(',') if ',' in line else line.split()
        if len(cells) != 2:
            raise ValueError('Enter two columns per line: x, y (no header).')
        rows.append([number(v.strip()) for v in cells])
    rows.sort()
    if len(rows) < minimum or len({r[0] for r in rows}) != len(rows):
        raise ValueError(f'Use at least {minimum} points with distinct x-values.')
    return np.array(rows, dtype=float)


def finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('A nonfinite result occurred. Check the formula domain and step size.')
    return value


def curve(label, x, y, style='line'):
    return {'label': label, 'x': [finite(v) for v in x], 'y': [finite(v) for v in y], 'style': style}


def report(title, columns, rows):
    return {'title': title, 'columns': columns, 'rows': rows}


def ode_step(f, x, y, h, method, iterations):
    # Slopes k_i are derivatives; the step size is applied in the update.
    k1 = f(x, y)
    if method == 'Euler':
        return y + h*k1
    if method == 'Predictor-Corrector':
        predicted = y + h*k1
        corrected = predicted
        for _ in range(iterations):
            corrected = y + h/2*(k1 + f(x+h, corrected))
        return corrected
    if method == 'RK2.1 (Heun)':
        k2 = f(x+h, y+h*k1)
        return y + h/2*(k1+k2)
    if method == 'RK2.2 (lesson Ralston)':
        k2 = f(x+3*h/4, y+3*h*k1/4)
        return y + h*(k1/3+2*k2/3)
    if method == 'RK3':
        k2 = f(x+h/2, y+h*k1/2)
        k3 = f(x+h, y-h*k1+2*h*k2)
        return y + h/6*(k1+4*k2+k3)
    if method == 'RK4':
        k2 = f(x+h/2, y+h*k1/2)
        k3 = f(x+h/2, y+h*k2/2)
        k4 = f(x+h, y+h*k3)
        return y + h/6*(k1+2*k2+2*k3+k4)
    raise ValueError('Choose a supported ODE method.')


def node_value(xs, ys, target):
    index = int(np.argmin(np.abs(xs-target)))
    tolerance = 1e-9 * max(1, abs(target))
    return float(ys[index]) if abs(xs[index]-target) <= tolerance else None


def sample_value(xs,ys,target,mode='Grid nodes only'):
    value=node_value(xs,ys,target)
    if value is not None:return value
    if target<xs[0] or target>xs[-1]:return None
    if mode=='Linear interpolation':return float(np.interp(target,xs,ys))
    if mode=='Polynomial interpolation':
        if len(xs)>31:raise ValueError('Polynomial resampling supports at most 31 nodes. Choose linear interpolation or fewer nodes.')
        poly,info=Polynomial.fit(xs,ys,len(xs)-1,full=True)
        if info[1]<len(xs):raise ValueError('Interpolation is rank deficient; choose fewer nodes.')
        return finite(poly(target))
    return None


def solve(p):
    x0, xf, y0 = (number(p[k]) for k in ('x0', 'xf', 'y0'))
    if xf <= x0:
        raise ValueError('xf must be greater than x0.')
    mode=p.get('grid_mode','Step sizes h')
    if mode=='Interval counts n':
        hs=sorted(set((xf-x0)/integer(v.strip()) for v in p['steps'].split(',')),reverse=True)
    elif mode=='Refinement h0 / 2^m':
        base=number(p['steps']); last=integer(p.get('max_m','3'),0)
        if last>20:raise ValueError('Use m at most 20 and at most 100,000 total steps.')
        hs=[base/2**m for m in range(last+1)]
    else:
        hs = sorted(set(number(v.strip()) for v in p['steps'].split(',')), reverse=True)
    methods = p['methods']
    if not methods or any(m not in METHODS for m in methods):
        raise ValueError('Select at least one method.')
    counts = []
    for h in hs:
        if h <= 0:
            raise ValueError('Step sizes must be positive.')
        n = round((xf-x0)/h)
        if n < 1 or not math.isclose(n*h, xf-x0, rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError(f'h={h:g} does not divide the interval. Use (xf-x0)/n.')
        counts.append(n)
    if sum(counts)*len(methods) > 100000:
        raise ValueError('Use at most 100,000 total steps per run.')
    iterations = integer(p.get('iterations', '1'))
    if iterations > 1000:
        raise ValueError('Use at most 1000 corrector passes.')
    f = sp.lambdify(sp.symbols('x y'), expr(p['formula'], 'x y'), 'math')
    finite(f(x0, y0))
    exact_text = p.get('exact', '').strip()
    exact = sp.lambdify(sp.Symbol('x'), expr(exact_text), 'math') if exact_text else None
    references = table(p['reference']) if p.get('reference', '').strip() else None
    output = {'tables': [], 'charts': [], 'datasets': [], 'notes': []}
    curves = []
    error_curves = []
    summaries=[]
    requested=[number(v.strip()) for v in p['requested'].split(',')] if p.get('requested','').strip() else []
    sampling=p.get('sampling','Grid nodes only')
    output['notes'].append(f"ODE: y' = {p['formula']}; y({x0:g}) = {y0:g}; interval [{x0:g}, {xf:g}].")
    if 'Predictor-Corrector' in methods:
        output['notes'].append(f'Predictor-Corrector: Euler predictor, trapezoidal corrector, {iterations} pass(es) per step.')
    if exact and not math.isclose(finite(exact(x0)), y0, rel_tol=1e-6, abs_tol=1e-8):
        raise ValueError('The exact-solution formula does not match y(x0). A given final value belongs in the reference table as xf, y(xf).')
    for method in methods:
        previous = None
        for level,(h, n) in enumerate(zip(hs, counts)):
            xs = np.linspace(x0, xf, n+1)
            ys = np.empty(n+1)
            ys[0] = y0
            for i in range(n):
                try:
                    ys[i+1] = finite(ode_step(f, float(xs[i]), float(ys[i]), (xf-x0)/n, method, iterations))
                except (ValueError, OverflowError, ZeroDivisionError) as error:
                    raise ValueError(f'{method}, h={h:g}, step {i+1}: {error}') from error
            label = f'{method}, h={h:g}'
            output['datasets'].append({'label':label, 'x':xs.tolist(), 'y':ys.tolist()})
            curves.append(curve(label, xs, ys, 'nodes'))
            output['tables'].append(report(label, ['i', 'x', 'y numerical'], [[i,float(x),float(y)] for i,(x,y) in enumerate(zip(xs,ys))]))
            if requested:
                output['tables'].append(report(label+' — requested x-values ('+sampling+')',['x','y','source'],
                    [[x,sample_value(xs,ys,x,sampling),'grid node' if node_value(xs,ys,x) is not None else sampling] for x in requested]))
            reference_sets = []
            if exact:
                reference_sets.append(('Exact formula', [(float(x),finite(exact(float(x)))) for x in xs]))
            if references is not None:
                reference_sets.append(('Given exact values', references.tolist()))
            for source, pairs in reference_sets:
                rows, ex, ey = [], [], []
                for x, reference in pairs:
                    numerical = sample_value(xs, ys, x, sampling)
                    absolute = abs(reference-numerical) if numerical is not None else None
                    percent = absolute/abs(reference)*100 if absolute is not None and reference != 0 else None
                    rows.append([x,numerical,reference,absolute,percent,
                                 'not available on grid/domain' if numerical is None else ('relative error undefined: exact=0' if reference==0 else ('interpolated: '+sampling if node_value(xs,ys,x) is None else ''))])
                    if percent is not None:
                        ex.append(x); ey.append(percent)
                output['tables'].append(report(label+' — '+source,
                    ['x','y numerical','y exact','absolute error','absolute error %','note'], rows))
                if ex:
                    error_curves.append(curve(label+' / '+source, ex, ey, 'nodes'))
            if previous is not None:
                coarse_x, coarse_y, coarse_h = previous
                rows = []
                comparison_x=requested if requested else coarse_x
                for x in comparison_x:
                    old=sample_value(coarse_x,coarse_y,x,sampling)
                    current = sample_value(xs, ys, x,sampling)
                    difference = abs(current-old) if current is not None and old is not None else None
                    percent = difference/abs(current)*100 if current not in (None,0) and difference is not None else None
                    rows.append([float(x),old,current,difference,percent])
                output['tables'].append(report(f'{method} refinement: h={coarse_h:g} to {h:g}',
                    ['x (coarse grid)','coarse y','fine y','absolute difference','approx. relative error %'],rows))
            true_percent=[]
            final_error=None
            if exact:
                for x,y in zip(xs,ys):
                    ref=finite(exact(float(x)))
                    if ref!=0:true_percent.append(abs((ref-y)/ref)*100)
                ref=finite(exact(xf))
                if ref!=0:final_error=abs((ref-ys[-1])/ref)*100
            elif references is not None:
                for x,ref in references:
                    val=sample_value(xs,ys,x,sampling)
                    if val is not None and ref!=0:true_percent.append(abs((ref-val)/ref)*100)
                    if math.isclose(x,xf) and ref!=0:final_error=abs((ref-ys[-1])/ref)*100
            approx=[]
            final_successive=None
            if previous is not None:
                if ys[-1]!=0:final_successive=abs((ys[-1]-previous[1][-1])/ys[-1])*100
                for x,old in zip(previous[0],previous[1]):
                    val=sample_value(xs,ys,x,sampling)
                    if val not in (None,0):approx.append(abs((val-old)/val)*100)
            summaries.append([method,level,n,h,float(ys[-1]),finite(max(true_percent)) if true_percent else None,
                              finite(final_error) if final_error is not None else None,finite(max(approx)) if approx else None])
            previous = (xs,ys,h)
            criterion=p.get('stop_metric','No early stop')
            if criterion!='No early stop':
                threshold=number(p.get('stop_tolerance','1'))
                if threshold<=0:raise ValueError('Error tolerance must be positive (percent).')
                value={'Maximum true error %':summaries[-1][5],'Final true error %':summaries[-1][6],
                       'Maximum successive error %':summaries[-1][7], 'Final successive error %':final_successive}.get(criterion)
                if value is not None and (value<threshold if p.get('comparison','<=')=='<' else value<=threshold):
                    output['notes'].append(f'{method}: first tested grid meeting {criterion} {p.get("comparison","<=")} {threshold:g}% has n={n}, h={h:g}. This does not establish a minimum n among untested grids.')
                    break
        else:
            if p.get('stop_metric','No early stop')!='No early stop':output['notes'].append(method+': requested tolerance was not reached or its reference values were unavailable on the tested grids.')
    if exact:
        xp = np.linspace(x0,xf,400)
        curves.append(curve('Exact formula',xp,[exact(float(x)) for x in xp]))
    if references is not None:
        curves.append(curve('Given exact values', references[:,0],references[:,1], 'points'))
    output['charts'].append({'title':'ODE solutions','ylabel':'y','curves':curves})
    if error_curves:
        output['charts'].append({'title':'Errors against supplied references','ylabel':'Absolute relative error (%)','curves':error_curves})
    output['tables'].append(report('Grid convergence summary',['method','grid level m (if halving)','n','h','final y','max true error %','final true error %','max successive error %'],summaries))
    for metric,column in [('Maximum true error (%)',5),('Final true error (%)',6),('Maximum successive error (%)',7)]:
        cc=[]
        for method in methods:
            rows=[r for r in summaries if r[0]==method and r[column] is not None]
            if rows:cc.append(curve(method,[r[1] if mode=='Refinement h0 / 2^m' else r[2] for r in rows],[r[column] for r in rows],'nodes'))
        if cc:output['charts'].append({'title':metric+' vs '+('m' if mode=='Refinement h0 / 2^m' else 'n'),'xlabel':'m' if mode=='Refinement h0 / 2^m' else 'n','ylabel':metric,'curves':cc})
    output['notes'].append('Off-grid evaluation policy: '+sampling+'. Percent error is undefined when its denominator is zero. Maximum true error uses only supplied reference points (or all grid nodes with an exact formula), excluding zero references.')
    if len(hs)>1:
        output['notes'].append('Refinement error = |(fine - coarse)/fine| × 100 at shared nodes. These are successive-grid differences, not exact errors or Richardson estimates.')
    return output


def secant(g, a, b, tolerance, limit, criterion='Residual'):
    rows = []
    for k, x in enumerate((a,b)):
        rows.append([k,x,finite(g(x))])
    if rows[0][2]==0 or (criterion=='Residual' and abs(rows[0][2]) < tolerance):
        return a, rows[:1], True
    if rows[1][2]==0 or (criterion=='Residual' and abs(rows[1][2]) < tolerance):
        return b, rows, True
    for k in range(2,limit+2):
        ga, gb = finite(g(a)), finite(g(b))
        if gb == ga:
            return b, rows, False
        c = finite(b-gb*(b-a)/(gb-ga))
        gc = finite(g(c))
        rows.append([k,c,gc])
        if gc==0 or (abs(gc) < tolerance if criterion=='Residual' else c!=0 and abs((c-b)/c)*100<tolerance):
            return c, rows, True
        a,b = b,c
    return b,rows,False


def analyze(p):
    if p.get('data_source')=='Sample formula':
        low,high=number(p['data_min']),number(p['data_max']);n=integer(p['data_n'])
        if low>=high or n>100000:raise ValueError('Use increasing bounds and at most 100,000 intervals.')
        fun=sp.lambdify(sp.Symbol('x'),expr(p['data_formula']),'math')
        xx=np.linspace(low,high,n+1)
        points=np.array([[x,finite(fun(float(x)))] for x in xx])
    else:
        points = table(p['data'], 2)
    xs,ys = points[:,0],points[:,1]
    operation = p['operation']
    out = {'tables':[], 'charts':[], 'notes':[], 'datasets':[]}
    out['notes'].append(f'Analysis uses the {len(xs)} points currently entered in the Data table, over [{xs[0]:g}, {xs[-1]:g}].')
    out['tables'].append(report('Analysis input data',['x','y'],points.tolist()))
    reference_text=p.get('analysis_exact','').strip()
    reference_expr=expr(reference_text) if reference_text else None
    if operation=='Linearized least squares':
        return transformed_fit(p,xs,ys,out)
    if operation == 'Finite differences':
        scheme=p.get('scheme','Central + higher-order endpoints')
        if len(xs)<4 and scheme=='Central + higher-order endpoints':
            raise ValueError('Use at least 4 points for second-order endpoint second derivatives.')
        rows=[]
        for i,x in enumerate(xs):
            # Three-point interior/first derivative; four-point endpoint second derivative.
            def derivative(order):
                if scheme in ('Forward','Backward'):
                    width=order+1;start=i if scheme=='Forward' else i-width+1
                elif scheme in ('Forward higher-order','Backward higher-order'):
                    width=order+2;start=i if scheme.startswith('Forward') else i-width+1
                elif scheme=='Central only':
                    width=3;start=i-1
                else:
                    width = 4 if order==2 and i in (0,len(xs)-1) else 3
                    start = min(max(i-width//2,0),len(xs)-width)
                if start<0 or start+width>len(xs):return None
                offsets = xs[start:start+width]-x
                scale = max(abs(offsets))
                matrix = np.array([(offsets/scale)**k for k in range(width)])
                target = np.zeros(width); target[order]=math.factorial(order)/scale**order
                weights = np.linalg.solve(matrix,target)
                return finite(weights @ ys[start:start+width])
            rows.append([float(x),float(ys[i]),derivative(1),derivative(2)])
        out['tables'].append(report('Finite-difference derivatives',['x','y',"dy/dx","d²y/dx²"],rows))
        combined=[]
        for title,j in [('First derivative',2),('Second derivative',3)]:
            valid=[r for r in rows if r[j] is not None]
            c=curve(title,[r[0] for r in valid],[r[j] for r in valid],'nodes');combined.append(c)
            cc=[c]
            if reference_expr is not None:
                ref=sp.lambdify(sp.Symbol('x'),sp.diff(reference_expr,sp.Symbol('x'),j-1),'math')
                comparisons=[]
                for row in valid:
                    actual=finite(ref(row[0]));absolute=abs(actual-row[j]);percent=absolute/abs(actual)*100 if actual else None
                    comparisons.append([row[0],row[j],actual,absolute,percent])
                out['tables'].append(report(title+' errors',['x','numerical','exact','absolute error','absolute error %'],comparisons))
                cc.append(curve('Exact '+title,xs,[ref(float(x)) for x in xs]))
                errors=[r for r in comparisons if r[-1] is not None]
                if errors:out['charts'].append({'title':title+' errors','ylabel':'Error %','curves':[curve(title,[r[0] for r in errors],[r[-1] for r in errors],'nodes')]})
            out['charts'].append({'title':title,'ylabel':title,'curves':cc})
        out['charts'].insert(0,{'title':'First and second derivatives on one graph','ylabel':'Derivative','curves':combined})
        out['notes'].append('Finite-difference scheme: '+scheme+'. Unavailable stencils are N/A. Local polynomial weights reproduce the lesson formulas on uniform grids.')
        return out
    if operation in ('Trapezoidal rule', "Simpson's 1/3 rule"):
        h=np.diff(xs)
        if operation=='Trapezoidal rule':
            integral=float(np.sum(h*(ys[:-1]+ys[1:])/2))
        else:
            if (len(xs)-1)%2 or not np.allclose(h,h[0],rtol=1e-9,atol=1e-12):
                raise ValueError("Simpson's 1/3 rule requires equally spaced x-values and an even number of intervals.")
            integral=float(h[0]/3*(ys[0]+ys[-1]+4*ys[1:-1:2].sum()+2*ys[2:-1:2].sum()))
        out['tables'].append(report(operation,['lower x','upper x','integral'],[[float(xs[0]),float(xs[-1]),integral]]))
        if p.get('exact_integral','').strip():
            ref=number(p['exact_integral'])
            out['tables'].append(report('Integral error',['numerical integral','given exact integral','absolute error','absolute error %'],[[integral,ref,abs(integral-ref),abs((integral-ref)/ref)*100 if ref else None]]))
        out['charts'].append({'title':operation+' — input data','ylabel':'y','curves':[curve('Input data',xs,ys,'nodes')]})
        return out
    fit_kind=p.get('fit','Exact interpolation')
    degree=len(xs)-1 if fit_kind.startswith('Exact') else integer(p['degree'],0)
    if degree>=len(xs):
        raise ValueError('Polynomial degree must be smaller than the number of data points.')
    if degree>30:
        raise ValueError('Use degree ≤ 30. Select fewer points or use a lower-degree least-squares fit.')
    poly, diagnostics=Polynomial.fit(xs,ys,degree,full=True)
    if diagnostics[1] < degree+1:
        raise ValueError('The polynomial system is numerically rank deficient. Use fewer points or a lower degree.')
    if fit_kind=='Exact Lagrange':
        center=(xs[0]+xs[-1])/2;scale=(xs[-1]-xs[0])/2
        nodes=(xs-center)/scale
        coefficients=Polynomial([0.])
        for i in range(len(xs)):
            basis=Polynomial([1.])
            for j in range(len(xs)):
                if i!=j:basis*=Polynomial([-nodes[j],1.])/(nodes[i]-nodes[j])
            coefficients+=ys[i]*basis
        poly=Polynomial(coefficients.coef,domain=[xs[0],xs[-1]])
    elif fit_kind=='Exact Vandermonde':
        nodes=(2*xs-xs[0]-xs[-1])/(xs[-1]-xs[0])
        coefficients=np.linalg.solve(np.vander(nodes,N=len(xs),increasing=True),ys)
        poly=Polynomial(coefficients,domain=[xs[0],xs[-1]])
    coeff=poly.convert().coef
    out['tables'].append(report(f'{fit_kind}: degree {degree} — expanded coefficients',['power of x','coefficient'],[[i,float(v)] for i,v in enumerate(coeff)]))
    terms=' + '.join(f'({c:.16g})*x^{i}' for i,c in enumerate(coeff))
    out['notes'].append('f(x) = '+terms)
    out['notes'].append(f'Maximum residual at supplied points: {max(abs(poly(xs)-ys)):.6g}. Calculations use the scaled polynomial, not the rounded expanded coefficients.')
    out['tables'].append(report('Fit residuals',['x','given y','fitted y','residual'],[[float(x),float(y),finite(poly(x)),finite(poly(x)-y)] for x,y in zip(xs,ys)]))
    out['tables'].append(report('Fit quality',['sum squared residuals'],[[finite(np.sum((poly(xs)-ys)**2))]]))
    lower=number(p['analysis_min']) if p.get('analysis_min','').strip() else float(xs[0])
    upper=number(p['analysis_max']) if p.get('analysis_max','').strip() else float(xs[-1])
    if lower>=upper:raise ValueError('Plot/integration lower bound must be less than upper bound.')
    if lower<xs[0] or upper>xs[-1]:out['notes'].append('Requested domain extends beyond the data: the polynomial is extrapolated there.')
    xp=np.linspace(lower,upper,500)
    out['charts'].append({'title':fit_kind,'ylabel':'y','curves':[curve('Data',xs,ys,'points'),curve('Polynomial',xp,poly(xp))]})
    for order in (0,1,2):
        model=poly.deriv(order);title=['Function','First derivative','Second derivative'][order]
        cc=[curve('Fitted '+title,xp,model(xp))]
        if reference_expr is not None:
            ref=sp.lambdify(sp.Symbol('x'),sp.diff(reference_expr,sp.Symbol('x'),order),'math')
            cc.append(curve('Exact '+title,xp,[ref(float(x)) for x in xp]))
        out['charts'].append({'title':title+' comparison','ylabel':title,'curves':cc})
    query=p.get('query','').strip()
    if query:
        values=[number(v.strip()) for v in query.split(',')]
        if any(x<xs[0] or x>xs[-1] for x in values):out['notes'].append('Some requested evaluations extrapolate beyond the data domain.')
        out['tables'].append(report('Polynomial evaluation',['x','f(x)',"f′(x)","f″(x)"],[[x,finite(poly(x)),finite(poly.deriv()(x)),finite(poly.deriv(2)(x))] for x in values]))
    anti=poly.integ()
    integral=finite(anti(upper)-anti(lower))
    out['tables'].append(report('Integral of fitted polynomial',['lower x','upper x','integral'],[[lower,upper,integral]]))
    from polynomial_working import append_working
    append_working(out,poly,xs,ys,lower,upper,fit_kind=='Least squares')
    if p.get('exact_integral','').strip():
        ref=number(p['exact_integral']);out['tables'].append(report('Integral error',['numerical','exact','absolute error','absolute error %'],[[integral,ref,abs(integral-ref),abs((integral-ref)/ref)*100 if ref else None]]))
    if operation=='Stationary point (secant)':
        a,b,tol=number(p['a']),number(p['b']),number(p['tolerance'])
        if tol<=0: raise ValueError('Tolerance must be positive.')
        root,rows,ok=secant(poly.deriv(),a,b,tol,integer(p['limit']))
        out['tables'].append(report("Secant: g(x) = f′(x) = 0",['k','x','g(x)'],rows))
        second=finite(poly.deriv(2)(root))
        inside=xs[0]<=root<=xs[-1]
        classification=('local minimum' if second>0 else 'local maximum' if second<0 else 'inconclusive') if ok and inside else 'not verified within data domain'
        out['tables'].append(report('Stationary-point result',['x','f(x)',"f″(x)",'converged','classification'],[[root,finite(poly(root)),second,ok,classification]]))
        out['notes'].append(f'Stopping criterion: |f′(x)| < {tol:g}. A local minimum requires convergence, a point in the data domain, and f″(x)>0.')
        if ok and inside:
            out['charts'][0]['curves'].append(curve(classification,[root],[poly(root)],'points'))
    return out


def transformed_fit(p,xs,ys,out):
    model=p.get('transformed_model','Exponential: a*exp(b*x)')
    if model.startswith('Exponential'):
        if np.any(ys<=0):raise ValueError('Log-linear fitting requires positive y.')
        tx,ty=xs,np.log(ys)
    elif model.startswith('Power'):
        if np.any(xs<=0) or np.any(ys<=0):raise ValueError('Power-law log fitting requires positive x and y.')
        tx,ty=np.log(xs),np.log(ys)
    elif model.startswith('Reciprocal'):
        if np.any(ys==0):raise ValueError('Reciprocal fitting requires nonzero y.')
        tx,ty=xs,1/ys
    else:
        if np.any(xs==0):raise ValueError('The model a/x+b requires nonzero x.')
        tx,ty=1/xs,ys
    design=np.column_stack([np.ones(len(tx)),tx])
    coefficients,_,rank,_=np.linalg.lstsq(design,ty,rcond=None)
    if rank<2:raise ValueError('Insufficient independent transformed coordinates.')
    a,b=map(float,coefficients);x=sp.Symbol('x')
    if model.startswith('Exponential'):formula=sp.Float(math.exp(a))*sp.exp(b*x)
    elif model.startswith('Power'):formula=sp.Float(math.exp(a))*x**b
    elif model.startswith('Reciprocal'):formula=1/(a+b*x)
    else:formula=a+b/x
    f=sp.lambdify(x,formula,'math')
    fitted=[finite(f(float(v))) for v in xs]
    out['notes'].append('f(x) = '+str(formula))
    out['notes'].append('Linear least squares minimizes residuals in the transformed variables. SSE below is measured in original y for model comparison.')
    out['tables'].append(report('Fit quality',['SSE in original y'],[[finite(sum((a-b)**2 for a,b in zip(fitted,ys)))]]))
    out['tables'].append(report('Fitted values',['x','given y','fitted y'],[[float(x),float(y),z] for x,y,z in zip(xs,ys,fitted)]))
    xp=np.linspace(xs[0],xs[-1],400)
    out['charts'].append({'title':model,'ylabel':'y','curves':[curve('Data',xs,ys,'points'),curve('Fit',xp,[f(float(x)) for x in xp])]})
    return out


def roots(p):
    expression=expr(p['function'])
    g=sp.lambdify(sp.Symbol('x'),expression,'math')
    a,b=number(p['a']),number(p['b'])
    tol=number(p['tolerance']); limit=integer(p['limit'])
    if tol<=0: raise ValueError('Tolerance must be positive.')
    method=p['method']
    criterion=p.get('criterion','Residual')
    rows=[]; ok=False
    if method=='Secant':
        root,rows,ok=secant(g,a,b,tol,limit,criterion)
    elif method=='Bisection':
        if a>=b: raise ValueError('Use a < b for the bracket.')
        ga,gb=finite(g(a)),finite(g(b))
        if ga==0 or (criterion=='Residual' and abs(ga)<tol): root=a; rows=[[0,a,ga]];ok=True
        elif gb==0 or (criterion=='Residual' and abs(gb)<tol): root=b;rows=[[0,b,gb]];ok=True
        elif ga*gb>0: raise ValueError('The bracket must have opposite function signs.')
        else:
            for k in range(limit):
                root=(a+b)/2; gr=finite(g(root));rows.append([k,root,gr])
                if gr==0 or (abs(gr)<tol if criterion=='Residual' else len(rows)>1 and root!=0 and abs((root-rows[-2][1])/root)*100<tol): ok=True;break
                if ga*gr<0: b=root
                else: a=root;ga=gr
    elif method=='Newton-Raphson':
        dg=sp.lambdify(sp.Symbol('x'),sp.diff(expression,sp.Symbol('x')),'math')
        root=a
        for k in range(limit+1):
            gr=finite(g(root)); rows.append([k,root,gr])
            if gr==0 or (abs(gr)<tol if criterion=='Residual' else len(rows)>1 and root!=0 and abs((root-rows[-2][1])/root)*100<tol): ok=True;break
            if k==limit: break
            slope=finite(dg(root))
            if slope==0: break
            root=finite(root-gr/slope)
    else: raise ValueError('Choose a root-finding method.')
    low,high=number(p['plot_min']),number(p['plot_max'])
    if low>=high: raise ValueError('Plot lower bound must be smaller than upper bound.')
    xp=np.linspace(low,high,400)
    # Domain gaps in a plot do not invalidate an otherwise valid root calculation.
    yp=[]
    for x in xp:
        try: yp.append(finite(g(float(x))))
        except (ValueError,OverflowError,ZeroDivisionError): yp.append(None)
    relative=[]
    for i,row in enumerate(rows):
        relative.append(abs((row[1]-rows[i-1][1])/row[1])*100 if i and row[1]!=0 else None)
    return {'tables':[report(method,['k','x','g(x)','approximate relative error %'],[r+[e] for r,e in zip(rows,relative)]),report('Root result',['x','g(x)','converged'],[[root,finite(g(root)),ok]])],
            'charts':[{'title':'Root-finding function','ylabel':'g(x)','curves':[{'label':p['function'],'x':xp.tolist(),'y':yp,'style':'line'},curve('Root estimate',[root],[g(root)],'points')]},
                      {'title':'Root convergence','xlabel':'Iteration k','ylabel':'|g(x)|','curves':[curve(method,[r[0] for r in rows],[abs(r[2]) for r in rows],'nodes')]}],
            'notes':[f'g(x) = {p["function"]}; stopping criterion '+(f'|g(x)| < {tol:g}.' if criterion=='Residual' else f'absolute relative step change < {tol:g}%. A small step alone does not guarantee a small residual; inspect g(x).'), 'Converged.' if ok else 'NOT CONVERGED: inspect starting points, derivative or iteration limit.'], 'datasets':[]}


def search_n(p):
    """Bracket and bisect integer n under an explicit monotone-error assumption."""
    x0,xf,y0=number(p['x0']),number(p['xf']),number(p['y0'])
    if xf<=x0:raise ValueError('Use xf > x0.')
    f=sp.lambdify(sp.symbols('x y'),expr(p['formula'],'x y'),'math')
    if p.get('exact','').strip():
        exact=sp.lambdify(sp.Symbol('x'),expr(p['exact']),'math');reference=finite(exact(xf))
    else:
        refs=table(p.get('reference',''))
        reference=node_value(refs[:,0],refs[:,1],xf)
        if reference is None:raise ValueError('Supply the exact final value as an xf,y row or an exact formula.')
    if reference==0:raise ValueError('Relative-error search requires a nonzero exact final value.')
    start=integer(p['steps']);maximum=integer(p.get('max_n','100000'))
    if start>maximum or maximum>100000:raise ValueError('Use initial n ≤ maximum n ≤ 100,000.')
    tolerance=number(p['stop_tolerance'])
    if tolerance<=0:raise ValueError('Tolerance must be positive percent.')
    passes=integer(p.get('iterations','1'));all_results={'tables':[],'charts':[],'notes':[],'datasets':[]}
    if not p['methods']:raise ValueError('Select at least one method.')
    if passes>1000:raise ValueError('Use at most 1000 PC passes.')
    for method in p['methods']:
        tried={}
        def error(n):
            if n not in tried:
                h=(xf-x0)/n;y=y0
                for i in range(n):y=finite(ode_step(f,x0+i*h,y,h,method,passes))
                tried[n]=(y,abs((y-reference)/reference)*100)
            return tried[n][1]
        high=start
        from final_values import meets
        accepted=lambda n:meets(error(n),tolerance,p.get('comparison','<='))
        while not accepted(high) and high<maximum:high=min(maximum,2*high)
        if not accepted(high):
            all_results['notes'].append(f'{method}: no qualifying grid found up to n={maximum}.')
        else:
            low=1
            while low<high:
                middle=(low+high)//2
                if accepted(middle):high=middle
                else:low=middle+1
            if high>1:error(high-1)
            q=dict(p,grid_mode='Interval counts n',steps=str(high),methods=[method],stop_metric='No early stop')
            result=solve(q)
            for key in all_results:all_results[key].extend(result[key])
            all_results['notes'].append(f'{method}: candidate minimum n={high}, h={(xf-x0)/high:.12g}, final error={error(high):.8g}%.')
        rows=[[n,(xf-x0)/n,y,e] for n,(y,e) in sorted(tried.items())]
        all_results['tables'].append(report(method+' — n-search trials',['n','h','final y','final error %'],rows))
        all_results['charts'].append({'title':method+' final error vs n','xlabel':'n','ylabel':'Error %','curves':[curve(method,[r[0] for r in rows],[r[3] for r in rows],'nodes')]})
    all_results['notes'].append('n-search assumes final error decreases monotonically as n increases. The candidate is not a guaranteed global minimum when errors are nonmonotonic. Trial grids and the predecessor are shown for inspection.')
    return all_results


def dispatch(request):
    if request['action'] in ('systems','fdm_ivp','comparison'):
        import systems_methods
        return getattr(systems_methods,request['action'])(request['inputs'])
    if request['action'] in ('bvp','matrix'):
        import bvp_methods
        result=getattr(bvp_methods,request['action'])(request['inputs'])
        result['notes'].extend('WARNING: '+v for v in result.get('warnings',[]))
        return result
    if request['action']=='solve' and (request['inputs'].get('output_mode')=='Final values only' or request['inputs'].get('grid_mode')=='Powers of ten n=10^p'):
        from final_values import calculate
        return calculate(request['inputs'])
    if request['action']=='solve' and request['inputs'].get('grid_mode')=='Search n (final error)':
        return search_n(request['inputs'])
    if request['action']=='solve' and ',' in request['inputs'].get('iterations','1'):
        params=request['inputs'];combined={'tables':[],'charts':[],'datasets':[],'notes':[]}
        passes=sorted(set(integer(v.strip()) for v in params['iterations'].split(',')))
        for count in passes:
            p=dict(params,iterations=str(count))
            # Other methods do not use the corrector-pass setting.
            if count!=passes[0]:p['methods']=[m for m in p['methods'] if m=='Predictor-Corrector']
            if not p['methods']:continue
            result=solve(p)
            suffix=f' / PC passes={count}'
            for group in ('tables','datasets'):
                for item in result[group]:item['title' if group=='tables' else 'label']+=suffix
                combined[group].extend(result[group])
            for chart in result['charts']:
                for c in chart['curves']:c['label']+=suffix
                existing=next((c for c in combined['charts'] if c['title']==chart['title']),None)
                if existing:existing['curves'].extend(chart['curves'])
                else:combined['charts'].append(chart)
            combined['notes'].extend(result['notes'])
        return combined
    return {'solve':solve,'analyze':analyze,'roots':roots}[request['action']](request['inputs'])


def execute(request):
    result=dispatch(request)
    # Keep existing absolute-error columns and add the signed lesson convention.
    for item in result['tables']:
        columns=item['columns']
        if 'y exact' in columns and 'y numerical' in columns:
            ri,yi=columns.index('y exact'),columns.index('y numerical')
            columns.append('signed true error %')
            for row in item['rows']:
                ref,value=row[ri],row[yi]
                row.append((ref-value)/ref*100 if ref not in (None,0) and value is not None else None)
    return result


if __name__=='__main__':
    try:
        result={'ok':True,'result':execute(json.load(sys.stdin))}
    except Exception as error:
        result={'ok':False,'error':f'{type(error).__name__}: {error}'}
    print(json.dumps(result,allow_nan=False))
