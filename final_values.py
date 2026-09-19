"""Memory-constant final-value calculations and explicit reference handling."""
import math
from decimal import Decimal, InvalidOperation
import sympy as sp


def meets(error, tolerance, comparison):
    if error is None:
        return False
    return error < tolerance if comparison == '<' else error <= tolerance


def precision_warning(literal, tolerance):
    """Flag possible rounding, without claiming the supplied number is rounded."""
    try:
        decimal=Decimal(literal)
        if not decimal.is_finite() or decimal==0:return None
        uncertainty=Decimal('0.5')*Decimal(10)**decimal.as_tuple().exponent
        percent=float(uncertainty/abs(decimal)*100)
        if percent>=tolerance:
            return (f'If the supplied value {literal} is rounded to its displayed last digit, its possible rounding uncertainty is {percent:.8g}%, at least the requested {tolerance:g}%. '
                    'More iterations cannot certify that reference precision. It is being used unchanged; explicitly select a more precise value or formula if available.')
    except (InvalidOperation,ValueError,OverflowError):pass
    return None


def calculate(p):
    from numerics import number, integer, expr, table, node_value, ode_step, finite, METHODS, report, curve
    x0,xf,y0=(number(p[k]) for k in ('x0','xf','y0'))
    if xf<=x0:raise ValueError('Use final x > initial x.')
    methods=p['methods']
    if not methods or any(m not in METHODS for m in methods):raise ValueError('Select a supported method.')
    comparison=p.get('comparison','<')
    if comparison not in ('<','<='):raise ValueError('Select < or <= for the stopping comparison.')
    mode=p.get('grid_mode','Step sizes h')
    if mode=='Powers of ten n=10^p':
        first=integer(p.get('first_p','1'),0);last=integer(p.get('last_p','6'),0)
        if first>last or last>7:raise ValueError('Use 0 <= first p <= last p <= 7.')
        counts=[10**i for i in range(first,last+1)]
    elif mode=='Interval counts n':
        counts=sorted(set(integer(v.strip()) for v in p['steps'].split(',')))
    elif mode in ('Step sizes h','Refinement h0 / 2^m'):
        if mode=='Step sizes h':hs=[number(v.strip()) for v in p['steps'].split(',')]
        else:
            last=integer(p.get('max_m','3'),0)
            if last>24:raise ValueError('Use last m <= 24.')
            hs=[number(p['steps'])/2**m for m in range(last+1)]
        counts=[]
        for h in hs:
            if h<=0:raise ValueError('Use positive step sizes.')
            n=round((xf-x0)/h)
            if n<1 or not math.isclose(n*h,xf-x0,rel_tol=1e-9,abs_tol=1e-12):raise ValueError('Each h must divide the interval.')
            counts.append(n)
        counts=sorted(set(counts))
    else:raise ValueError('Final-values mode supports h, n, halving, or powers-of-ten grids. Use full tables for candidate n-search.')
    if max(counts)>10000000:raise ValueError('Use at most 10,000,000 steps per grid in final-values mode.')
    passes=sorted(set(integer(v.strip()) for v in p.get('iterations','1').split(',')))
    if max(passes)>1000:raise ValueError('Use at most 1000 PC passes.')
    metric=p.get('stop_metric','No early stop')
    if metric not in ('No early stop','Final true error %','Final successive error %'):
        raise ValueError('Final-values mode supports final true or final successive error, not maximum errors over a grid.')
    tol=number(p.get('stop_tolerance','1'))
    if tol<=0:raise ValueError('Tolerance must be positive percent.')
    output={'tables':[],'charts':[],'datasets':[],'notes':[], 'warnings':[]}
    reference=None;literal=None
    exact_text=p.get('exact','').strip();source=p.get('final_reference_source','Exact formula if provided')
    if exact_text and source!='Given final reference table':
        exact=sp.lambdify(sp.Symbol('x'),expr(exact_text),'math')
        if not math.isclose(finite(exact(x0)),y0,rel_tol=1e-6,abs_tol=1e-8):raise ValueError('Exact formula does not match the initial condition.')
        reference=finite(exact(xf))
        output['notes'].append('Reference: user-entered exact formula '+exact_text+f'; y({xf:g})={reference:.17g}.')
    elif p.get('reference','').strip():
        refs=table(p['reference']);reference=node_value(refs[:,0],refs[:,1],xf)
        for line in p['reference'].splitlines():
            cells=line.split(',') if ',' in line else line.split()
            if len(cells)==2 and math.isclose(number(cells[0]),xf,rel_tol=1e-9,abs_tol=1e-10):literal=cells[1].strip();break
        output['notes'].append(f'Reference: given table, final value {literal}. Used unchanged.')
    if metric=='Final true error %' and reference in (None,0):
        raise ValueError('Final true percentage-error stopping needs a nonzero reference at xf. Enter an exact formula or an xf,y reference row.')
    if literal and metric=='Final true error %':
        warning=precision_warning(literal,tol)
        if warning:output['warnings'].append(warning)
    f=sp.lambdify(sp.symbols('x y'),expr(p['formula'],'x y'),'math');finite(f(x0,y0))
    output['notes'].append(f"ODE: y'={p['formula']}; y({x0:g})={y0:g}; xf={xf:g}. Final values only: no node arrays are retained.")
    output['notes'].append(f'Stopping: {metric}, comparison {comparison}, tolerance {tol:g}%. Signed true error = (reference - final y)/reference * 100; decisions use its absolute value.')
    rows=[];charts=[];total=0
    for method in methods:
        for count in (passes if method=='Predictor-Corrector' else [1]):
            previous=None;points=[];converged=False
            label=method+(f' ({count} PC passes)' if method=='Predictor-Corrector' else '')
            for n in counts:
                total+=n
                if total>50000000:raise ValueError('This run exceeds 50 million total steps; split the methods or grids into separate runs.')
                h=(xf-x0)/n;y=y0;failure=''
                try:
                    for i in range(n):y=finite(ode_step(f,x0+i*h,y,h,method,count))
                except (OverflowError,ValueError,ZeroDivisionError) as error:
                    y=None;failure='Nonfinite/domain failure: '+str(error)
                signed=(reference-y)/reference*100 if reference not in (None,0) and y is not None else None
                if signed is not None and not math.isfinite(signed):signed=None;failure+=' Relative error overflow.'
                absolute=abs(signed) if signed is not None else None
                successive=abs((y-previous)/y)*100 if y not in (None,0) and previous is not None else None
                if successive is not None and not math.isfinite(successive):successive=None
                metric_value=absolute if metric=='Final true error %' else successive
                passed=meets(metric_value,tol,comparison) if metric!='No early stop' else None
                rows.append([label,n,h,y,signed,absolute,successive,passed,failure])
                if signed is not None:points.append((n,signed))
                previous=y
                if passed:
                    output['notes'].append(f'{label}: first qualifying tested n={n}, final y={y:.16g}.')
                    converged=True;break
            if metric!='No early stop' and not converged:output['notes'].append(f'{label}: NOT CONVERGED within the requested grid range.')
            if points:charts.append(curve(label,[x for x,y in points],[y for x,y in points],'nodes'))
    output['tables'].append(report('Final-value trials',['method','n','h','final y','signed true error %','absolute true error %','final successive error %','meets tolerance','note'],rows))
    if charts:output['charts'].append({'title':'Final true error vs n','xlabel':'n','ylabel':'Signed true error (%)','xscale':'log','yscale':'symlog','error_plot':True,'tolerance':tol if metric=='Final true error %' else None,'curves':charts})
    output['notes'].extend(output['warnings'])
    output['notes'].append('No analytical reference was inferred or substituted. A small error relative to a rounded reference alone does not certify accuracy relative to the mathematical solution.')
    return output
