"""Symbolic weak forms for regular scalar linear second-order equations in 1D.

No mesh, numerical solve, or automatic integrating-factor transformation.
Boundary derivatives use increasing x, not outward-normal derivatives.
"""
import sympy as sp
from sympy.calculus.util import continuous_domain
from numerics import expr

x=sp.Symbol('x')
Y,DY,V,DV=sp.symbols("y y' v v'")


def derive(p):
    out=dict(notes=[],tables=[],charts=[],datasets=[],equations=[],warnings=[])
    def add(title,explanation,latex):
        out['equations'].append(dict(title=title,text=explanation,latex=latex))
        out['notes'].append(title+'\n'+explanation+'\nLaTeX: '+latex)
    def constant(s):
        q=expr(s,'')
        if q.is_real is not True or q.is_finite is not True:
            raise ValueError('Endpoints and boundary coefficients must be finite real constants.')
        return q
    left,right=constant(p['left']),constant(p['right'])
    if not bool(left<right):raise ValueError('Left endpoint must be less than right endpoint.')
    interval=sp.Interval(left,right)
    mode=p.get('form','Standard')
    if mode not in ('Standard','Divergence'):raise ValueError('Choose Standard or Divergence form.')
    coeff=[expr(p[k]) for k in ('a','b','c','f')]
    if any(e.has(sp.Piecewise,sp.DiracDelta,sp.Heaviside,sp.sign,sp.Abs) for e in coeff):
        raise ValueError('This derivation requires regular smooth coefficients. Piecewise/distributional coefficients and interface terms are outside this workflow.')
    aa,bb,cc,ff=coeff
    # Check the domain, including endpoints where boundary traces are evaluated.
    for e in coeff+[sp.diff(aa,x)]:
        try:
            domain=continuous_domain(e,x,sp.S.Reals)
            covered=interval.is_subset(domain)
        except (NotImplementedError,ValueError):covered=None
        if covered is False:raise ValueError('A coefficient or its required derivative is singular/nonreal on the closed interval.')
        if covered is not True:
            out['warnings'].append('Symbolic domain checking was inconclusive. The formulas are conditional on regular real coefficients over the closed interval.')
        for pt in (left,right):
            value=sp.simplify(e.subs(x,pt))
            if value.is_finite is not True or value.is_real is not True:
                raise ValueError('A coefficient or derivative is not finite and real at an endpoint.')
    try:
        zeros=sp.solveset(aa,x,domain=interval)
        if zeros.is_empty is False:raise ValueError('The leading coefficient vanishes on this interval. Degenerate problems require a separate weak-space analysis.')
        if zeros.is_empty is None:out['warnings'].append('Could not establish that the leading coefficient is nonzero throughout the interval; verify this regularity assumption.')
    except NotImplementedError:out['warnings'].append('Verify that the leading coefficient is nonzero throughout the interval.')
    L=sp.latex
    def integral(e):return rf'\int_{{{L(left)}}}^{{{L(right)}}}\left({L(sp.expand(e))}\right)\,dx'
    def bracket(e):return rf'\left[{L(e)}\right]_{{{L(left)}}}^{{{L(right)}}}'
    hcoef=aa if mode=='Standard' else -aa
    drift=sp.simplify(bb-sp.diff(aa,x)) if mode=='Standard' else bb
    integrand=sp.expand(-hcoef*DY*DV+drift*DY*V+cc*Y*V)
    unknown=sp.Function('y')(x)
    strong=aa*sp.diff(unknown,x,2)+bb*sp.diff(unknown,x)+cc*unknown if mode=='Standard' else -sp.Derivative(aa*sp.diff(unknown,x),x,evaluate=False)+bb*sp.diff(unknown,x)+cc*unknown
    add('1. Strong form','The entered equation is used as written; no equation scaling or integrating factor is introduced.',L(strong)+'='+L(ff))
    boundaries=[];essential=[]
    for side,point,normal in [('left',left,-1),('right',right,1)]:
        al,be,ga=[constant(p[side+'_'+k]) for k in ('alpha','beta','gamma')]
        if al.is_zero and be.is_zero:raise ValueError('Each boundary needs a nonzero alpha or beta.')
        boundaries.append((point,normal,al,be,ga))
        kind='Dirichlet (essential)' if be.is_zero else ('Neumann (natural)' if al.is_zero else 'Robin (natural)')
        yp=sp.Symbol('y('+L(point)+')');dp=sp.Symbol("y'("+L(point)+')')
        add('Boundary at '+str(point)+' - '+kind,'Boundary convention: alpha*y + beta*y\' = gamma; derivatives point in increasing x.',L(al*yp+be*dp)+'='+L(ga))
        if be.is_zero:essential.append((point,sp.simplify(ga/al)))
    restrictions=r',\quad '.join(rf'y({L(pt)})={L(val)}' for pt,val in essential)
    tests=r',\quad '.join(rf'v({L(pt)})=0' for pt,val in essential)
    base=rf'H^1({L(left)},{L(right)})'
    add('2. Trial functions','Only essential (Dirichlet) values are imposed on trial functions. Natural conditions enter through the boundary term.',
        rf'\mathcal{{U}}=\left\{{y\in {base}\;:\;{restrictions}\right\}}' if essential else rf'\mathcal{{U}}={base}')
    add('Test functions','Variations vanish at Dirichlet endpoints, including where the prescribed solution value is nonzero. They remain free at natural endpoints.',
        rf'\mathcal{{V}}=\left\{{v\in {base}\;:\;{tests}\right\}}' if essential else rf'\mathcal{{V}}={base}')
    if mode=='Standard':
        add('3. Weighted residual','Multiply by an arbitrary admissible test function and integrate.',integral(V*(aa*sp.Symbol("y''")+bb*DY+cc*Y))+'='+integral(ff*V))
        add('4. Integration by parts','Transfer one derivative from y to a*v. The derivative of a must be included.',
            integral(V*aa*sp.Symbol("y''"))+'='+bracket(aa*DY*V)+'-'+integral(DY*(aa*DV+sp.diff(aa,x)*V)))
        add('Coefficient derivative','This coefficient derivative contributes to the first-derivative term.',r'a^\prime(x)='+L(sp.diff(aa,x)))
    else:
        add('3. Weighted residual','Multiply the divergence-form equation by an arbitrary admissible test function and integrate.',
            '-'+integral(V*sp.Symbol("(p y')'"))+'+'+integral(bb*DY*V+cc*Y)+'='+integral(ff*V))
        add('4. Integration by parts','Integrate the entire flux derivative. Do not omit its boundary term.',
            '-'+integral(V*sp.Symbol("(p y')'"))+'=-'+bracket(aa*DY*V)+'+'+integral(aa*DY*DV))
        add('Flux coefficient','The p in the flux derivative denotes this entered coefficient.',r'p(x)='+L(aa))
    add('5. Weak identity before boundary substitution','The volume integrals contain only first derivatives of y and v.',integral(integrand)+'+'+bracket(hcoef*DY*V)+'='+integral(ff*V))
    bilinear=sp.Integer(0);load=sp.Integer(0);boundary_rows=[]
    for point,normal,al,be,ga in boundaries:
        vp=sp.Symbol('v('+L(point)+')');yp=sp.Symbol('y('+L(point)+')')
        H=sp.simplify(hcoef.subs(x,point))
        if be.is_zero:
            boundary_rows.append([str(point),'Dirichlet','0','0','v=0'])
            add('6. Boundary contribution at '+str(point),'The test function is zero here, so this endpoint term vanishes.',rf'v({L(point)})=0\quad\Longrightarrow\quad B_{{{L(point)}}}=0')
        else:
            dp=sp.Symbol("y'("+L(point)+')')
            replacement=sp.simplify((ga-al*yp)/be)
            add('6. Boundary substitution at '+str(point),'Use the prescribed increasing-x derivative condition; the left endpoint has a minus sign from evaluation of the bracket.',L(dp)+'='+L(replacement))
            bc=sp.expand(normal*H*replacement*vp)
            add('Endpoint contribution at '+str(point),'Keep the part involving y on the left and move the known part to the right.',rf'B_{{{L(point)}}}='+L(bc))
            bilinear+=-normal*H*al/be*yp*vp
            load+=-normal*H*ga/be*vp
            boundary_rows.append([str(point),'Neumann' if al.is_zero else 'Robin',str(sp.simplify(-normal*H*al/be)),str(sp.simplify(-normal*H*ga/be)),'v free'])
    def extra(e):return '' if sp.simplify(e)==0 else r'+\left('+L(sp.simplify(e))+r'\right)'
    add('7. Bilinear form','A(y,v) collects the volume terms and the unknown Robin boundary terms.',r'A(y,v)='+integral(integrand)+extra(bilinear))
    add('Load functional','L(v) collects the volume forcing and the known natural-boundary contributions.',r'L(v)='+integral(ff*V)+extra(load))
    add('8. Final weak statement','Find y in the trial space such that this identity holds for every admissible test function.',r'y\in\mathcal{U}:\quad A(y,v)=L(v)\quad\forall v\in\mathcal{V}')
    add('Expanded final weak form','This is the boundary-substituted weak equation for the original input formulation.',integral(integrand)+extra(bilinear)+'='+integral(ff*V)+extra(load))
    out['notes'].append('Scope: regular scalar linear second-order equations in one space dimension. H^1 means square-integrable functions with square-integrable weak first derivatives. This derivation does not prove existence, uniqueness or coercivity, construct finite elements, or solve the equation. Pure Neumann problems may need compatibility and normalization conditions. Equivalent weak forms can differ by multiplying the whole identity by a nonzero constant.')
    out['notes'].extend('WARNING: '+w for w in out['warnings'])
    out['tables'].append({'title':'Weak boundary contributions','columns':['endpoint','type','A coefficient of y*v','L coefficient of v','test trace'],'rows':boundary_rows})
    out['weak_form']={'volume':str(integrand),'boundary_bilinear':str(sp.simplify(bilinear)),'boundary_load':str(sp.simplify(load)),'forcing':str(ff*V),'volume_coefficients':[str(-hcoef),str(drift),str(cc)],'rhs':str(ff)}
    return out
