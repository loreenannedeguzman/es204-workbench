"""General 1D linear Galerkin elements with an explicit manual-partition stage.

Canonical convention: -(p y')' + r y = f; outward flux is n*p*y'.
Input/output are JSON-compatible. This module does not access files or networks.
"""
import ast
import re
import numpy as np
import sympy as sp

x = sp.Symbol('x')
xi = sp.Symbol('xi')
FUNCTIONS = {'sin':sp.sin, 'cos':sp.cos, 'tan':sp.tan, 'exp':sp.exp,
             'log':sp.log, 'ln':sp.log, 'sqrt':sp.sqrt, 'abs':sp.Abs}


def expression(source):
    """Restricted AST evaluator: no Python eval, attributes, indexing or imports."""
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int,float):
            return sp.Rational(str(node.value))
        if isinstance(node, ast.Name) and node.id in ('x','pi','e','E'):
            return {'x':x,'pi':sp.pi,'e':sp.E,'E':sp.E}[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op,(ast.USub,ast.UAdd)):
            v=visit(node.operand)
            return -v if isinstance(node.op,ast.USub) else v
        if isinstance(node, ast.BinOp):
            a,b=visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
            if isinstance(node.op,ast.Pow):
                if b.is_number and abs(float(b))>100:raise ValueError('Exponent exceeds 100.')
                return a**b
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id in FUNCTIONS and len(node.args)==1 and not node.keywords:
            return FUNCTIONS[node.func.id](visit(node.args[0]))
        raise ValueError('Use arithmetic, x and supported mathematical functions only.')
    source=str(source).strip().replace('^','**')
    if not source or len(source)>1000:raise ValueError('Enter a mathematical expression of at most 1000 characters.')
    try:return visit(ast.parse(source,mode='eval').body)
    except (SyntaxError,TypeError,ZeroDivisionError) as error:raise ValueError('Invalid mathematical expression.') from error


def number(value):
    try:v=float(sp.N(value,17))
    except (TypeError,ValueError,OverflowError) as error:raise ValueError('Expected a finite real number; check coefficient domains.') from error
    if not np.isfinite(v):raise ValueError('Expected finite real coefficients and values.')
    return v


def values(source):
    return [number(expression(t)) for t in re.split(r'[,\s]+',str(source).strip()) if t]


def equation(out,title,text,latex,**extra):
    out['equations'].append(dict(title=title,text=text,latex=latex,**extra))
    out['notes'].append(title+': '+text)


def table(out,title,columns,rows):
    out['tables'].append(dict(title=title,columns=columns,rows=rows))


def matrix_equation(out,title,text,label,array):
    a=np.asarray(array,dtype=float)
    if a.ndim==1:a=a[:,None]
    latex=label+r'=\begin{bmatrix}'+r'\\'.join(' & '.join(f'{v:.10g}' for v in row) for row in a)+r'\end{bmatrix}'
    equation(out,title,text,latex,matrix=a.tolist(),matrix_label=label)


def region_rows(source,width,nodes):
    result=[]
    for line in str(source).splitlines():
        if not line.strip():continue
        parts=[s.strip() for s in line.split(';')]
        if len(parts)!=width:raise ValueError(f'Region rows need {width} semicolon-separated entries.')
        a,b=map(lambda s:number(expression(s)),parts[:2])
        if b<=a:raise ValueError('Region right endpoint must exceed its left endpoint.')
        result.append((a,b,*[expression(s) for s in parts[2:]]))
    if result:
        result.sort(key=lambda q:q[0])
        ends=[result[0][0]]+[q[1] for q in result]
        if not np.isclose(ends[0],nodes[0],rtol=0,atol=1e-12) or not np.isclose(ends[-1],nodes[-1],rtol=0,atol=1e-12):raise ValueError('Regions must cover the full mesh domain.')
        for a,b in zip(result,result[1:]):
            if abs(a[1]-b[0])>1e-12:raise ValueError('Regions must be contiguous without gaps or overlaps.')
        for end in ends:
            if not np.any(np.isclose(nodes,end,rtol=0,atol=1e-12)):raise ValueError('Every region interface must be a mesh node.')
    return result


def primitive_supported(g):
    # Bounded symbolic family, avoiding unbounded attempts on arbitrary inputs.
    # Failure to support an expression does not imply it has no antiderivative.
    if g.is_polynomial(x):return True
    if not g.atoms(sp.Function):return False
    return all(fn.func in (sp.sin,sp.cos,sp.exp) and sp.diff(fn.args[0],x,2)==0 for fn in g.atoms(sp.Function)) and not g.has(sp.Abs)


def integrate_entry(out,g,left,right,mode,title):
    g=sp.expand(g)
    l,r=sp.Rational(str(left)),sp.Rational(str(right))
    primitive=None
    if mode!='Gauss 2' and primitive_supported(g):
        primitive=sp.integrate(g,x)
        if primitive.has(sp.Integral):primitive=None
    if primitive is not None:
        exact=primitive.subs(x,r)-primitive.subs(x,l)
        value=number(exact)
        equation(out,title,'Analytical integration using the displayed antiderivative.',
                 r'\int_{'+sp.latex(l)+'}^{'+sp.latex(r)+'} '+sp.latex(g)+r'\,dx=\left['+sp.latex(primitive)+r'\right]_{'+sp.latex(l)+'}^{'+sp.latex(r)+'}='+f'{value:.10g}')
        return value
    if mode=='Analytical only':raise ValueError(f'{title}: outside the supported analytical family. Select Auto or Gauss 2.')
    if mode=='Auto':
        # Split additive terms after expansion: an unsupported reaction/load
        # must not force an analytically integrable diffusion term into Gauss.
        exact_terms=[];primitives=[];remaining=[]
        for term in sp.Add.make_args(g):
            anti=sp.integrate(term,x) if primitive_supported(term) else None
            if anti is not None and not anti.has(sp.Integral):
                exact_terms.append(term);primitives.append(anti)
            else:
                remaining.append(term)
        if exact_terms:
            exact_part=sp.Add(*exact_terms);anti=sp.Add(*primitives)
            exact_value=number(anti.subs(x,r)-anti.subs(x,l))
            equation(out,title+' — analytical terms',
                     'Auto integrates the supported additive terms analytically. Only the remaining terms use quadrature.',
                     r'\int_{'+sp.latex(l)+'}^{'+sp.latex(r)+'}('+sp.latex(exact_part)+r')\,dx=\left['+sp.latex(anti)+r'\right]_{'+sp.latex(l)+'}^{'+sp.latex(r)+'}='+f'{exact_value:.10g}')
            remainder=sp.Add(*remaining)
            gauss_value=integrate_entry(out,remainder,left,right,'Gauss 2',title+' — remaining terms') if remaining else 0.0
            value=exact_value+gauss_value
            equation(out,title+' — combined integral',
                     'Sum the analytical contribution and the two-point Gauss contribution. Gauss fallback reflects the supported symbolic family, not proof that no antiderivative exists.',
                     r'I=I_{\mathrm{analytical}}+I_{\mathrm{Gauss}}=('+f'{exact_value:.10g}'+r')+('+f'{gauss_value:.10g}'+r')='+f'{value:.10g}')
            return value
    jac=(r-l)/2;mapping=(r+l)/2+jac*xi
    mapped=g.subs(x,mapping)
    rows=[]
    for q in (-1/sp.sqrt(3),1/sp.sqrt(3)):
        xp=mapping.subs(xi,q);gv=number(g.subs(x,xp));contribution=number(jac)*gv
        rows.append([number(q),1.0,number(xp),number(jac),gv,contribution])
    value=sum(row[-1] for row in rows)
    table(out,title+' — Gauss points',['xi','weight','x','Jacobian','integrand','weighted contribution'],rows)
    equation(out,title,'Two-point Gauss quadrature; the accompanying table shows both weighted contributions. Auto fallback means the analytical evaluator did not handle this expression, not that no antiderivative exists.',
             r'\int_{'+sp.latex(l)+'}^{'+sp.latex(r)+'} '+sp.latex(g)+r'\,dx\approx '+sp.latex(jac)+r'\sum_{k=1}^{2}w_k\left('+sp.latex(mapped)+r'\right)_{\xi=\xi_k}='+f'{value:.10g}')
    return value


def calculate(inputs):
    p=inputs
    nodes=np.array(values(p.get('nodes','')))
    if not 2<=len(nodes)<=101 or np.any(np.diff(nodes)<=0):raise ValueError('Enter 2 to 101 strictly increasing mesh nodes.')
    n=len(nodes)
    mode=p.get('integration','Auto')
    if mode not in ('Auto','Gauss 2','Analytical only'):raise ValueError('Unknown integration method.')
    out=dict(fea=True,tables=[],charts=[],notes=[],datasets=[],equations=[])
    forcing=expression(p.get('f','0'))
    # The region table overrides the scalar forcing in both calculation and report.
    loads=region_rows(p.get('loads',''),3,nodes)
    form=p.get('form','Self-adjoint')
    multiplier=sp.Integer(1)
    if form=='Standard':
        a,b,c=[expression(p.get(k,v)) for k,v in [('a','1'),('b','0'),('c','0')]]
        if a==0:raise ValueError('The leading coefficient cannot be zero.')
        domain=sp.Interval(sp.Rational(str(nodes[0])),sp.Rational(str(nodes[-1])))
        if sp.calculus.util.continuous_domain(a,x,domain)!=domain or sp.solveset(a,x,domain=domain)!=sp.EmptySet:
            raise ValueError('Standard-form leading coefficient must be continuous and nonzero on the closed domain.')
        log_mu=sp.integrate(sp.cancel((b-sp.diff(a,x))/a),x)
        if log_mu.has(sp.Integral):raise ValueError('Unable to derive an explicit integrating factor; enter a self-adjoint equation instead.')
        mu=sp.simplify(sp.exp(log_mu))
        cp,cr=sp.simplify(mu*a),sp.simplify(-mu*c)
        multiplier=-mu
        rhs_latex=r'f_{\mathrm{in}}(x)' if loads else sp.latex(forcing)
        equation(out,'Original standard equation',"The entered equation is a*y'' + b*y' + c*y = f. "+('The piecewise forcing below overrides the scalar right-side field.' if loads else ''),'('+sp.latex(a)+r")y''+("+sp.latex(b)+r")y'+("+sp.latex(c)+')y='+rhs_latex)
        for j,(lo,hi,load) in enumerate(loads,1):
            upper=r'\leq' if j==len(loads) else '<'
            equation(out,f'Original piecewise forcing: region {j}',
                     'Original right side, before the integrating-factor transformation.',
                     r'f_{\mathrm{in}}(x)='+sp.latex(load)+r',\quad '+sp.latex(sp.Rational(str(lo)))+r'\leq x'+upper+sp.latex(sp.Rational(str(hi))))
        equation(out,'Integrating factor','Multiply the original equation by -mu to use the canonical negative-divergence convention. A nonzero constant scaling does not change the solution.',r'\mu=\exp\left(\int\frac{b-a\prime}{a}\,dx\right)='+sp.latex(mu))
    elif form=='Self-adjoint':cp,cr=expression(p.get('p','1')),expression(p.get('r','0'))
    else:raise ValueError('Unknown equation form.')
    cf=sp.simplify(multiplier*forcing)
    regions=region_rows(p.get('regions',''),5,nodes)
    if regions and (form!='Self-adjoint' or loads):raise ValueError('Coefficient regions require Self-adjoint form and cannot be combined with forcing regions.')
    if not regions:
        regions=[(lo,hi,cp,cr,sp.simplify(multiplier*f)) for lo,hi,f in loads] if loads else [(nodes[0],nodes[-1],cp,cr,cf)]
    equation(out,'Canonical equation','All assembly uses this sign convention, including transformed or piecewise coefficients.',r'-(p y\prime)\prime+r y=f')
    for j,(lo,hi,pp,rr,ff) in enumerate(regions,1):
        equation(out,f'Region {j}: coefficients',f'Coefficients on [{lo:g}, {hi:g}]. Interface nodes share one continuous solution value.',r'p='+sp.latex(pp)+r',\quad r='+sp.latex(rr)+r',\quad f='+sp.latex(ff))
    equation(out,'Weighted residual','Multiply by a test function that vanishes at prescribed-value boundaries and integrate.',r'\int_a^b v\left[-(p y\prime)\prime+r y-f\right]dx=0')
    equation(out,'Integration by parts','The outward flux is q = n*p*y prime; n=-1 at the left endpoint and +1 at the right endpoint.',r'\int_a^b(p v\prime y\prime+r v y)\,dx=\int_a^b v f\,dx+[v p y\prime]_a^b')
    equation(out,'Galerkin element equations','Use the same local linear functions for approximation and test functions.',r'K^{(e)}_{ij}=\int_{x_L}^{x_R}(p N_i\prime N_j\prime+r N_i N_j)\,dx,\quad F^{(e)}_i=\int_{x_L}^{x_R}N_i f\,dx')
    equation(out,'Local-to-global assembly','A_e selects the two global node values listed in each element connectivity row. Overlapping entries are summed.',r'K=\sum_e A_e^T K^{(e)} A_e,\quad F=\sum_e A_e^T F^{(e)},\quad y_h=\sum_j y_j\phi_j')
    equation(out,'Quadrature mapping','When used, Gauss 2 evaluates two points with unit weights.',r'x=\frac{x_L+x_R}{2}+\frac{h_e}{2}\xi,\quad J=\frac{h_e}{2},\quad \xi_{1,2}=\mp\frac{1}{\sqrt{3}},\quad w_{1,2}=1')
    K=np.zeros((n,n));F=np.zeros(n);elements=[];p_signs=set()
    for e,(lo,hi) in enumerate(zip(nodes[:-1],nodes[1:])):
        mid=(lo+hi)/2
        reg=next(q for q in regions if q[0]<=mid<=q[1])
        pp,rr,ff=reg[2:]
        # Exclude poles, zeros of the leading coefficient, and non-real domains.
        for coefficient in (pp,rr,ff):
            domain=sp.calculus.util.continuous_domain(coefficient,x,sp.Interval(sp.Rational(str(lo)),sp.Rational(str(hi))))
            if domain!=sp.Interval(sp.Rational(str(lo)),sp.Rational(str(hi))):raise ValueError(f'Element {e+1}: a coefficient is not continuous and real on the closed element. Put discontinuities at interfaces; singular equations are unsupported.')
        zeros=sp.solveset(pp,x,domain=sp.Interval(sp.Rational(str(lo)),sp.Rational(str(hi))))
        if zeros!=sp.EmptySet:raise ValueError(f'Element {e+1}: flux coefficient may vanish; regular nonzero coefficients are required.')
        p_signs.add(np.sign(number(pp.subs(x,mid))))
        l,r=sp.Rational(str(lo)),sp.Rational(str(hi));h=r-l
        shape=[(r-x)/h,(x-l)/h]
        equation(out,f'Element {e+1}: shapes and connectivity',f'Local nodes 1,2 map to global nodes {e+1},{e+2}; length {float(h):g}.',r'N_1='+sp.latex(shape[0])+r',\quad N_2='+sp.latex(shape[1])+r',\quad N_1\prime='+sp.latex(-1/h)+r',\quad N_2\prime='+sp.latex(1/h))
        ke=np.zeros((2,2));fe=np.zeros(2)
        for i in range(2):
            for j in range(i,2):
                ke[i,j]=integrate_entry(out,pp*sp.diff(shape[i],x)*sp.diff(shape[j],x)+rr*shape[i]*shape[j],lo,hi,mode,f'Element {e+1}: K{i+1}{j+1}')
                ke[j,i]=ke[i,j]
            fe[i]=integrate_entry(out,ff*shape[i],lo,hi,mode,f'Element {e+1}: F{i+1}')
        K[e:e+2,e:e+2]+=ke;F[e:e+2]+=fe
        elements.append(dict(nodes=[e+1,e+2],matrix=ke.tolist(),rhs=fe.tolist()))
        matrix_equation(out,f'Element {e+1}: matrix','Symmetry gives K21 = K12.',r'K^{('+str(e+1)+')}',ke)
        matrix_equation(out,f'Element {e+1}: load','Add the local entries at the corresponding global node IDs.',r'F^{('+str(e+1)+')}',fe)
        table(out,f'Element {e+1}: matrix and load',['local row','global node','K1','K2','F'],[[i+1,e+i+1,*ke[i].tolist(),fe[i]] for i in range(2)])
    if len(p_signs)>1:raise ValueError('Flux coefficient must have a consistent nonzero sign across all regions.')
    if -1 in p_signs:out['notes'].append('WARNING: Negative flux coefficient. The matrix is not assumed positive definite; the direct solver supports nonsingular indefinite systems.')
    table(out,'Mesh and connectivity',['element','left node','right node','left x','right x','length'],[[e+1,e+1,e+2,lo,hi,hi-lo] for e,(lo,hi) in enumerate(zip(nodes[:-1],nodes[1:]))])
    out['charts'].append(dict(title='Finite-element mesh (nodes numbered left to right)',xlabel='x',ylabel='Mesh',curves=[dict(label='Mesh nodes',x=nodes.tolist(),y=[0.0]*n,style='nodes')]))
    matrix_equation(out,'Global assembly before boundary terms','Shared-node contributions are added, not averaged.', 'K',K)
    matrix_equation(out,'Global load before boundary terms','Sum all element load vectors.', 'F',F)
    prescribed={};bcrows=[]
    for side,index,normal,reg in [('left',0,-1,regions[0]),('right',n-1,1,regions[-1])]:
        alpha,beta,gamma=[number(expression(p.get(side+'_'+key,default))) for key,default in [('alpha','1'),('beta','0'),('gamma','0')]]
        if alpha==beta==0:raise ValueError('Boundary alpha and beta cannot both be zero.')
        pv=number(reg[2].subs(x,sp.Rational(str(nodes[index]))))
        if beta==0:
            prescribed[index]=gamma/alpha;bcrows.append([side,index+1,'prescribed y',gamma/alpha,0,0])
            equation(out,side.title()+' essential boundary','This value is kept for manual partitioning; no row is silently eliminated.',r'y_{'+str(index+1)+'}='+f'{gamma/alpha:.12g}')
        else:
            dk=normal*pv*alpha/beta;df=normal*pv*gamma/beta
            K[index,index]+=dk;F[index]+=df
            bcrows.append([side,index+1,'natural / Robin',None,dk,df])
            equation(out,side.title()+' natural / Robin boundary',f'Outward normal {normal:+d}. Add {dk:.12g} to the endpoint diagonal and {df:.12g} to its load.',r'q=n p\frac{\gamma-\alpha y}{\beta},\quad K_{jj}\leftarrow K_{jj}+\frac{n p\alpha}{\beta},\quad F_j\leftarrow F_j+\frac{n p\gamma}{\beta}')
    for line in str(p.get('point_loads','')).splitlines():
        if not line.strip():continue
        pair=values(line)
        if len(pair)!=2 or pair[0]!=int(pair[0]) or not 1<=pair[0]<=n:raise ValueError('Point loads require a valid integer node ID and a load.')
        F[int(pair[0])-1]+=pair[1]
        out['notes'].append(f'Added nodal load {pair[1]:.12g} at node {int(pair[0])}.')
    free=[i for i in range(n) if i not in prescribed]
    table(out,'Boundary contributions',['side','node','type','prescribed y','diagonal addition','load addition'],bcrows)
    table(out,'Assembled system for manual partition',['node']+[f'K{j+1}' for j in range(n)]+['F','prescribed y'],[[i+1,*K[i].tolist(),F[i],prescribed.get(i)] for i in range(n)])
    matrix_equation(out,'Assembled matrix including boundary terms','Use this matrix when forming your reduced system.', 'K',K)
    matrix_equation(out,'Assembled right side including boundary and point loads','Prescribed-value contributions have not yet been subtracted.', 'F',F)
    equation(out,'Manual partition','Free nodes: '+', '.join(str(i+1) for i in free)+'. Prescribed nodes: '+(', '.join(str(i+1) for i in prescribed) or 'none')+'. Form the reduced matrix and right side yourself, then choose Check partition and solve.',r'K_{FF}y_F=F_F-K_{FD}y_D')
    out['fea_system']=dict(nodes=nodes.tolist(),matrix=K.tolist(),rhs=F.tolist(),prescribed={str(i+1):v for i,v in prescribed.items()},free_nodes=[i+1 for i in free],elements=elements)
    if p.get('workflow','Assemble only')=='Assemble only' or p.get('equations_only'):
        out['notes'].append('Assembly complete. No reduced system was generated or solved. Enter your manual partition to continue.')
        return out
    if p.get('workflow')!='Check partition and solve':raise ValueError('Unknown FEA workflow.')
    order=values(p.get('free_order',''))
    if any(v!=int(v) for v in order) or sorted(order)!=[i+1 for i in free]:raise ValueError('Free-node order must contain every free node ID exactly once and no prescribed nodes.')
    order=[int(v)-1 for v in order]
    y=np.zeros(n)
    for i,v in prescribed.items():y[i]=v
    if free:
        rows=[values(line) for line in str(p.get('manual_matrix','')).splitlines() if line.strip()]
        if len(rows)!=len(free) or any(len(row)!=len(free) for row in rows):raise ValueError('Enter a square reduced matrix with one row per free node.')
        entered=np.array(rows);rhs=np.array(values(p.get('manual_rhs','')))
        if len(rhs)!=len(free):raise ValueError('Enter one right-side value per free node.')
        tol=number(expression(p.get('partition_tolerance','1e-6')))
        if not 0<tol<=0.01:raise ValueError('Use an absolute partition-check tolerance greater than zero and at most 0.01.')
        expected=K[np.ix_(order,order)]
        expected_rhs=F[order]-K[order,:]@y
        deltaK=float(np.max(abs(entered-expected)));deltaF=float(np.max(abs(rhs-expected_rhs)))
        if max(deltaK,deltaF)>tol:raise ValueError(f'Manual partition does not match assembly: maximum matrix difference {deltaK:.8g}, RHS difference {deltaF:.8g}; tolerance {tol:g}. Check free-node order, boundary contributions and subtraction of K_FD*y_D.')
        condition=float(np.linalg.cond(entered))
        if not np.isfinite(condition) or condition>1/np.finfo(float).eps:raise ValueError('Reduced system is singular or numerically singular; check boundary conditions.')
        try:y[order]=np.linalg.solve(entered,rhs)
        except np.linalg.LinAlgError as error:raise ValueError('Reduced system is singular; check boundary conditions.') from error
        if condition>1e10:out['notes'].append('WARNING: Reduced matrix is poorly conditioned; displayed digits may not be reliable.')
        table(out,'Manual partition validation',['maximum K difference','maximum RHS difference','absolute tolerance','condition number'],[[deltaK,deltaF,tol,condition]])
        matrix_equation(out,'User-entered reduced matrix','Validated against assembly; the solution uses your entered values, including their rounding.',r'K_{FF}',entered)
        matrix_equation(out,'User-entered reduced right side','Validated against the prescribed-value elimination.',r'F_F-K_{FD}y_D',rhs)
    residual=K@y-F
    if not np.all(np.isfinite(y)):raise ValueError('Solution contains nonfinite values.')
    table(out,'FEA nodal solution',['node','x','y','assembled residual','prescribed'],[[i+1,nodes[i],y[i],residual[i],i in prescribed] for i in range(n)])
    slopes=np.diff(y)/np.diff(nodes)
    table(out,'Element slopes',['element','left x','right x','dy/dx'],[[i+1,nodes[i],nodes[i+1],slopes[i]] for i in range(n-1)])
    equation(out,'Reconstruction','The continuous solution is linear inside each element. Its slope is elementwise constant and need not agree across an interface.',r'y_h=N_1 y_L+N_2 y_R,\quad y_h\prime=\frac{y_R-y_L}{h_e}')
    out['notes'].append('Maximum free-node assembled residual: '+f'{max((abs(residual[i]) for i in free),default=0):.12g}'+'. Prescribed-node residuals are reactions, not solver errors. A small algebraic residual does not establish mesh accuracy.')
    out['charts'].append(dict(title='Finite-element solution',xlabel='x',ylabel='y',curves=[dict(label='Linear FEA',x=nodes.tolist(),y=y.tolist(),style='nodes')]))
    out['charts'].append(dict(title='Element slopes (jumps at shared nodes)',xlabel='x',ylabel="y prime",curves=[dict(label='Element slope',x=[v for pair in zip(nodes[:-1],nodes[1:]) for v in pair],y=[v for s in slopes for v in (s,s)],style='line')]))
    out['datasets'].append(dict(label='FEA nodal solution',x=nodes.tolist(),y=y.tolist()))
    out['fea_solution']=y.tolist()
    return out
