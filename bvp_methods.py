"""Linear BVPs and iterative matrix methods, with reproducible equation reports.

No lecture-specific solutions are embedded. All equations come from input data.
"""
import math
import numpy as np
import sympy as sp
from numerics import expr, number, integer, finite, report, curve, ode_step


def output():
    return dict(notes=[], tables=[], charts=[], datasets=[], equations=[], warnings=[])


def equation(out, title, text, latex):
    out['notes'].append(title+':\n'+text)
    out['equations'].append({'title':title,'text':text,'latex':latex})


def num(v):
    return sp.Float(float(v), 10)


def emit_system(out, A, b, names):
    out['tables'].append(report('Linear system — augmented matrix', ['row']+names+['right side'],
        [[i+1]+[float(v) for v in row]+[float(b[i])] for i,row in enumerate(A)]))
    symbols=[sp.Symbol(n) for n in names]
    equation(out,'Matrix form','A u = b; unknown order: '+', '.join(names),r'A\mathbf{u}=\mathbf{b}')
    if len(b)<=12:
        latex=sp.latex(sp.Matrix(A).applyfunc(num))+sp.latex(sp.Matrix(symbols))+'='+sp.latex(sp.Matrix(b).applyfunc(num))
        equation(out,'Assembled matrix equation','See the augmented-matrix table and component equations below.',latex)
    for i,(row,rhs) in enumerate(zip(A,b)):
        lhs=sp.Add(*[num(v)*s for v,s in zip(row,symbols) if v!=0])
        equation(out,f'Equation {i+1}',f'{lhs} = {rhs:.12g}',sp.latex(lhs)+'='+sp.latex(num(rhs)))


def iterative(A,b,p,out,names):
    method=p.get('solver','Direct')
    n=len(b)
    if not p.get('equations_only') and np.linalg.matrix_rank(A)<n:
        raise ValueError('The matrix is singular or numerically rank deficient. Boundary conditions may be inconsistent or may not determine a unique solution.')
    if method=='Direct':
        if p.get('equations_only'):return None,False
        x=np.linalg.solve(A,b)
        out['notes'].append('Direct linear solve. Residual = max |A u - b|.')
        return x,True
    if method not in ('Jacobi','Gauss-Seidel','SOR'):
        raise ValueError('Choose Direct, Jacobi, Gauss-Seidel, or SOR.')
    if np.any(np.diag(A)==0):raise ValueError('A diagonal coefficient is zero. Reorder equations/unknowns or use Direct; no silent reordering is performed.')
    initial=p.get('initial','').strip()
    x=np.zeros(n) if not initial else np.array([number(v) for v in initial.replace(',',' ').split()])
    if len(x)!=n:raise ValueError(f'Enter {n} initial values (one per matrix unknown), or leave blank for zeros.')
    tol=number(p.get('tolerance','1e-6'));limit=integer(p.get('limit','1000'))
    if tol<=0 or limit>20000:raise ValueError('Use a positive tolerance and at most 20000 iterations.')
    if n*limit>2000000:raise ValueError('Too many iteration table cells; reduce the matrix size or iteration limit.')
    omega=number(p.get('omega','1')) if method=='SOR' else 1.
    if not 0<omega<2:raise ValueError('SOR requires 0 < omega < 2; this alone does not guarantee convergence.')
    metric=p.get('criterion','Residual');comparison=p.get('comparison','<')
    if metric not in ('Residual','Absolute update','Relative update %'):raise ValueError('Unknown stopping criterion.')
    if comparison not in ('<','<='):raise ValueError('Choose < or <=.')
    dominant=bool(np.all(abs(np.diag(A))>np.sum(abs(A),axis=1)-abs(np.diag(A))))
    out['notes'].append(f'Strict row diagonal dominance: {dominant}. Absence does not prove divergence. SOR convergence also depends on omega.')
    for i in range(n):
        old=lambda j:sp.Symbol(names[j]+'^{(k)}')
        new=lambda j:sp.Symbol(names[j]+'^{(k+1)}')
        rhs=(num(b[i])-sum(num(A[i,j])*(new(j) if method!='Jacobi' and j<i else old(j)) for j in range(n) if j!=i))/num(A[i,i])
        if method=='SOR':rhs=(1-num(omega))*old(i)+num(omega)*rhs
        equation(out,f'{method} update for {names[i]}',f'{names[i]}(k+1) = {rhs}',sp.latex(new(i))+'='+sp.latex(rhs))
    if p.get('equations_only'):return None,False
    rows=[];ok=False
    for k in range(1,limit+1):
        old=x.copy()
        with np.errstate(over='ignore',invalid='ignore'):
            for i in range(n):
                used=old if method=='Jacobi' else x
                candidate=(b[i]-A[i]@used+A[i,i]*used[i])/A[i,i]
                x[i]=(1-omega)*old[i]+omega*candidate
        if not np.all(np.isfinite(x)) or np.max(abs(x))>1e150:
            x=old;out['warnings'].append('Iteration diverged/overflowed; last finite vector retained.');break
        update=float(max(abs(x-old)));res=float(max(abs(A@x-b)))
        relative=np.divide(abs(x-old),abs(x),out=np.full(n,np.inf),where=x!=0)*100
        relative[(x==0)&(old==0)]=0
        rel=float(max(relative));value={'Residual':res,'Absolute update':update,'Relative update %':rel}[metric]
        ok=value<tol if comparison=='<' else value<=tol
        rows.append([k]+x.tolist()+[update,rel if math.isfinite(rel) else None,res,ok])
        if ok:break
    out['tables'].append(report(method+' iterations',['k']+names+['max absolute update','max relative update %','max residual','criterion met'],rows))
    if rows:out['charts'].append({'title':method+' convergence','xlabel':'iteration','ylabel':'max residual','curves':[curve('Residual',[r[0] for r in rows],[r[-2] for r in rows],'nodes')]})
    out['notes'].append(f'Stopping: {metric} {comparison} {tol:g}. Relative updates are undefined for changed components whose new value is zero (N/A). Converged: {ok}.')
    if not ok:out['warnings'].append('NOT CONVERGED under the selected criterion. Inspect residuals and adjust the method, initial vector, or iteration limit.')
    return x,ok


def matrix(p):
    out=output()
    rows=[[number(v) for v in line.replace(',',' ').split()] for line in p.get('matrix','').splitlines() if line.strip()]
    n=len(rows)
    if n<1 or n>250 or any(len(r)!=n for r in rows):raise ValueError('Enter a square matrix A (1–250 rows), one row per line.')
    A=np.array(rows);b=np.array([number(v) for v in p.get('rhs','').replace(',',' ').split()])
    if len(b)!=n:raise ValueError('The right-side vector must contain one value per matrix row.')
    names=[f'u_{i+1}' for i in range(n)]
    emit_system(out,A,b,names)
    # Equation-only requests still generate the selected iteration formulas.
    if p.get('equations_only'):
        q=dict(p,initial='')
        iterative(A,b,q,out,names)
        out['tables']=[t for t in out['tables'] if not t['title'].endswith(' iterations')]
        out['charts']=[];out['warnings']=[]
        out['notes']=[v for v in out['notes'] if not v.startswith('Stopping:')]
        out['notes'].append('Equation preview only; no solution is claimed.')
        return out
    x,ok=iterative(A,b,p,out,names)
    out['tables'].append(report('Matrix solution',['unknown','value'],[[n,float(v)] for n,v in zip(names,x)]))
    out['tables'].append(report('Matrix checks',['converged','max residual'],[[ok,float(max(abs(A@x-b)))]]))
    return out


def grid(p):
    a,b=number(p['left']),number(p['right'])
    if a>=b:raise ValueError('The left endpoint must be less than the right endpoint.')
    mode=p.get('grid','Points n (lesson)');value=p['size']
    if mode=='Step h':
        h=number(value)
        if h<=0:raise ValueError('h must be positive.')
        count=(b-a)/h
        if not math.isclose(count,round(count),rel_tol=1e-9,abs_tol=1e-9):raise ValueError('h must divide the interval; no silent endpoint adjustment is made.')
        intervals=round(count)
    else:intervals=integer(value)-(1 if mode=='Points n (lesson)' else 0)
    if not 2<=intervals<=250:raise ValueError('Use 2–250 intervals (3–251 points).')
    return np.linspace(a,b,intervals+1)


def bvp(p):
    out=output();xs=grid(p);a,b=xs[0],xs[-1];h=xs[1]-xs[0];N=len(xs)
    coeff=[expr(p[k]) for k in ('a','b','c','d')]
    fs=[sp.lambdify(sp.Symbol('x'),v,'math') for v in coeff]
    cv=np.array([[finite(f(x)) for f in fs] for x in xs])
    if np.any(cv[:,0]==0):raise ValueError('The coefficient of y\'\' must be nonzero at every grid node.')
    bc=[]
    for side in ('l','r'):
        alpha,beta,gamma=[number(p[side+k]) for k in ('alpha','beta','gamma')]
        if alpha==beta==0:raise ValueError('Each boundary must specify a nonzero alpha or beta.')
        bc.append((alpha,beta,gamma))
    x=sp.Symbol('x');y=sp.Function('y')(x)
    lhs=coeff[0]*sp.diff(y,x,2)+coeff[1]*sp.diff(y,x)+coeff[2]*y
    equation(out,'Governing equation',f'{coeff[0]} y\'\' + ({coeff[1]}) y\' + ({coeff[2]}) y = {coeff[3]}',sp.latex(lhs)+'='+sp.latex(coeff[3]))
    for side,point,(al,be,ga) in zip(('Left','Right'),(a,b),bc):
        equation(out,side+' boundary',f'{al:g} y({point:g}) + {be:g} y\'({point:g}) = {ga:g}',rf'{al:g}y({point:g})+{be:g}y^\prime({point:g})={ga:g}')
    out['notes'].append(f'Grid: {N} points, {N-1} intervals, h={h:.12g}. Lesson indexing: i=1..{N}. Boundary derivatives are with respect to increasing x, not outward normals.')
    equation(out,'Grid equations',f'h=({b:g}-{a:g})/({N}-1)={h:.12g}; x_i={a:g}+(i-1)*{h:.12g}',
        rf'h=\frac{{b-a}}{{n-1}}={h:.12g},\quad x_i={a:g}+(i-1)h,\quad i=1,\ldots,{N}')
    if p.get('method')=='Shooting':return shooting(p,out,xs,coeff,fs,bc)
    return differences(p,out,xs,cv,coeff,bc)


def differences(p,out,xs,cv,coeff,bc):
    h=xs[1]-xs[0];N=len(xs);A=np.zeros((N,N));b=np.zeros(N)
    equation(out,'Central interior formulas',"y'_i ≈ (y_(i+1)-y_(i-1))/(2h); y''_i ≈ (y_(i-1)-2y_i+y_(i+1))/h^2",
        r"y'_i\approx\frac{y_{i+1}-y_{i-1}}{2h},\quad y''_i\approx\frac{y_{i-1}-2y_i+y_{i+1}}{h^2}")
    equation(out,'Interior equation multiplied by h^2',
        '(a_i-h*b_i/2)y_(i-1)+(-2*a_i+h^2*c_i)y_i+(a_i+h*b_i/2)y_(i+1)=h^2*d_i',
        r'(a_i-hb_i/2)y_{i-1}+(-2a_i+h^2c_i)y_i+(a_i+hb_i/2)y_{i+1}=h^2d_i')
    for i in range(1,N-1):
        aa,bb,cc,dd=cv[i];A[i,i-1:i+2]=[aa-h*bb/2,-2*aa+h*h*cc,aa+h*bb/2];b[i]=h*h*dd
    for index,(al,be,ga) in zip((0,N-1),bc):
        A[index,index]=al;b[index]=ga
        if index==0:A[index,:3]+=be*np.array([-3,4,-1])/(2*h)
        else:A[index,-3:]+=be*np.array([1,-4,3])/(2*h)
    equation(out,'Second-order boundary derivatives',
        "y'_1 ≈ (-3y_1+4y_2-y_3)/(2h); y'_n ≈ (3y_n-4y_(n-1)+y_(n-2))/(2h)",
        r"y'_1\approx\frac{-3y_1+4y_2-y_3}{2h},\quad y'_n\approx\frac{3y_n-4y_{n-1}+y_{n-2}}{2h}")
    # Eliminate the boundary values, retaining the interior matrix taught in class.
    edges=[0,N-1];inner=list(range(1,N-1));B=A[np.ix_(edges,edges)]
    if np.linalg.matrix_rank(B)==2:
        C=np.linalg.solve(B,A[np.ix_(edges,inner)]);d=np.linalg.solve(B,b[edges])
        M=A[np.ix_(inner,inner)]-A[np.ix_(inner,edges)]@C
        rhs=b[inner]-A[np.ix_(inner,edges)]@d
        names=[f'y_{i+1}' for i in inner]
        symbols=[sp.Symbol(n) for n in names]
        for j,i in enumerate(edges):
            val=num(d[j])-sum(num(C[j,k])*symbols[k] for k in range(len(inner)) if C[j,k]!=0)
            equation(out,f'Boundary elimination y_{i+1}',f'y_{i+1} = {val}',rf'y_{{{i+1}}}='+sp.latex(val))
        reconstruct=lambda u:np.concatenate(([float(d[0]-C[0]@u)],u,[float(d[1]-C[1]@u)]))
    else:
        M,rhs=A,b;names=[f'y_{i+1}' for i in range(N)];reconstruct=lambda u:u
        out['notes'].append('Boundary elimination is singular; using the full system including boundary unknowns.')
    emit_system(out,M,rhs,names)
    if p.get('equations_only'):
        if p.get('solver','Direct')!='Direct':
            q=dict(p,initial='');iterative(M,rhs,q,out,names)
            out['tables']=[t for t in out['tables'] if not t['title'].endswith(' iterations')];out['charts']=[];out['warnings']=[]
            out['notes']=[v for v in out['notes'] if not v.startswith('Stopping:')]
        out['notes'].append('Equation preview only; no solution is claimed.');return out
    u,ok=iterative(M,rhs,p,out,names);ys=reconstruct(u)
    zs=np.gradient(ys,h,edge_order=2)
    out['tables'].append(report('Discrete system checks',['converged','max full-system residual'],[[ok,float(max(abs(A@ys-b)))]]))
    return solution(p,out,xs,ys,zs,'Finite differences')


def shooting(p,out,xs,coeff,fs,bc):
    method=p.get('ivp','RK4');passes=integer(p.get('passes','1'));h=xs[1]-xs[0]
    if passes>100:raise ValueError('Use at most 100 corrector passes.')
    al,be,ga=bc[0];ar,br,gr=bc[1]
    start=lambda s:np.array([s,(ga-al*s)/be]) if be!=0 else np.array([ga/al,s])
    variable='y(left)' if be!=0 else "y'(left)"
    equation(out,'First-order system',f"y'=z; z'=({coeff[3]}-({coeff[1]})z-({coeff[2]})y)/({coeff[0]})",
        r"y'=z,\quad z'="+sp.latex((coeff[3]-coeff[1]*sp.Symbol('z')-coeff[2]*sp.Symbol('y'))/coeff[0]))
    equation(out,'Shooting parameter and residual',f's={variable}; R(s)={ar:g} y(right;s)+{br:g} z(right;s)-({gr:g}); solve R(s)=0.',
        rf'R(s)={ar:g}y(b;s)+{br:g}z(b;s)-({gr:g})=0')
    initial=(f'y(left)=s, z(left)=({ga:g}-({al:g})s)/({be:g})' if be else f'y(left)={ga/al:g}, z(left)=s')
    out['notes'].append('Trial initial conditions: '+initial)
    equation(out,'Secant update', 's_new = s_1 - R(s_1)*(s_1-s_0)/(R(s_1)-R(s_0))',r's_{k+1}=s_k-R(s_k)\frac{s_k-s_{k-1}}{R(s_k)-R(s_{k-1})}')
    step_formulas={
      'Euler':r'U_{i+1}=U_i+hK_1',
      'Predictor-Corrector':r'U^p=U_i+hF(x_i,U_i),\quad U^{(j+1)}=U_i+\frac h2[F(x_i,U_i)+F(x_{i+1},U^{(j)})]',
      'RK2.1 (Heun)':r'K_2=F(x_i+h,U_i+hK_1),\quad U_{i+1}=U_i+h(K_1+K_2)/2',
      'RK2.2 (lesson Ralston)':r'K_2=F(x_i+3h/4,U_i+3hK_1/4),\quad U_{i+1}=U_i+h(K_1/3+2K_2/3)',
      'RK3':r'K_2=F(x_i+h/2,U_i+hK_1/2),\quad K_3=F(x_i+h,U_i-hK_1+2hK_2),\quad U_{i+1}=U_i+h(K_1+4K_2+K_3)/6',
      'RK4':r'K_2=F(x_i+h/2,U_i+hK_1/2),\quad K_3=F(x_i+h/2,U_i+hK_2/2),\quad K_4=F(x_i+h,U_i+hK_3),\quad U_{i+1}=U_i+h(K_1+2K_2+2K_3+K_4)/6'}
    if method not in step_formulas:raise ValueError('Choose a supported IVP step method.')
    equation(out,method+' vector stepping (U=[y,z], K1=F(x_i,U_i))',step_formulas[method],step_formulas[method])
    aa,bb,cc,dd=coeff
    yi,zi,hh=sp.symbols('y_i z_i h');xx=sp.Symbol('x')
    def field_at(at,u):
        return sp.Matrix([u[1],(dd.subs(xx,at)-bb.subs(xx,at)*u[1]-cc.subs(xx,at)*u[0])/aa.subs(xx,at)])
    u=sp.Matrix([yi,zi]);k1=field_at(xx,u)
    equation(out,'Problem-specific first stage',str(k1),r'K_1='+sp.latex(k1))
    # Later stage equations use named vector components to remain readable.
    k1s=sp.Matrix(sp.symbols('k_1y k_1z'));k2s=sp.Matrix(sp.symbols('k_2y k_2z'));k3s=sp.Matrix(sp.symbols('k_3y k_3z'))
    stages=[]
    if method in ('RK3','RK4'):stages.append(('K_2',field_at(xx+hh/2,u+hh*k1s/2)))
    elif method=='RK2.1 (Heun)':stages.append(('K_2',field_at(xx+hh,u+hh*k1s)))
    elif method=='RK2.2 (lesson Ralston)':stages.append(('K_2',field_at(xx+3*hh/4,u+3*hh*k1s/4)))
    if method=='RK3':stages.append(('K_3',field_at(xx+hh,u-hh*k1s+2*hh*k2s)))
    if method=='RK4':stages.extend([('K_3',field_at(xx+hh/2,u+hh*k2s/2)),('K_4',field_at(xx+hh,u+hh*k3s))])
    for name,v in stages:equation(out,'Problem-specific '+name,str(v),name+'='+sp.latex(v))
    if method=='Predictor-Corrector':
        up=sp.Matrix(sp.symbols('y_c z_c'));v=u+hh/2*(k1+field_at(xx+hh,up))
        equation(out,'Problem-specific corrector',str(v),r'U^{(j+1)}='+sp.latex(v))
    if p.get('equations_only'):return out
    guesses=[number(p.get('guess0','0')),number(p.get('guess1','1'))]
    if guesses[0]==guesses[1]:raise ValueError('Use two distinct shooting guesses.')
    tol=number(p.get('shoot_tol','1e-8'));limit=integer(p.get('shoot_limit','20'))
    if tol<=0 or limit>200:raise ValueError('Use a positive shooting tolerance and at most 200 updates.')
    comparison=p.get('comparison','<');stop=p.get('shoot_stop','Boundary residual')
    if stop not in ('Boundary residual','Absolute parameter change'):raise ValueError('Unknown shooting stopping criterion.')
    def f(x,u):
        aa,bb,cc,dd=[finite(fun(x)) for fun in fs]
        if aa==0:raise ValueError("Coefficient of y'' is zero at an IVP stage.")
        return np.array([u[1],(dd-bb*u[1]-cc*u[0])/aa])
    def trial(s):
        uv=[start(s)]
        for x in xs[:-1]:
            u=ode_step(f,x,uv[-1],h,method,passes)
            if not np.all(np.isfinite(u)):raise ValueError('Shooting overflowed; change grid or guesses.')
            uv.append(u)
        uv=np.array(uv);res=float(ar*uv[-1,0]+br*uv[-1,1]-gr)
        return uv,res
    first_trials=[trial(s) for s in guesses]
    if abs(first_trials[1][1]-first_trials[0][1])<1e-14*max(1.,abs(first_trials[0][1]),abs(first_trials[1][1])):
        raise ValueError('The boundary residual is insensitive to the shooting parameter: the discrete BVP may have no unique solution. Check boundary conditions and the grid.')
    records=[];trials=[];ok=False
    for k in range(limit+2):
        if k>=2:
            s0,r0=records[-2][1],records[-2][4];s1,r1=records[-1][1],records[-1][4]
            den=r1-r0
            if abs(den)<1e-14*max(1.,abs(r0),abs(r1)):
                out['warnings'].append('Shooting residual has negligible sensitivity; solution may be nonunique or incompatible.');break
            guesses.append(s1-r1*(s1-s0)/den)
        s=guesses[k];uv,res=first_trials[k] if k<2 else trial(s);change=abs(s-records[-1][1]) if records else None
        metric=abs(res) if stop=='Boundary residual' else change
        ok=metric is not None and (metric<tol if comparison=='<' else metric<=tol)
        records.append([k,s,float(uv[-1,0]),float(uv[-1,1]),res,change,ok]);trials.append(uv)
        if ok:break
    out['tables'].append(report('Shooting trials',['trial','parameter s','y(right)','z(right)','boundary residual','parameter change','criterion met'],records))
    out['notes'].append(f'Shooting stop: {stop} {comparison} {tol:g}. Criterion met: {ok}. Check boundary residual separately when stopping on parameter change.')
    if not ok:out['warnings'].append('Shooting did not converge under the selected criterion.')
    for k,uvv in enumerate(trials):
        out['tables'].append(report(f'Shooting trajectory {k}', ['i','x','y','z'],[[i+1,float(x),float(u[0]),float(u[1])] for i,(x,u) in enumerate(zip(xs,uvv))]))
    out['charts'].append({'title':'Shooting trajectories','ylabel':'y','curves':[curve(f'Trial {k}',xs,u[:,0],'nodes') for k,u in enumerate(trials)]})
    uv=trials[-1]
    return solution(p,out,xs,uv[:,0],uv[:,1],'Shooting '+method)


def solution(p,out,xs,ys,zs,label):
    rows=[[i+1,float(x),float(y),float(z)] for i,(x,y,z) in enumerate(zip(xs,ys,zs))]
    out['tables'].append(report(label+' solution',['i','x','y','z = y prime (approx.)'],rows))
    curves=[curve(label,xs,ys,'nodes')]
    if p.get('exact','').strip():
        exact=expr(p['exact']);ef=sp.lambdify(sp.Symbol('x'),exact,'math');df=sp.lambdify(sp.Symbol('x'),sp.diff(exact,sp.Symbol('x')),'math')
        errors=[]
        for x,y,z in zip(xs,ys,zs):
            ref=finite(ef(x));dref=finite(df(x))
            errors.append([float(x),float(y),ref,(ref-y)/ref*100 if ref else None,float(z),dref,(dref-z)/dref*100 if dref else None])
        out['tables'].append(report('Exact-reference comparison',['x','y numerical','y reference','signed error %','z numerical','z reference','signed derivative error %'],errors))
        curves.append(curve('Given exact formula',xs,[ef(x) for x in xs]))
    out['charts'].append({'title':'BVP solution','ylabel':'y','curves':curves})
    out['charts'].append({'title':'BVP derivative','ylabel':"y'",'curves':[curve(label,xs,zs,'nodes')]})
    out['datasets'].append({'label':label,'x':list(map(float,xs)),'y':list(map(float,ys))})
    return out
