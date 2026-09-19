"""Study reports for polynomial least squares and geometric area."""
import numpy as np


def table(title, columns, rows):
    return {'title': title, 'columns': columns, 'rows': rows}


def append_working(out, poly, xs, ys, lower, upper, least_squares):
    degree = len(poly.coef)-1
    if least_squares:
        # Display the conventional equations; retain the stable scaled fit.
        design = np.vander(xs, degree+1, increasing=True)
        normal = design.T @ design
        rhs = design.T @ ys
        out['notes'].append(
            f'Least-squares working: f(x)=sum(a_j*x^j, j=0..{degree}). '
            'Minimize S=sum_i(f(x_i)-y_i)^2. Setting each partial derivative '
            'to zero gives sum_j a_j sum_i x_i^(k+j)=sum_i x_i^k*y_i, '
            f'k=0..{degree}; equivalently (A^T A)a=A^T y, A_ij=x_i^j. '
            'The matrices below document these equations; coefficients are computed '
            'by the scaled least-squares solver, not by inverting A^T A.')
        out['tables'].append(table('Least squares — design matrix A',
            ['x']+[f'x^{j}' for j in range(degree+1)]+['y'],
            [[float(x)]+list(map(float,row))+[float(y)] for x,row,y in zip(xs,design,ys)]))
        out['tables'].append(table('Least squares — normal equations (A^T A)a = A^T y',
            ['k']+[f'a{j} multiplier' for j in range(degree+1)]+['right side'],
            [[k]+list(map(float,row))+[float(rhs[k])] for k,row in enumerate(normal)]))
    # Isolate sign-changing real roots recursively between derivative roots.
    # Including near-zero critical points also captures tangential roots.
    def roots_between(p, a, b):
        if len(p.coef)<=1:
            return []
        critical=roots_between(p.deriv(),a,b)
        nodes=[a]+critical+[b]
        roots=[]
        scale=max(1.,sum(abs(p.coef)))
        for x in nodes:
            if abs(p(x))<1e-12*scale:
                roots.append(x)
        for left,right in zip(nodes,nodes[1:]):
            fl,fr=p(left),p(right)
            if fl*fr<0:
                for _ in range(80):
                    mid=(left+right)/2;fm=p(mid)
                    if fm==0:
                        left=right=mid;break
                    if fl*fm<0:right=mid
                    else:left=mid;fl=fm
                roots.append((left+right)/2)
        unique=[]
        for x in sorted(roots):
            if not unique or abs(x-unique[-1])>1e-10*max(1.,b-a):unique.append(x)
        return unique
    roots=roots_between(poly,lower,upper)
    cuts=[lower]+[r for r in roots if lower<r<upper]+[upper]
    anti=poly.integ()
    out['notes'].append('Area working: find zeros of f in the selected interval, '
        'split there, and sum |F(b)-F(a)| over the resulting regions, where F\'=f. '
        'Roots are numerical estimates isolated between derivative roots; '
        'near-multiple roots are subject to floating-point precision.')
    out['tables'].append(table('Antiderivative F — expanded coefficients (constant chosen by solver)',
        ['power of x','coefficient'],[[i,float(c)] for i,c in enumerate(anti.convert().coef)]))
    out['tables'].append(table('Polynomial roots in area interval',['x','f(x)'],
        [[float(r),float(poly(r))] for r in roots]))
    rows=[]
    for a,b in zip(cuts,cuts[1:]):
        value=float(anti(b)-anti(a))
        rows.append([float(a),float(b),'positive' if poly((a+b)/2)>0 else 'negative' if poly((a+b)/2)<0 else 'zero',float(anti(a)),float(anti(b)),value,abs(value)])
    out['tables'].append(table('Area — split integrals',
        ['lower x','upper x','sign inside','F(lower)','F(upper)','signed integral','region area'],rows))
    out['tables'].append(table('Total geometric area',['lower x','upper x','area'],
        [[lower,upper,sum(r[-1] for r in rows)]]))
    xp=np.unique(np.concatenate([np.linspace(lower,upper,600),cuts]))
    out['charts'].append({'title':'Polynomial area — split at roots','ylabel':'f(x)',
        'curves':[{'label':'Fitted polynomial','x':xp.tolist(),'y':poly(xp).tolist(),'style':'line'},
                  {'label':'x-axis','x':[lower,upper],'y':[0,0],'style':'line'},
                  {'label':'Roots','x':list(map(float,roots)),'y':[0]*len(roots),'style':'points'}]})
