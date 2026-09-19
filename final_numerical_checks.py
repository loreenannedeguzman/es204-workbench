from pathlib import Path
import ast
p=Path('numerical_methods_program.py')
s=p.read_text(encoding='utf-8')

def edit(name, change):
    global s
    t=ast.parse(s);n=next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name==name)
    lines=s.splitlines(keepends=True);a=sum(map(len,lines[:n.lineno-1]));b=sum(map(len,lines[:n.end_lineno]))
    s=s[:a]+change(s[a:b])+s[b:]

edit('bisection_method',lambda b:b.replace('"root": root if root is not None else x_new,','"root": root,'))
def fd(b):
    a=b.index('    # --------------------------------------------------------',b.index('    h = get_number('))
    return b[:a]+'''    if h <= 0:
        print("Error: h must be positive.")
        return

'''+b[a:]
edit('finite_difference_approximations',fd)
def trap(b):
    old='''                print(
                    "Error: x-values must have uniform spacing."
                )

                return'''
    assert old in b
    b=b.replace(old,'''                # A-2 also allows unequal intervals for the trapezoidal rule.
                if not (all(v > 0 for v in h_values) or all(v < 0 for v in h_values)):
                    print("Error: x-values must be strictly ordered.")
                    return
                areas = [h_values[i] * (y_values[i] + y_values[i+1]) / 2
                         for i in range(len(h_values))]
                integral = sum(areas)
                print("\\nUnequal-spacing trapezoids: A_i = h_i*(y_i + y_(i+1))/2")
                for i, area in enumerate(areas):
                    print(f"i = {i}, h_i = {h_values[i]:.10f}, A_i = {area:.10f}")
                print(f"Integral = {integral:.10f}")
                reference = input("Enter exact integral (ENTER if unavailable): ").strip()
                if reference:
                    exact_value = float(reference)
                    error = abs(integral - exact_value)
                    print(f"Absolute error = {error:.10f}")
                    if exact_value != 0:
                        print(f"Absolute percentage error = {100*error/abs(exact_value):.6f}%")
                    else:
                        print("Percentage error is undefined for an exact integral of zero.")
                return''')
    return b
edit('trapezoidal_rule',trap)
edit('plot_root',lambda b:b.replace('label=f"Root = {root:.6f}"','label=f"Estimate = {root:.6f}"'))
ast.parse(s)
p.write_bytes(s.replace('\n','\r\n').encode('utf-8'))
print('Final targeted input and unequal-spacing integration checks applied.')
