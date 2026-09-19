from pathlib import Path
import ast

p=Path('numerical_methods_program.py')
s=p.read_text(encoding='utf-8')

def edit(name, transform):
    global s
    tree=ast.parse(s)
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
    lines=s.splitlines(keepends=True)
    a=sum(map(len,lines[:node.lineno-1]));b=sum(map(len,lines[:node.end_lineno]))
    s=s[:a]+transform(s[a:b])+s[b:]

s=s.replace('import sympy as sp\n','''import sympy as sp

# Options 27-33 save results for plotting during this program session.
# Numerical formulas remain inside their own method blocks.
saved_results = {}
''',1)

for name,key in [('vandermonde_approximation','vandermonde'),('lagrangian_polynomial','lagrangian')]:
    def interpolation(block):
        block=block.replace('(show_plot=False):','():',1)
        block=block.replace('    print(',f'    saved_results.pop("{key}", None)\n    print(',1)
        a=block.index('    if show_plot:')
        block=block[:a]+f'''    # Retain the full-precision polynomial, not the rounded printed formula.
    saved_results["{key}"] = {{
        "x": list(x_values), "y": list(y_values), "polynomial": polynomial,
        "function": sp.lambdify(x, polynomial, "math")
    }}
'''
        return block
    edit(name,interpolation)

def fd(block):
    block=block.replace('    print(', '    saved_results.pop("finite_difference", None)\n    print(',1)
    a=block.index('    plot_choice = input(')
    # Calculation block stores values; option 36 handles plotting separately.
    block=block[:a]+'''    saved_results["finite_difference"] = {
        "x": list(x_values), "y": list(y_values),
        "first": list(fd_first), "second": list(fd_second),
        "lagrange_first": fprime, "lagrange_second": fdoubleprime
    }
    print("Results saved. Use option 36 to plot the derivative comparisons.")
'''
    return block
edit('lagrangian_finite_difference',fd)

def least(block):
    block=block.replace('    print(', '    saved_results.pop("least_squares", None)\n    print(',1)
    a=block.index('        x_min = np.min(x_data)')
    block=block[:a]+'''    if results:
        saved_results["least_squares"] = {
            "x": x_data.copy(), "y": y_data.copy(), "models": results
        }
        print("Results saved. Use option 37 to plot the fitted models.")
'''
    return block
edit('least_squares',least)

for name,key in [('bisection_method','bisection'),('newton_raphson_method','newton'),('secant_method','secant')]:
    def roots(block):
        block=block.replace('    print(',f'    saved_results.pop("root", None)\n    saved_results.pop("{key}", None)\n    print(',1)
        # Remove the extra pause: the main menu already asks whether to continue.
        block='\n'.join(l for l in block.split('\n') if 'input("\\nPress ENTER to return to the menu...")' not in l)
        if key=='bisection':
            for endpoint in ['a','b']:
                needle=f'        print(f"\\nExact root found at {endpoint} = {{{endpoint}:.10f}}")'
                block=block.replace(needle,needle+f'''
        result = {{"method": "Bisection", "function": f, "root": {endpoint},
                  "iterations": [{endpoint}], "converged": True}}
        saved_results["bisection"] = saved_results["root"] = result''')
            estimates='[row["c"] for row in rows]'
            converged='rows[-1]["converged"]'
            method='Bisection'
        elif key=='newton':
            # Save the actual starting point as well as each updated estimate.
            block=block.replace('    rows = []','    initial_x = x_current\n    rows = []',1)
            estimates='[initial_x] + [row["x_new"] for row in rows if "x_new" in row]'
            converged='rows[-1]["error_percent"] < tolerance'
            method='Newton-Raphson'
        else:
            block=block.replace('    root = None','    initial_points = [x_previous, x_current]\n    root = None',1)
            estimates='initial_points + [row["x"] for row in rows]'
            converged='root is not None'
            method='Secant'
        block+=f'''
    result = {{
        "method": "{method}", "function": f,
        "root": root if root is not None else x_new,
        "iterations": {estimates}, "converged": {converged}
    }}
    saved_results["{key}"] = saved_results["root"] = result
    print("Results saved. Options 38-39 use this root-finding run.")
'''
        return block
    edit(name,roots)

edit('plot_root_finding',lambda block:'''def plot_root_finding():
    result = saved_results.get("root")
    if result is None:
        print("Run option 31, 32, or 33 first to calculate a root.")
        return
    estimates = result["iterations"]
    left, right = min(estimates), max(estimates)
    margin = max((right - left) * 0.2, 0.5)
    status = "converged" if result["converged"] else "last estimate; not converged"
    plot_root(result["function"], result["root"], left-margin, right+margin,
              title=f"{result['method']}: {status}")
''')
edit('plot_root_convergence',lambda block:'''def plot_root_convergence():
    result = saved_results.get("root")
    if result is None:
        print("Run option 31, 32, or 33 first to calculate root iterations.")
        return
    plot_root_iterations(result["iterations"], title=result["method"] + " iterations")
''')

pos=s.index('def plot_root_finding():')
s=s[:pos]+'''def plot_vandermonde():
    result = saved_results.get("vandermonde")
    if result is None:
        print("Run option 27 first to calculate the Vandermonde polynomial.")
        return
    plot_data_and_fit(result["x"], result["y"], result["function"],
                      title="Vandermonde approximation")


def plot_lagrangian():
    result = saved_results.get("lagrangian")
    if result is None:
        print("Run option 28 first to calculate the Lagrangian polynomial.")
        return
    plot_data_and_fit(result["x"], result["y"], result["function"],
                      title="Lagrangian approximation")


def plot_finite_difference():
    result = saved_results.get("finite_difference")
    if result is None:
        print("Run option 29 first to calculate the derivative comparison.")
        return
    # Both sets of derivatives were calculated in option 29.
    xs = result["x"]
    dense_x = np.linspace(min(xs), max(xs), 500)
    for order, label in [("first", "First"), ("second", "Second")]:
        plt.figure()
        derivative = result["lagrange_" + order]
        plt.plot(dense_x, [derivative(x) for x in dense_x], label="Lagrangian")
        valid = [(x, y) for x, y in zip(xs, result[order]) if y is not None]
        if valid:
            plt.plot([v[0] for v in valid], [v[1] for v in valid], "o--", label="Finite difference")
        plt.xlabel("x")
        plt.ylabel(label + " derivative")
        plt.title(label + " derivative comparison")
        plt.legend()
        plt.grid(True)
        plt.show()


def plot_least_squares():
    result = saved_results.get("least_squares")
    if result is None:
        print("Run option 30 first to calculate the least-squares models.")
        return
    xs = np.linspace(min(result["x"]), max(result["x"]), 500)
    plt.figure(figsize=(9, 6))
    plt.scatter(result["x"], result["y"], label="Data points")
    for model in result["models"]:
        plt.plot(xs, [model["function"](x) for x in xs], label=model["name"])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Least-squares approximation")
    plt.legend()
    plt.grid(True)
    plt.show()


'''+s[pos:]
s=s.replace('            vandermonde_approximation(show_plot=True)','            plot_vandermonde()')
s=s.replace('            lagrangian_polynomial(show_plot=True)','            plot_lagrangian()')
s=s.replace('elif choice == "36":\n            lagrangian_finite_difference()', 'elif choice == "36":\n            plot_finite_difference()')
s=s.replace('elif choice == "37":\n            least_squares()', 'elif choice == "37":\n            plot_least_squares()')
ast.parse(s)
p.write_bytes(s.replace('\n','\r\n').encode('utf-8'))
print('Options 27-33 cache calculations; options 34-39 consume cached results.')
