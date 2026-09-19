from pathlib import Path
import ast
import re

source = Path(r'F:\ES204\notebooks\numerical_methods_program.py')
original = source.read_bytes()
Path('numerical_methods_program.before_validation.py').write_bytes(original)
s = original.decode('utf-8').replace('\r\n', '\n')

def change(old, new, count=1):
    global s
    assert s.count(old) == count, (old[:100], s.count(old), count)
    s = s.replace(old, new)

change('        f"{f\\"f\'(x^k)\\":>18} "', '        "f\'(x^k)".rjust(18) + " "')
# The concatenated header needs an explicit + on both sides of this expression.
change('        f"{\'f(x^k)\':>18} "\n        "f\'(x^k)".rjust(18) + " "\n        f"{\'x^(k+1)\':>16} "',
       '        f"{\'f(x^k)\':>18} "\n        + "f\'(x^k)".rjust(18) + " "\n        + f"{\'x^(k+1)\':>16} "')
change('format_formula(', 'formula_string(', 4)
change('proceeded_to_menu()', 'proceed_to_menu()', 4)
change('def plot_error_comparison(\n    error_results,', 'def plot_error_data(\n    error_results,')
change('import math\n', 'import math\nimport numpy as np\n')
change('import numpy as np\nimport matplotlib.pyplot as plt\n', '')

# Keep each solver independent. k values here are slopes; the notes use h*k.
a = s.index('# RK2.2 — MIDPOINT METHOD')
b = s.index('# RK3 —', a)
block = s[a:b].replace('MIDPOINT METHOD', "RALSTON METHOD (LESSON RK2.2)")
block = block.replace('x_i+h/2', 'x_i+3*h/4').replace('y_i+h*k1/2', 'y_i+3*h*k1/4')
block = block.replace('y_i + h*k2', 'y_i + h*(k1/3 + 2*k2/3)')
block = block.replace('x + h / 2', 'x + 3 * h / 4').replace('y + h * k1 / 2', 'y + 3 * h * k1 / 4')
block = block.replace('h * k2', 'h * (k1 / 3 + 2 * k2 / 3)')
s = s[:a] + block + s[b:]

# Reject invalid meshes instead of dividing by zero or silently changing h.
needle = '    n = round(\n        (xf - x0) / h\n    )'
change(needle, '''    if not math.isfinite(h) or h <= 0 or xf <= x0:
        raise ValueError("Require finite h > 0 and xf > x0.")

    n = round((xf - x0) / h)
    if n < 1 or not math.isclose(n * h, xf - x0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("h must divide the interval: use h = (xf - x0) / n.")''', 6)

# A zero reference cannot establish a percentage accuracy requirement.
change('    if exact_at_xf is None:\n\n        return None', '''    if exact_at_xf is None:

        return None
    if exact_at_xf == 0:
        return float("inf")  # Relative percentage error is undefined.''')
change('            error = 0.0', '            error = 0.0 if yi == 0 else float("inf")')
change('                error = 0\n', '                error = 0.0 if y == 0 else float("inf")\n')
change('''                        change = (
                            abs(
                                fine_value
                                -
                                coarse_value
                            )
                            *
                            100
                        )''', '                        change = float("inf")')

# Replace only the incorrect h_max block with the lesson's quadratic calculation.
a = s.index('def lagrange_hmax(')
b = s.index('\n# ============================================================', a)
s = s[:a] + '''def lagrange_hmax(method, f, exact, x0, y0, xf, initial_h=None, tolerance=1.0):
    """Lesson A-2 Example 2: interpolate p(m), then solve p(m)=tolerance.

    p(m) is maximum absolute percentage error, and h_m = h0 / 2**m.
    The three lesson samples are m = 6, 7, 8; the estimate is local.
    """
    if exact is None:
        print("Exact solution is required for Lagrange h_max.")
        return
    if initial_h is None:
        initial_h = h0
    if tolerance <= 0:
        raise ValueError("Tolerance must be positive.")

    m_values = [6, 7, 8]
    p_values = []
    print("\\nLAGRANGE INTERPOLATION FOR h_max")
    for m in m_values:
        h = initial_h / (2 ** m)
        xs, ys, actual_h = method(f, x0, y0, xf, h)
        errors = []
        for xi, yi in zip(xs, ys):
            reference = exact(xi)
            if reference == 0:
                if yi != 0:
                    print("Percentage error is undefined at a zero exact value.")
                    return
                continue  # The exactly imposed initial condition needs no test.
            errors.append(abs((reference - yi) / reference) * 100)
        if not errors:
            print("No nonzero exact values are available for percentage errors.")
            return
        p_values.append(max(errors))
        print(f"m = {m}, h = {actual_h:.10f}, p(m) = {p_values[-1]:.10f}%")

    # Same three equations a + b*m + c*m**2 = p(m) as the notes.
    A = np.array([[1, m, m**2] for m in m_values], dtype=float)
    a, b, c = np.linalg.solve(A, np.array(p_values))
    print(f"p(m) = {a:.10f} + ({b:.10f})m + ({c:.10f})m^2")
    roots = np.roots([c, b, a - tolerance])
    candidates = [r.real for r in roots if abs(r.imag) < 1e-9
                  and 6 <= r.real <= 8 and b + 2*c*r.real < 0]
    if not candidates:
        print("No decreasing threshold crossing lies within m = 6..8.")
        print("Use m-refinement to locate the threshold for this problem.")
        return
    m_estimate = min(candidates)
    h_estimate = initial_h / (2 ** m_estimate)
    # Equal spacing requires integer n. Verify the rounded-up choice directly.
    n = math.ceil((xf - x0) / h_estimate)
    h_used = (xf - x0) / n
    xs, ys, _ = method(f, x0, y0, xf, h_used)
    checked_error = max(abs((exact(xi) - yi) / exact(xi)) * 100
                        for xi, yi in zip(xs, ys) if exact(xi) != 0)
    print(f"Estimated m = {m_estimate:.10f}, estimated h_max = {h_estimate:.10f}")
    print(f"Equal-grid choice: n = {n}, h = {h_used:.10f}")
    print(f"Verified maximum error = {checked_error:.10f}%")
    print("Criterion satisfied." if checked_error <= tolerance else
          "Interpolation is only an estimate: refine further to meet the tolerance.")
    return m_estimate, h_estimate, checked_error

''' + s[b:]

# Avoid teaching an unconditional minimum guarantee from a binary search.
change('    n = 1\n\n    while True:', '''    print("This search assumes final-point error decreases as n increases.")
    n = 1

    while True:''')
change('f"is the minimum."', 'f"is the minimum under the decreasing-error assumption."')
for old,new in [(22,17),(23,18),(24,19),(25,20),(26,21)]:
    s = s.replace(f'    OPTION {old}\n', f'    OPTION {new}\n')

# Guard repeated nodes and too-short Simpson tables.
change('\n    h = h_values[0]\n', '''
    h = h_values[0]
    if h == 0:
        print("Error: x-values must be distinct.")
        return
''', 3)
change('\n        h = h_values[0]\n', '''
        h = h_values[0]
        if h == 0:
            print("Error: x-values must be distinct.")
            return
''', 2)
change('        n = len(x_values) - 1\n\n        if n % 2 != 0:', '''        n = len(x_values) - 1
        if n < 2:
            print("Error: Simpson's rule requires at least three data points.")
            return

        if n % 2 != 0:''')
for func in ['vandermonde_approximation','lagrangian_polynomial','lagrangian_finite_difference']:
    a=s.index('def '+func+'('); b=s.index('\ndef ',a+5)
    block=s[a:b]
    anchor='        y_values.append(y_i)\n'
    block=block.replace(anchor,anchor+'''
    if n < 1 or len(set(x_values)) != n:
        print("At least one point and distinct x-values are required.")
        return
''')
    s=s[:a]+block+s[b:]

# The menu's plot entries reuse the existing calculations instead of duplicating them.
change('def vandermonde_approximation():', 'def vandermonde_approximation(show_plot=False):')
change('def lagrangian_polynomial():', 'def lagrangian_polynomial(show_plot=False):')
for func in ['vandermonde_approximation','lagrangian_polynomial']:
    a=s.index('def '+func+'('); b=s.index('\ndef ',a+5)
    end=s.rfind('\n\n#',a,b)
    addition='''
    if show_plot:
        polynomial_function = sp.lambdify(x, polynomial, "math")
        plot_data_and_fit(x_values, y_values, polynomial_function)

'''
    s=s[:end]+addition+s[end:]
change('            plot_vandermonde()', '            vandermonde_approximation(show_plot=True)')
change('            plot_lagrangian()', '            lagrangian_polynomial(show_plot=True)')
change('            plot_finite_difference()', '            lagrangian_finite_difference()')
change('            plot_least_squares()', '            least_squares()')

# Small input adapters for the existing generic plotting functions; no solver duplication.
a=s.index('def menu():')
s=s[:a]+'''def plot_root_finding():
    expression = convert_expression(input("Enter f(x): "))
    function = sp.lambdify(x_symbol, expression, "math")
    root = get_number("Enter the root estimate from your calculation: ")
    left = get_number("Plot from x = ")
    right = get_number("Plot to x = ")
    print(f"Residual f(root) = {function(root):.10g}")
    plot_root(function, root, left, right, title="Root estimate and function")


def plot_root_convergence():
    values = [float(v) for v in input("Enter successive root estimates, comma-separated: ").split(",")]
    plot_root_iterations(values)


''' + s[a:]

# Secant can stop before making its first update (equal function values).
change('    f_current = f(x_current)\n\n    #', '    f_current = f(x_current)\n    x_new, f_new = x_current, f_current\n\n    #')
change('The requested tolerance was not achieved.', 'The requested tolerance was not achieved.', 3)

ast.parse(s)
Path('numerical_methods_program.py').write_bytes(s.replace('\n','\r\n').encode('utf-8'))
print('Staged corrections; Python syntax valid.')

