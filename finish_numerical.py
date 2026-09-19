from pathlib import Path
import ast
import re
p=Path('numerical_methods_program.py')
s=p.read_text(encoding='utf-8')

def block_edit(name, edit):
    global s
    a=s.index('def '+name+'(')
    b=s.find('\ndef ',a+5)
    if b<0: b=len(s)
    s=s[:a]+edit(s[a:b])+s[b:]

def pc(block):
    block=block.replace('    h\n):','    h,\n    corrector_iterations=1\n):',1)
    a=block.index('        # Corrector')
    b=block.index('        x = x + h',a)
    block=block[:a]+'''        # Keep y_i fixed while iterating the corrector (A-3, Example 10).
        y_corrected = y_predictor
        for k in range(corrector_iterations):
            y_corrected = y + (h / 2) * (f(x, y) + f(x + h, y_corrected))
        y = y_corrected

'''+block[b:]
    block=block.replace('    if not math.isfinite(h)', '''    if corrector_iterations < 1 or int(corrector_iterations) != corrector_iterations:
        raise ValueError("The number of corrector iterations must be a positive integer.")
    if not math.isfinite(h)''')
    return block
block_edit('predictor_corrector',pc)
def display(block):
    block=block.replace('    h\n):','    h,\n    corrector_iterations=1\n):',1)
    a=block.index('    x_values, y_values, actual_h = method(')
    b=block.index('    table =',a)
    block=block[:a]+'''    if method_name == "PC":
        x_values, y_values, actual_h = predictor_corrector(
            f, x0, y0, xf, h, corrector_iterations
        )
        print(f"Corrector iterations per step: {corrector_iterations}")
    else:
        x_values, y_values, actual_h = method(f, x0, y0, xf, h)

'''+block[b:]
    return block
block_edit('display_solution',display)
s=s.replace('''            display_solution(
                "PC",
                h
            )''','''            iterations = int(get_number("Corrector iterations per step (default 1): ", 1))
            display_solution("PC", h, corrector_iterations=iterations)''')

# Percentage errors at zero are undefined; retain the absolute error output.
for name in ['trapezoidal_rule','simpson_one_third_rule','compare_integration_methods']:
    def zero_errors(block):
        return re.sub(r'(error_percent = \([\s\S]*?\) \* 100)',
                      lambda m:m[0]+' if exact_value != 0 else float("nan")',block)
    block_edit(name,zero_errors)

# Independent root routines retain the lesson's stopping criteria.
for name in ['bisection_method','newton_raphson_method','secant_method']:
    def root_guard(block):
        needle='    # --------------------------------------------------------\n    # 5.'
        # Place before the iteration display, after tolerance and max_iterations.
        start=block.index('    max_iterations =')
        pos=block.index('    # --------------------------------------------------------',start)
        block=block[:pos]+'''    if tolerance <= 0 or max_iterations < 1:
        print("Tolerance and maximum iterations must be positive.")
        return

'''+block[pos:]
        return block
    block_edit(name,root_guard)

# Recognize an already exact Newton root before testing its derivative.
def newton(block):
    block=block.replace('        dfx = df(x_current)', '''        if fx == 0:
            root = x_current
            rows.append({"error_percent": 0.0})
            break
        dfx = df(x_current)''')
    block=block.replace('            error_percent = abs(x_new - x_current) * 100',
                        '            error_percent = 0.0 if f(x_new) == 0 else float("inf")')
    return block
block_edit('newton_raphson_method',newton)
def secant(block):
    block=block.replace('    for k in range(1, max_iterations + 1):', '''    for k in range(1, max_iterations + 1):
        if abs(f_current) < tolerance:
            root = x_current
            break
        if abs(f_previous) < tolerance:
            root = x_previous
            break
''')
    block=block.replace('            break\n\n        #', '            return\n\n        #',1) if False else block
    block=block.replace('            print("The Secant Method cannot continue.")\n            break',
                        '            print("The Secant Method cannot continue.")\n            return')
    return block
block_edit('secant_method',secant)

# Label the retained rounded reference explicitly instead of changing the user's IVP.
s=s.replace('exact_input = "-3*exp(x)*(x - 2.16667)"\n',
'''# The lecture rounds 13/6 to 2.16667. For an exact IVP reference, use
# exact_input = "exp(x)*(6.5 - 3*x)"; the rounded version reproduces its tables.
exact_input = "-3*exp(x)*(x - 2.16667)"
''',1)
# Above can match a commented example; fix placement explicitly if needed.
ast.parse(s)
p.write_bytes(s.replace('\n','\r\n').encode('utf-8'))
print('Additional focused guards and optional PC corrector iteration added.')
