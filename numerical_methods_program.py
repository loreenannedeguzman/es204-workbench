import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sympy as sp

# Options 27-33 save results for plotting during this program session.
# Numerical formulas remain inside their own method blocks.
saved_results = {}


# ============================================================
# NUMERICAL METHODS PROGRAM
# ============================================================
#
# This is a menu-driven numerical methods program.
#
# To solve a NEW problem, normally change ONLY the
# PROBLEM DEFINITION section.
#
# Enter:
#   1. x0
#   2. xf
#   3. y0
#   4. h0
#   5. f_input
#   6. exact_input OR exact_at_xf_given
#
# The ODE MUST be rearranged into:
#
#                   y' = f(x,y)
#
# ============================================================


# ============================================================
# PROBLEM DEFINITION
# ============================================================


# ------------------------------------------------------------
# GIVEN IVP
# ------------------------------------------------------------

x0 = 0
xf = 2
y0 = 1


# ------------------------------------------------------------
# INITIAL STEP SIZE
# ------------------------------------------------------------

h0 = .4


# ------------------------------------------------------------
# DIFFERENTIAL EQUATION
# ------------------------------------------------------------
#
# Example 3:
#
# Original:
#
#       y' = y - 3e^x
#
# Therefore enter:
#
#       f_input = "y - 3*exp(x)"
#
# ------------------------------------------------------------

# Rearrange y' + y**2 + exp(-x) - x = 0 before entering f(x,y).
f_input = "x - y**2 - exp(-x)"

# ------------------------------------------------------------
# EXACT SOLUTION
# ------------------------------------------------------------
#
# If the exact solution formula is available, enter it here.
#
# Example:
#
#       y(x) = -3e^x(x - 2.16667)
#
# Enter:
#
#       exact_input = "-3*exp(x)*(x - 2.16667)"
#
#
# If NO exact solution formula is available:
#
#       exact_input = ""
#
# Then enter the known exact value at xf below.
#
# ------------------------------------------------------------

# The lecture rounds 13/6 to 2.16667. For an exact IVP reference, use
# exact_input = "exp(x)*(6.5 - 3*x)"; the rounded version reproduces its tables.
exact_input = "" 


# ------------------------------------------------------------
# GIVEN EXACT VALUE AT xf
# ------------------------------------------------------------
#
# Only used when exact_input = ""
#
# Example:
#
#       exact_input = ""
#       exact_at_xf_given = -2.481808
#
# ------------------------------------------------------------

exact_at_xf_given = None


# ============================================================
# END OF PROBLEM DEFINITION
# ============================================================


# ============================================================
# SYMBOLIC SETUP
# ============================================================

x_symbol, y_symbol = sp.symbols("x y")


symbol_dictionary = {

    "x": x_symbol,
    "y": y_symbol,

    "pi": sp.pi,

    "sin": sp.sin,
    "cos": sp.cos,
    "tan": sp.tan,

    "asin": sp.asin,
    "acos": sp.acos,
    "atan": sp.atan,

    "sqrt": sp.sqrt,

    "exp": sp.exp,

    "log": sp.log,

    "ln": sp.log

}


# ============================================================
# FORMULA CONVERSION
# ============================================================


def convert_expression(expression):

    """
    Converts a text expression into a SymPy expression.

    Example:
        "y - 3*exp(x)"

    becomes:
        y - 3*exp(x)
    """

    return sp.sympify(
        expression.replace("^", "**"),
        locals=symbol_dictionary
    )


# ------------------------------------------------------------
# Differential equation
# ------------------------------------------------------------

# Only the right-hand side of y' = f(x,y) belongs in f_input.
try:
    f_expression = convert_expression(f_input)
    if not isinstance(f_expression, sp.Expr):
        raise ValueError("Enter an expression, not a complete equation.")
    unknown_symbols = f_expression.free_symbols - {x_symbol, y_symbol}
    if unknown_symbols:
        names = ", ".join(sorted(str(symbol) for symbol in unknown_symbols))
        raise ValueError(f"Unknown symbol(s): {names}. Only x and y are allowed.")
except (sp.SympifyError, ValueError, TypeError) as error:
    raise SystemExit(
        "Invalid f_input: " + str(error) + "\n"
        "First rearrange the ODE into y' = f(x,y), then enter only f(x,y).\n"
        "For y' + y**2 + exp(-x) - x = 0, use: x - y**2 - exp(-x)"
    ) from None


f = sp.lambdify(
    (x_symbol, y_symbol),
    f_expression,
    "math"
)


# ------------------------------------------------------------
# Exact solution
# ------------------------------------------------------------

if exact_input.strip() != "":

    exact_expression = convert_expression(
        exact_input
    )

    exact = sp.lambdify(
        x_symbol,
        exact_expression,
        "math"
    )

    HAS_EXACT_FUNCTION = True

else:

    exact_expression = None

    exact = None

    HAS_EXACT_FUNCTION = False


# ------------------------------------------------------------
# Exact value at xf
# ------------------------------------------------------------

if HAS_EXACT_FUNCTION:

    exact_at_xf = exact(xf)

elif exact_at_xf_given is not None:

    exact_at_xf = exact_at_xf_given

else:

    exact_at_xf = None


# ============================================================
# FORMATTING FUNCTIONS
# ============================================================


def formula_string(expression):

    """
    Converts a SymPy expression into a readable
    one-line mathematical expression.
    """

    return sp.sstr(
        sp.simplify(expression)
    )


def formula_pretty(expression):

    """
    Prints a SymPy expression in a more mathematical
    terminal format.
    """

    return sp.pretty(
        sp.simplify(expression)
    )


# ============================================================
# SHOW PROBLEM
# ============================================================


def show_problem():

    print()
    print("=" * 70)
    print("                       PROBLEM DEFINITION")
    print("=" * 70)

    print()

    print("Initial condition:")

    print(
        f"    y({x0}) = {y0}"
    )

    print()

    print("Interval:")

    print(
        f"    {x0} <= x <= {xf}"
    )

    print()

    print("Initial step size:")

    print(
        f"    h0 = {h0}"
    )

    print()

    print("Differential equation:")

    print(
        f"    y' = {formula_string(f_expression)}"
    )

    print()

    print("Differential equation in mathematical form:")

    print(
        "    y' ="
    )

    print(
        formula_pretty(f_expression)
    )

    # --------------------------------------------------------
    # Exact solution exists
    # --------------------------------------------------------

    if HAS_EXACT_FUNCTION:

        print()

        print("Exact solution:")

        print(
            f"    y(x) = "
            f"{formula_string(exact_expression)}"
        )

        print()

        print(
            "Exact solution in mathematical form:"
        )

        print(
            "    y(x) ="
        )

        print(
            formula_pretty(exact_expression)
        )

        print()

        print("Exact value at xf:")

        print(
            f"    y({xf}) = "
            f"{exact_at_xf:.12f}"
        )

    # --------------------------------------------------------
    # Only exact value is available
    # --------------------------------------------------------

    elif exact_at_xf is not None:

        print()

        print("Exact solution:")

        print(
            "    Not available"
        )

        print()

        print("Given exact value at xf:")

        print(
            f"    y({xf}) = "
            f"{exact_at_xf}"
        )

    # --------------------------------------------------------
    # No exact information
    # --------------------------------------------------------

    else:

        print()

        print("Exact solution:")

        print(
            "    Not available"
        )

    print()

    print("=" * 70)


# ============================================================
# NUMERICAL INPUT
# ============================================================


def get_number(prompt, default=None):

    """
    Allows the user to enter:

        0.25

    or expressions such as:

        5*pi/24
        pi/4
        2*pi

    Pressing ENTER uses the default value.
    """

    while True:

        user_input = input(prompt).strip()

        if user_input == "":

            if default is not None:

                return default

            print(
                "A value is required."
            )

            continue

        try:

            value = float(
                eval(
                    user_input,
                    {
                        "__builtins__": {}
                    },
                    {
                        "math": math,
                        "pi": math.pi
                    }
                )
            )

            return value

        except Exception:

            print(
                "Invalid input."
            )

            print(
                "You may enter a number or "
                "an expression such as 5*pi/24."
            )


# ============================================================
# EULER METHOD
# ============================================================


def euler(f, x0, y0, xf, h):

    if not math.isfinite(h) or h <= 0 or xf <= x0:
        raise ValueError("Require finite h > 0 and xf > x0.")

    n = round((xf - x0) / h)
    if n < 1 or not math.isclose(n * h, xf - x0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("h must divide the interval: use h = (xf - x0) / n.")

    h = (
        (xf - x0) / n
    )

    x = x0
    y = y0

    x_values = [x]
    y_values = [y]

    for i in range(n):

        k1 = f(x, y)

        y = y + h * k1

        x = x + h

        x_values.append(x)
        y_values.append(y)

    return (
        x_values,
        y_values,
        h
    )


# ============================================================
# PREDICTOR-CORRECTOR METHOD
# ============================================================
#
# Predictor:
#
# y_p = y_i + h f(x_i,y_i)
#
# Corrector:
#
# y_(i+1) =
# y_i + h/2 [
#       f(x_i,y_i)
#       +
#       f(x_(i+1),y_p)
# ]
#
# ============================================================


def predictor_corrector(
    f,
    x0,
    y0,
    xf,
    h,
    corrector_iterations=1
):

    if corrector_iterations < 1 or int(corrector_iterations) != corrector_iterations:
        raise ValueError("The number of corrector iterations must be a positive integer.")
    if not math.isfinite(h) or h <= 0 or xf <= x0:
        raise ValueError("Require finite h > 0 and xf > x0.")

    n = round((xf - x0) / h)
    if n < 1 or not math.isclose(n * h, xf - x0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("h must divide the interval: use h = (xf - x0) / n.")

    h = (
        (xf - x0) / n
    )

    x = x0
    y = y0

    x_values = [x]
    y_values = [y]

    for i in range(n):

        # Predictor
        y_predictor = (
            y
            +
            h * f(x, y)
        )

        # Keep y_i fixed while iterating the corrector (A-3, Example 10).
        y_corrected = y_predictor
        for k in range(corrector_iterations):
            y_corrected = y + (h / 2) * (f(x, y) + f(x + h, y_corrected))
        y = y_corrected

        x = x + h

        x_values.append(x)
        y_values.append(y)

    return (
        x_values,
        y_values,
        h
    )


# ============================================================
# RK2.1 — HEUN METHOD
# ============================================================
#
# This is the RK2.1 form used for the examples discussed.
#
# k1 = f(x_i,y_i)
#
# k2 = f(x_i+h, y_i+h*k1)
#
# y_(i+1) =
# y_i + h/2 (k1+k2)
#
# ============================================================


def rk2_1(
    f,
    x0,
    y0,
    xf,
    h
):

    if not math.isfinite(h) or h <= 0 or xf <= x0:
        raise ValueError("Require finite h > 0 and xf > x0.")

    n = round((xf - x0) / h)
    if n < 1 or not math.isclose(n * h, xf - x0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("h must divide the interval: use h = (xf - x0) / n.")

    h = (
        (xf - x0) / n
    )

    x = x0
    y = y0

    x_values = [x]
    y_values = [y]

    for i in range(n):

        # First slope
        k1 = f(
            x,
            y
        )

        # Second slope
        k2 = f(
            x + h,
            y + h * k1
        )

        # RK2.1 / Heun update
        y = (
            y
            +
            (h / 2)
            *
            (
                k1 + k2
            )
        )

        x = x + h

        x_values.append(x)
        y_values.append(y)

    return (
        x_values,
        y_values,
        h
    )


# ============================================================
# RK2.2 — RALSTON METHOD (LESSON RK2.2)
# ============================================================
#
# k1 = f(x_i,y_i)
#
# k2 =
# f(
#   x_i+3*h/4,
#   y_i+3*h*k1/4
# )
#
# y_(i+1) = y_i + h*(k1/3 + 2*k2/3)
#
# ============================================================


def rk2_2(
    f,
    x0,
    y0,
    xf,
    h
):

    if not math.isfinite(h) or h <= 0 or xf <= x0:
        raise ValueError("Require finite h > 0 and xf > x0.")

    n = round((xf - x0) / h)
    if n < 1 or not math.isclose(n * h, xf - x0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("h must divide the interval: use h = (xf - x0) / n.")

    h = (
        (xf - x0) / n
    )

    x = x0
    y = y0

    x_values = [x]
    y_values = [y]

    for i in range(n):

        k1 = f(
            x,
            y
        )

        k2 = f(
            x + 3 * h / 4,
            y + 3 * h * k1 / 4
        )

        y = (
            y
            +
            h * (k1 / 3 + 2 * k2 / 3)
        )

        x = x + h

        x_values.append(x)
        y_values.append(y)

    return (
        x_values,
        y_values,
        h
    )


# ============================================================
# RK3 — CLASSICAL THIRD-ORDER RK
# ============================================================


def rk3(
    f,
    x0,
    y0,
    xf,
    h
):

    if not math.isfinite(h) or h <= 0 or xf <= x0:
        raise ValueError("Require finite h > 0 and xf > x0.")

    n = round((xf - x0) / h)
    if n < 1 or not math.isclose(n * h, xf - x0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("h must divide the interval: use h = (xf - x0) / n.")

    h = (
        (xf - x0) / n
    )

    x = x0
    y = y0

    x_values = [x]
    y_values = [y]

    for i in range(n):

        k1 = f(
            x,
            y
        )

        k2 = f(
            x + h / 2,
            y + h * k1 / 2
        )

        k3 = f(
            x + h,
            y - h * k1 + 2 * h * k2
        )

        y = (
            y
            +
            (h / 6)
            *
            (
                k1
                +
                4 * k2
                +
                k3
            )
        )

        x = x + h

        x_values.append(x)
        y_values.append(y)

    return (
        x_values,
        y_values,
        h
    )


# ============================================================
# RK4 — CLASSICAL FOURTH-ORDER RK
# ============================================================


def rk4(
    f,
    x0,
    y0,
    xf,
    h
):

    if not math.isfinite(h) or h <= 0 or xf <= x0:
        raise ValueError("Require finite h > 0 and xf > x0.")

    n = round((xf - x0) / h)
    if n < 1 or not math.isclose(n * h, xf - x0, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("h must divide the interval: use h = (xf - x0) / n.")

    h = (
        (xf - x0) / n
    )

    x = x0
    y = y0

    x_values = [x]
    y_values = [y]

    for i in range(n):

        k1 = f(
            x,
            y
        )

        k2 = f(
            x + h / 2,
            y + h * k1 / 2
        )

        k3 = f(
            x + h / 2,
            y + h * k2 / 2
        )

        k4 = f(
            x + h,
            y + h * k3
        )

        y = (
            y
            +
            (h / 6)
            *
            (
                k1
                +
                2 * k2
                +
                2 * k3
                +
                k4
            )
        )

        x = x + h

        x_values.append(x)
        y_values.append(y)

    return (
        x_values,
        y_values,
        h
    )


# ============================================================
# METHOD DICTIONARY
# ============================================================
#
# IMPORTANT:
# Every method has the same argument order:
#
#     method(f, x0, y0, xf, h)
#
# ============================================================


methods = {

    "Euler": euler,

    "PC": predictor_corrector,

    "RK2.1": rk2_1,

    "RK2.2": rk2_2,

    "RK3": rk3,

    "RK4": rk4

}


# ============================================================
# ERROR FUNCTIONS
# ============================================================


def final_point_error(y_n):

    if exact_at_xf is None:

        return None
    if exact_at_xf == 0:
        return float("inf")  # Relative percentage error is undefined.

    return (
        abs(
            (exact_at_xf - y_n)
            /
            exact_at_xf
        )
        *
        100
    )


# ------------------------------------------------------------
# Error at every point
# ------------------------------------------------------------


def errors_at_all_points(
    x_values,
    y_values
):

    if not HAS_EXACT_FUNCTION:

        return None

    errors = []

    for xi, yi in zip(
        x_values,
        y_values
    ):

        exact_value = exact(xi)

        if exact_value == 0:

            error = 0.0 if yi == 0 else float("inf")

        else:

            error = (
                abs(
                    (exact_value - yi)
                    /
                    exact_value
                )
                *
                100
            )

        errors.append(error)

    return errors


# ------------------------------------------------------------
# Maximum error
# ------------------------------------------------------------


def maximum_error(
    x_values,
    y_values
):

    errors = errors_at_all_points(
        x_values,
        y_values
    )

    if errors is None:

        return None

    return max(errors)


# ============================================================
# DISPLAY ONE NUMERICAL SOLUTION
# ============================================================


def display_solution(
    method_name,
    h,
    corrector_iterations=1
):

    method = methods[
        method_name
    ]

    if method_name == "PC":
        x_values, y_values, actual_h = predictor_corrector(
            f, x0, y0, xf, h, corrector_iterations
        )
        print(f"Corrector iterations per step: {corrector_iterations}")
    else:
        x_values, y_values, actual_h = method(f, x0, y0, xf, h)

    table = pd.DataFrame({

        "i": range(
            len(x_values)
        ),

        "x": x_values,

        "y": y_values

    })

    print()

    print("=" * 70)

    print(
        f"                         {method_name}"
    )

    print("=" * 70)

    print()

    print(
        f"Step size h = {actual_h}"
    )

    print(
        f"Number of steps n = "
        f"{len(x_values) - 1}"
    )

    print()

    # --------------------------------------------------------
    # Table
    # --------------------------------------------------------

    if len(table) <= 40:

        print(
            table.to_string(
                index=False,
                float_format=lambda z:
                    f"{z:.10f}"
            )
        )

    else:

        print(
            "First 10 points:"
        )

        print(
            table.head(10).to_string(
                index=False,
                float_format=lambda z:
                    f"{z:.10f}"
            )
        )

        print()

        print("...")

        print()

        print(
            "Last 10 points:"
        )

        print(
            table.tail(10).to_string(
                index=False,
                float_format=lambda z:
                    f"{z:.10f}"
            )
        )

    # --------------------------------------------------------
    # Final value
    # --------------------------------------------------------

    print()

    print(
        f"Numerical y({xf}) = "
        f"{y_values[-1]:.12f}"
    )

    # --------------------------------------------------------
    # Final-point error
    # --------------------------------------------------------

    if exact_at_xf is not None:

        error = final_point_error(
            y_values[-1]
        )

        print(
            f"Exact y({xf}) = "
            f"{exact_at_xf:.12f}"
        )

        print(
            f"Final-point absolute % error = "
            f"{error:.10f}%"
        )

    # --------------------------------------------------------
    # Maximum error
    # --------------------------------------------------------

    if HAS_EXACT_FUNCTION:

        max_error = maximum_error(
            x_values,
            y_values
        )

        print(
            f"Maximum absolute % error = "
            f"{max_error:.10f}%"
        )

    print()

    print("=" * 70)


# ============================================================
# COMPARE ALL METHODS
# ============================================================


def compare_all_methods(h):

    print()

    print("=" * 90)

    print(
        "                         METHOD COMPARISON"
    )

    print("=" * 90)

    rows = []

    for method_name, method in methods.items():

        x_values, y_values, actual_h = method(
            f,
            x0,
            y0,
            xf,
            h
        )

        y_final = y_values[-1]

        error = final_point_error(
            y_final
        )

        if HAS_EXACT_FUNCTION:

            max_error = maximum_error(
                x_values,
                y_values
            )

        else:

            max_error = None

        rows.append({

            "Method": method_name,

            "n": len(x_values) - 1,

            "h": actual_h,

            "y(xf)": y_final,

            "Final % Error": error,

            "Maximum % Error": max_error

        })

    table = pd.DataFrame(rows)

    print()

    print(
        table.to_string(
            index=False,
            float_format=lambda z:
                f"{z:.10f}"
        )
    )

    print()

    print("=" * 90)


# ============================================================
# FINAL-POINT ERROR TABLE
# ============================================================


def final_point_error_table(
    method_name,
    n_values
):

    if exact_at_xf is None:

        print()

        print(
            "An exact value at xf is required "
            "to calculate percentage error."
        )

        return

    method = methods[
        method_name
    ]

    rows = []

    for n in n_values:

        h = (
            (xf - x0) / n
        )

        x_values, y_values, actual_h = method(
            f,
            x0,
            y0,
            xf,
            h
        )

        y_final = y_values[-1]

        error = final_point_error(
            y_final
        )

        rows.append({

            "n": n,

            "h": actual_h,

            "y_n": y_final,

            "Absolute % Error": error

        })

    table = pd.DataFrame(rows)

    print()

    print("=" * 80)

    print(
        f"              FINAL-POINT ERROR: {method_name}"
    )

    print("=" * 80)

    print()

    print(
        table.to_string(
            index=False,
            float_format=lambda z:
                f"{z:.10f}"
        )
    )

    print()

    print("=" * 80)


# ============================================================
# FIND MINIMUM n
# ============================================================
#
# Requirement:
#
# Find minimum n such that:
#
#       absolute percentage error <= tolerance
#
# ============================================================


def find_minimum_n(
    method_name,
    tolerance=0.001,
    maximum_n=1000000
):

    if exact_at_xf is None:

        print()

        print(
            "Cannot calculate minimum n because "
            "no exact value at xf is available."
        )

        return None

    method = methods[
        method_name
    ]

    print()

    print(
        f"Searching for minimum n "
        f"for {method_name}..."
    )

    # --------------------------------------------------------
    # Find upper bound
    # --------------------------------------------------------

    print("This search assumes final-point error decreases as n increases.")
    n = 1

    while True:

        h = (
            (xf - x0) / n
        )

        x_values, y_values, actual_h = method(
            f,
            x0,
            y0,
            xf,
            h
        )

        error = final_point_error(
            y_values[-1]
        )

        if error <= tolerance:

            break

        n *= 2

        if n > maximum_n:

            print()

            print(
                "Maximum n was exceeded."
            )

            return None

    upper = n

    lower = max(
        1,
        n // 2
    )

    # --------------------------------------------------------
    # Binary search
    # --------------------------------------------------------

    while lower < upper:

        middle = (
            lower
            +
            upper
        ) // 2

        h = (
            (xf - x0) / middle
        )

        x_values, y_values, actual_h = method(
            f,
            x0,
            y0,
            xf,
            h
        )

        error = final_point_error(
            y_values[-1]
        )

        if error <= tolerance:

            upper = middle

        else:

            lower = middle + 1

    # --------------------------------------------------------
    # Final result
    # --------------------------------------------------------

    n = lower

    h = (
        (xf - x0) / n
    )

    x_values, y_values, actual_h = method(
        f,
        x0,
        y0,
        xf,
        h
    )

    y_final = y_values[-1]

    error = final_point_error(
        y_final
    )

    print()

    print("=" * 80)

    print(
        "                         MINIMUM n"
    )

    print("=" * 80)

    print(
        f"Method:              {method_name}"
    )

    print(
        f"Minimum n:           {n}"
    )

    print(
        f"h:                   {actual_h:.12f}"
    )

    print(
        f"Numerical y(xf):     "
        f"{y_final:.12f}"
    )

    print(
        f"Exact y(xf):         "
        f"{exact_at_xf:.12f}"
    )

    print(
        f"Absolute % error:    "
        f"{error:.12f}%"
    )

    print(
        f"Required tolerance:  "
        f"{tolerance}%"
    )

    # --------------------------------------------------------
    # Verify n-1
    # --------------------------------------------------------

    if n > 1:

        previous_n = n - 1

        previous_h = (
            (xf - x0)
            /
            previous_n
        )

        x_previous, y_previous, _ = method(
            f,
            x0,
            y0,
            xf,
            previous_h
        )

        previous_error = final_point_error(
            y_previous[-1]
        )

        print()

        print(
            "Minimum-n verification:"
        )

        print(
            f"n = {previous_n}: "
            f"error = "
            f"{previous_error:.12f}%"
        )

        print(
            f"n = {n}: "
            f"error = "
            f"{error:.12f}%"
        )

    print()

    print(
        f"Therefore, n = {n} "
        f"is the minimum under the decreasing-error assumption."
    )

    print()

    print("=" * 80)

    return (
        n,
        actual_h,
        y_final,
        error
    )

# --------------------------------------------------------
# QUADRATIC LANGRAGIAN H-MAX
# --------------------------------------------------------


def lagrange_hmax(method, f, exact, x0, y0, xf, initial_h=None, tolerance=1.0):
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
    print("\nLAGRANGE INTERPOLATION FOR h_max")
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


# ============================================================
# EULER m-REFINEMENT
# ============================================================
#
# h_m = h0 / 2^m
#
# Continue until:
#
# maximum absolute percentage error <= tolerance
#
# ============================================================


def euler_m_refinement(
    tolerance=1.0,
    maximum_m=20
):

    if not HAS_EXACT_FUNCTION:

        print()

        print(
            "Euler m-refinement requires "
            "an exact solution formula."
        )

        return None

    rows = []

    answer = None

    for m in range(
        maximum_m + 1
    ):

        refined_h = (
            h0
            /
            (2 ** m)
        )

        x_values, y_values, actual_h = euler(
            f,
            x0,
            y0,
            xf,
            refined_h
        )

        max_error = maximum_error(
            x_values,
            y_values
        )

        n = (
            len(x_values) - 1
        )

        rows.append({

            "m": m,

            "n": n,

            "h_m": actual_h,

            "max % error": max_error

        })

        if max_error <= tolerance:

            answer = (
                m,
                n,
                actual_h,
                max_error
            )

            break

    table = pd.DataFrame(rows)

    print()

    print("=" * 80)

    print(
        "                       EULER m-REFINEMENT"
    )

    print("=" * 80)

    print()

    print(
        table.to_string(
            index=False,
            float_format=lambda z:
                f"{z:.10f}"
        )
    )

    if answer is not None:

        m, n, h, error = answer

        print()

        print(
            f"Required m = {m}"
        )

        print(
            f"h_m = {h:.12f}"
        )

        print(
            f"Number of steps n = {n}"
        )

        print(
            f"Maximum absolute % error = "
            f"{error:.10f}%"
        )

        print()

        print(
            f"Therefore:"
        )

        print(
            f"h_max = h0 / 2^{m}"
        )

    else:

        print()

        print(
            "Tolerance was not reached."
        )

    print()

    print("=" * 80)

    return (
        table,
        answer
    )


# ============================================================
# RK4 m-REFINEMENT
# ============================================================
#
# h_m = h0 / 2^m
#
# Compare two successive RK4 solutions.
#
# Because the fine mesh has twice as many points:
#
#       coarse point i
#
# corresponds to:
#
#       fine point 2i
#
# Stop when:
#
#       maximum % change < tolerance
#
# ============================================================


def rk4_m_refinement(
    tolerance=0.5,
    maximum_m=20
):

    rows = []

    previous_y = None

    answer = None

    for m in range(
        maximum_m + 1
    ):

        refined_h = (
            h0
            /
            (2 ** m)
        )

        x_values, y_values, actual_h = rk4(
            f,
            x0,
            y0,
            xf,
            refined_h
        )

        n = (
            len(x_values) - 1
        )

        # ----------------------------------------------------
        # First run
        # ----------------------------------------------------

        if previous_y is None:

            max_change = None

        # ----------------------------------------------------
        # Compare with previous run
        # ----------------------------------------------------

        else:

            changes = []

            for i in range(
                len(previous_y)
            ):

                coarse_value = (
                    previous_y[i]
                )

                fine_value = (
                    y_values[2 * i]
                )

                if fine_value == 0:

                    if coarse_value == 0:

                        change = 0.0

                    else:

                        change = float("inf")

                else:

                    change = (
                        abs(
                            (
                                fine_value
                                -
                                coarse_value
                            )
                            /
                            fine_value
                        )
                        *
                        100
                    )

                changes.append(change)

            max_change = max(
                changes
            )

        rows.append({

            "m": m,

            "n": n,

            "h_m": actual_h,

            "max % change": max_change

        })

        # ----------------------------------------------------
        # Check stopping condition
        # ----------------------------------------------------

        if (
            max_change is not None
            and
            max_change < tolerance
        ):

            answer = (
                m,
                n,
                actual_h,
                max_change
            )

            break

        previous_y = y_values

    table = pd.DataFrame(rows)

    print()

    print("=" * 80)

    print(
        "                        RK4 m-REFINEMENT"
    )

    print("=" * 80)

    print()

    print(
        table.to_string(
            index=False,
            float_format=lambda z:
                (
                    "N/A"
                    if pd.isna(z)
                    else f"{z:.10f}"
                )
        )
    )

    if answer is not None:

        m, n, h, change = answer

        print()

        print(
            f"Required m = {m}"
        )

        print(
            f"h_m = {h:.12f}"
        )

        print(
            f"Number of steps n = {n}"
        )

        print(
            f"Maximum % change = "
            f"{change:.10f}%"
        )

        print()

        print(
            "RK4 m-refinement requirement "
            "has been satisfied."
        )

    else:

        print()

        print(
            "The required tolerance was not "
            "reached."
        )

    print()

    print("=" * 80)

    return (
        table,
        answer
    )

# ============================================================
# RK4 successive m errors
# ============================================================


def rk4_successive_m_error(f, x0, y0, xf, h0):
    """
    RK4 successive-m refinement.

    Compares two successive RK4 solutions:
        m = k
        m = k + 1

    The finer solution is evaluated at the x-values
    belonging to the coarser solution.

    Stops when maximum approximate relative error
    is less than or equal to the user-specified tolerance.
    """

    # ---------------------------------------------------------
    # USER INPUTS
    # ---------------------------------------------------------

    start_m = int(input("Enter starting m: "))
    tolerance_percent = float(
        input("Enter maximum allowable error (%): ")
    )

    if start_m < 0:
        print("Starting m must be 0 or greater.")
        return

    if tolerance_percent <= 0:
        print("Error tolerance must be greater than 0.")
        return

    print("\nRK4 SUCCESSIVE m - APPROXIMATE ERROR CRITERION")
    print("-" * 75)

    m = start_m

    while True:

        # -----------------------------------------------------
        # Number of intervals
        #
        # Each increase in m doubles the number of intervals.
        #
        # m = 0 → n
        # m = 1 → 2n
        # m = 2 → 4n
        # -----------------------------------------------------

        n0 = round((xf - x0) / h0)

        n_coarse = n0 * (2 ** m)
        n_fine = n0 * (2 ** (m + 1))

        h_coarse = (xf - x0) / n_coarse
        h_fine = (xf - x0) / n_fine

        # -----------------------------------------------------
        # RK4 with m
        # -----------------------------------------------------

        result_coarse = rk4(
            f,
            x0,
            y0,
            xf,
            h_coarse
        )

        x_coarse = result_coarse[0]
        y_coarse = result_coarse[1]

        # -----------------------------------------------------
        # RK4 with m + 1
        # -----------------------------------------------------

        result_fine = rk4(
            f,
            x0,
            y0,
            xf,
            h_fine
        )

        x_fine = result_fine[0]
        y_fine = result_fine[1]

        # -----------------------------------------------------
        # Compare the finer solution at the coarse x-points
        #
        # Because the fine grid has twice as many intervals,
        # every second fine point corresponds to a coarse point.
        #
        # coarse point 0 → fine point 0
        # coarse point 1 → fine point 2
        # coarse point 2 → fine point 4
        # etc.
        # -----------------------------------------------------

        errors = []

        for i in range(len(y_coarse)):

            fine_index = 2 * i

            y_m = y_coarse[i]
            y_m_plus_1 = y_fine[fine_index]

            # -------------------------------------------------
            # Calculate approximate relative error
            #
            # εa = |(y^(m+1) - y^m) / y^(m+1)|
            # -------------------------------------------------

            if y_m_plus_1 == 0:

                if y_m == 0:
                    epsilon_a = 0.0
                else:
                    epsilon_a = float("inf")

            else:

                epsilon_a = abs(
                    (y_m_plus_1 - y_m) / y_m_plus_1
                )

            errors.append(epsilon_a)

        # -----------------------------------------------------
        # Maximum approximate relative error
        # -----------------------------------------------------

        maximum_error = max(errors)

        maximum_error_percent = maximum_error * 100

        # -----------------------------------------------------
        # Display results
        # -----------------------------------------------------

        print(
            f"m = {m}, m+1 = {m + 1}"
        )

        print(
            f"n = {n_coarse}, "
            f"n(next) = {n_fine}"
        )

        print(
            f"h = {h_coarse:.10f}, "
            f"h(next) = {h_fine:.10f}"
        )

        print(
            f"maximum εa = "
            f"{maximum_error_percent:.10f}%"
        )

        # -----------------------------------------------------
        # Check criterion
        # -----------------------------------------------------

        if maximum_error_percent <= tolerance_percent:

            print("\nCriterion satisfied!")
            print(
                f"Final m = {m + 1}"
            )
            print(
                f"Maximum εa = "
                f"{maximum_error_percent:.10f}%"
            )
            print(
                f"Required εa = "
                f"{tolerance_percent:.10f}%"
            )

            return (
                m + 1,
                maximum_error_percent
            )

        # -----------------------------------------------------
        # Try next refinement
        # -----------------------------------------------------

        m += 1


# ============================================================
# SHOW ERRORS AT ALL POINTS
# ============================================================


def show_all_point_errors():

    if not HAS_EXACT_FUNCTION:

        print()

        print(
            "This calculation requires "
            "an exact solution formula."
        )

        return

    print()

    print(
        "Available methods:"
    )

    for name in methods:

        print(
            f"    {name}"
        )

    print()

    method_name = input(
        "Enter method name exactly: "
    ).strip()

    if method_name not in methods:

        print(
            "Invalid method."
        )

        return

    h = get_number(
        f"Enter h "
        f"(press Enter for {h0}): ",
        h0
    )

    method = methods[
        method_name
    ]

    x_values, y_values, actual_h = method(
        f,
        x0,
        y0,
        xf,
        h
    )

    errors = errors_at_all_points(
        x_values,
        y_values
    )

    table = pd.DataFrame({

        "i":
            range(
                len(x_values)
            ),

        "x":
            x_values,

        "Numerical y":
            y_values,

        "Exact y":
            [
                exact(x)
                for x in x_values
            ],

        "Absolute % Error":
            errors

    })

    print()

    print("=" * 90)

    print(
        f"                    {method_name} ERROR TABLE"
    )

    print("=" * 90)

    print()

    if len(table) <= 40:

        print(
            table.to_string(
                index=False,
                float_format=lambda z:
                    f"{z:.10f}"
            )
        )

    else:

        print(
            "First 10 points:"
        )

        print(
            table.head(10).to_string(
                index=False,
                float_format=lambda z:
                    f"{z:.10f}"
            )
        )

        print()

        print("...")

        print()

        print(
            "Last 10 points:"
        )

        print(
            table.tail(10).to_string(
                index=False,
                float_format=lambda z:
                    f"{z:.10f}"
            )
        )

    print()

    print(
        f"Maximum absolute % error = "
        f"{max(errors):.10f}%"
    )

    print()

    print("=" * 90)


# ============================================================
# PLOT NUMERICAL SOLUTIONS
# ============================================================


def plot_solutions():

    default_n = round(
        (xf - x0) / h0
    )

    print()

    n = int(
        get_number(
            "Enter n for the plot "
            f"(Enter for {default_n}): ",
            default_n
        )
    )

    h = (
        (xf - x0) / n
    )

    plt.figure(
        figsize=(10, 6)
    )

    # --------------------------------------------------------
    # Numerical methods
    # --------------------------------------------------------

    for method_name, method in methods.items():

        x_values, y_values, actual_h = method(
            f,
            x0,
            y0,
            xf,
            h
        )

        plt.plot(
            x_values,
            y_values,
            label=method_name
        )

    # --------------------------------------------------------
    # Exact solution
    # --------------------------------------------------------

    if HAS_EXACT_FUNCTION:

        x_exact = [

            x0
            +
            i
            *
            (xf - x0)
            /
            500

            for i in range(501)

        ]

        y_exact_values = [

            exact(x)

            for x in x_exact

        ]

        plt.plot(
            x_exact,
            y_exact_values,
            label="Exact"
        )

    plt.xlabel("x")

    plt.ylabel("y")

    plt.title(
        "Numerical Solutions"
    )

    plt.legend()

    plt.grid(True)

    plt.show()


# ============================================================
# PLOT EULER MAXIMUM ERROR VS m
# ============================================================


def plot_euler_error_vs_m():

    if not HAS_EXACT_FUNCTION:

        print()

        print(
            "This plot requires "
            "an exact solution formula."
        )

        return

    table, answer = euler_m_refinement(
        tolerance=1.0,
        maximum_m=15
    )

    plt.figure(
        figsize=(10, 6)
    )

    plt.plot(
        table["m"],
        table["max % error"],
        marker="o"
    )

    plt.axhline(
        1.0,
        linestyle="--",
        label="1% tolerance"
    )

    plt.xlabel("m")

    plt.ylabel(
        "Maximum absolute % error"
    )

    plt.title(
        "Euler m-Refinement"
    )

    plt.legend()

    plt.grid(True)

    plt.show()


# ============================================================
# PLOT FINAL-POINT ERROR VS n
# ============================================================


def plot_final_error_vs_n():

    if exact_at_xf is None:

        print()

        print(
            "An exact value at xf is required."
        )

        return

    print()

    print(
        "Available methods:"
    )

    for name in methods:

        print(
            f"    {name}"
        )

    print()

    method_name = input(
        "Enter method name exactly: "
    ).strip()

    if method_name not in methods:

        print(
            "Invalid method."
        )

        return

    max_n = int(
        get_number(
            "Enter maximum n "
            "(default 100): ",
            100
        )
    )

    n_values = []

    n = 1

    while n <= max_n:

        n_values.append(n)

        n *= 2

    errors = []

    for n in n_values:

        h = (
            (xf - x0) / n
        )

        x_values, y_values, actual_h = methods[
            method_name
        ](
            f,
            x0,
            y0,
            xf,
            h
        )

        error = final_point_error(
            y_values[-1]
        )

        errors.append(error)

    plt.figure(
        figsize=(10, 6)
    )

    plt.plot(
        n_values,
        errors,
        marker="o"
    )

    plt.xlabel("n")

    plt.ylabel(
        "Final-point absolute % error"
    )

    plt.title(
        f"{method_name}: Final-Point Error vs n"
    )

    plt.xscale(
        "log",
        base=2
    )

    plt.yscale(
        "log"
    )

    plt.grid(True)

    plt.show()


# ============================================================
# PLOT COMPARISON OF ERRORS
# ============================================================


def plot_error_comparison(f, exact, x0, y0, xf, h, methods):
    """
    Plot absolute percentage error versus x
    for any selected numerical methods.
    """

    if exact is None:
        print("Exact solution is required for error comparison.")
        return

    # ---------------------------------------------------------
    # Ask which methods to compare
    # ---------------------------------------------------------

    print("\nAVAILABLE METHODS:")
    for name in methods:
        print(f"- {name}")

    method_input = input(
        "\nEnter methods to compare "
        "(separated by commas): "
    )

    selected_names = [
        name.strip()
        for name in method_input.split(",")
    ]

    # ---------------------------------------------------------
    # Validate method names
    # ---------------------------------------------------------

    selected_methods = []

    for name in selected_names:

        if name in methods:
            selected_methods.append(
                (name, methods[name])
            )

        else:
            print(
                f"Warning: '{name}' is not a valid method."
            )

    if not selected_methods:
        print("No valid methods were selected.")
        return

    # ---------------------------------------------------------
    # Plot
    # ---------------------------------------------------------

    import matplotlib.pyplot as plt

    plt.figure(figsize=(9, 6))

    # ---------------------------------------------------------
    # Calculate and plot error for each method
    # ---------------------------------------------------------

    for method_name, method in selected_methods:

        result = method(
            f,
            x0,
            y0,
            xf,
            h
        )

        x_values = result[0]
        y_values = result[1]

        errors = []

        for x, y in zip(x_values, y_values):

            exact_y = exact(x)

            if exact_y != 0:

                error = abs(
                    (exact_y - y) / exact_y
                ) * 100

            else:

                error = 0.0 if y == 0 else float("inf")

            errors.append(error)

        plt.plot(
            x_values,
            errors,
            marker="o",
            linestyle="--",
            label=method_name
        )

    # ---------------------------------------------------------
    # Graph formatting
    # ---------------------------------------------------------

    plt.xlabel("x")
    plt.ylabel("Absolute Percentage Error (%)")

    plt.title(
        "Comparison of Numerical Method Errors"
    )

    plt.legend()
    plt.grid(True)

    plt.show()

# ============================================================
# PROCEED TO MENU
# ============================================================


def proceed_to_menu():

    while True:

        choice = input(
            "\nProceed to menu? (y/n): "
        ).strip().lower()

        if choice == "y":

            return True

        elif choice == "n":

            print()

            print(
                "Program ended."
            )

            return False

        else:

            print(
                "Please enter y or n."
            )

# ============================================================
# FINITE DIFFERENCE APPROXIMATION
# ============================================================

def finite_difference_approximations():
    """
    Calculate first-derivative approximations using
    forward, backward, and central finite differences.

    The exact derivative is automatically generated
    from the user's input function using SymPy.
    """

    print("\nFINITE DIFFERENCE APPROXIMATIONS")
    print("-" * 60)

    # ---------------------------------------------------------
    # USER INPUTS
    # ---------------------------------------------------------

    function_input = input(
        "Enter f(x): "
    )

    x_value = get_number(
        "Enter x-value: "
    )

    h = get_number(
        "Enter h: "
    )

    if h <= 0:
        print("Error: h must be positive.")
        return

    # ---------------------------------------------------------
    # CONVERT USER'S FUNCTION INTO A SYMPY EXPRESSION
    # ---------------------------------------------------------

    expression = convert_expression(
        function_input
    )

    # ---------------------------------------------------------
    # CREATE NUMERICAL FUNCTION f(x)
    # ---------------------------------------------------------

    function = sp.lambdify(
        x_symbol,
        expression,
        "math"
    )

    # ---------------------------------------------------------
    # AUTOMATICALLY CALCULATE EXACT FIRST DERIVATIVE
    # ---------------------------------------------------------

    exact_derivative_expression = sp.diff(
        expression,
        x_symbol
    )

    # Convert exact derivative into a numerical function
    exact_derivative = sp.lambdify(
        x_symbol,
        exact_derivative_expression,
        "math"
    )

    # ---------------------------------------------------------
    # CALCULATE FUNCTION VALUES
    # ---------------------------------------------------------

    f_x = function(x_value)

    f_x_plus_h = function(
        x_value + h
    )

    f_x_minus_h = function(
        x_value - h
    )

    # ---------------------------------------------------------
    # FORWARD DIFFERENCE
    #
    # f'(x) ≈ [f(x+h) - f(x)] / h
    # ---------------------------------------------------------

    forward = (
        f_x_plus_h - f_x
    ) / h

    # ---------------------------------------------------------
    # BACKWARD DIFFERENCE
    #
    # f'(x) ≈ [f(x) - f(x-h)] / h
    # ---------------------------------------------------------

    backward = (
        f_x - f_x_minus_h
    ) / h

    # ---------------------------------------------------------
    # CENTRAL DIFFERENCE
    #
    # f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
    # ---------------------------------------------------------

    central = (
        f_x_plus_h - f_x_minus_h
    ) / (2 * h)

    # ---------------------------------------------------------
    # EXACT DERIVATIVE VALUE
    # ---------------------------------------------------------

    exact_value = exact_derivative(
        x_value
    )

    # ---------------------------------------------------------
    # ABSOLUTE ERRORS
    # ---------------------------------------------------------

    forward_error = abs(
        exact_value - forward
    )

    backward_error = abs(
        exact_value - backward
    )

    central_error = abs(
        exact_value - central
    )

    # ---------------------------------------------------------
    # PERCENTAGE ERRORS
    # ---------------------------------------------------------

    if exact_value != 0:

        forward_percent_error = (
            forward_error /
            abs(exact_value)
        ) * 100

        backward_percent_error = (
            backward_error /
            abs(exact_value)
        ) * 100

        central_percent_error = (
            central_error /
            abs(exact_value)
        ) * 100

    else:

        forward_percent_error = None
        backward_percent_error = None
        central_percent_error = None

    # ---------------------------------------------------------
    # DISPLAY FUNCTION AND DERIVATIVE
    # ---------------------------------------------------------

    print("\nFUNCTION")
    print("-" * 60)

    print(
        f"f(x) = {sp.sstr(expression)}"
    )

    print("\nEXACT DERIVATIVE")
    print("-" * 60)

    print(
        f"f'(x) = "
        f"{sp.sstr(exact_derivative_expression)}"
    )

    print(
        f"f'({x_value:.10f}) = "
        f"{exact_value:.10f}"
    )

    # ---------------------------------------------------------
    # DISPLAY FUNCTION VALUES
    # ---------------------------------------------------------

    print("\nFUNCTION VALUES")
    print("-" * 60)

    print(
        f"f(x)       = {f_x:.10f}"
    )

    print(
        f"f(x+h)     = {f_x_plus_h:.10f}"
    )

    print(
        f"f(x-h)     = {f_x_minus_h:.10f}"
    )

    # ---------------------------------------------------------
    # DISPLAY APPROXIMATIONS
    # ---------------------------------------------------------

    print("\nDERIVATIVE APPROXIMATIONS")
    print("-" * 60)

    print(
        f"Forward difference  = "
        f"{forward:.10f}"
    )

    print(
        f"Backward difference = "
        f"{backward:.10f}"
    )

    print(
        f"Central difference  = "
        f"{central:.10f}"
    )

    # ---------------------------------------------------------
    # DISPLAY ABSOLUTE ERRORS
    # ---------------------------------------------------------

    print("\nABSOLUTE ERRORS")
    print("-" * 60)

    print(
        f"Forward error  = "
        f"{forward_error:.10f}"
    )

    print(
        f"Backward error = "
        f"{backward_error:.10f}"
    )

    print(
        f"Central error  = "
        f"{central_error:.10f}"
    )

    # ---------------------------------------------------------
    # DISPLAY PERCENTAGE ERRORS
    # ---------------------------------------------------------

    if exact_value != 0:

        print("\nPERCENTAGE ERRORS")
        print("-" * 60)

        print(
            f"Forward % error  = "
            f"{forward_percent_error:.10f}%"
        )

        print(
            f"Backward % error = "
            f"{backward_percent_error:.10f}%"
        )

        print(
            f"Central % error  = "
            f"{central_percent_error:.10f}%"
        )

    else:

        print(
            "\nPercentage errors cannot be calculated "
            "because the exact derivative is zero."
        )


# ============================================================
# FINITE DIFFERENCE APPROXIMATION TABULATED
# ============================================================

def finite_difference_tabulated():
    """
    Finite-difference approximations using tabulated x and y values.

    First derivative:
        Left endpoint  -> Forward difference
        Interior       -> Central difference
        Right endpoint -> Backward difference

    Second derivative:
        Left endpoint  -> Forward difference
        Interior       -> Central difference
        Right endpoint -> Backward difference

    The analytical function is optional.
    If supplied, SymPy automatically generates the exact
    first and second derivatives for error analysis.
    """

    print("\nFINITE DIFFERENCE FROM TABULATED DATA")
    print("-" * 80)

    # =========================================================
    # 1. INPUT X-VALUES
    # =========================================================

    x_input = input(
        "Enter x-values separated by commas: "
    )

    try:
        x_values = [
            float(value.strip())
            for value in x_input.split(",")
        ]
    except ValueError:
        print("Error: Please enter valid numerical x-values.")
        return

    # =========================================================
    # 2. INPUT Y-VALUES
    # =========================================================

    y_input = input(
        "Enter y-values separated by commas: "
    )

    try:
        y_values = [
            float(value.strip())
            for value in y_input.split(",")
        ]
    except ValueError:
        print("Error: Please enter valid numerical y-values.")
        return

    # =========================================================
    # 3. CHECK DATA
    # =========================================================

    if len(x_values) != len(y_values):

        print(
            "\nError: The number of x-values and y-values "
            "must be the same."
        )

        return

    if len(x_values) < 2:

        print(
            "\nError: At least two data points are required."
        )

        return

    # =========================================================
    # 4. CHECK FOR UNIFORM SPACING
    # =========================================================

    h_values = []

    for i in range(len(x_values) - 1):

        h_values.append(
            x_values[i + 1] - x_values[i]
        )

    h = h_values[0]
    if h == 0:
        print("Error: x-values must be distinct.")
        return

    for value in h_values:

        if abs(value - h) > 1e-10:

            print(
                "\nError: The x-values must have "
                "uniform spacing."
            )

            return

    # =========================================================
    # 5. OPTIONAL ANALYTICAL FUNCTION
    # =========================================================

    function_input = input(
        "\nEnter analytical f(x) "
        "(press ENTER if unavailable): "
    )

    exact_first_derivative = None
    exact_second_derivative = None

    exact_first_expression = None
    exact_second_expression = None

    if function_input.strip() != "":

        try:

            # Convert original function into SymPy expression
            expression = convert_expression(
                function_input
            )

            # -------------------------------------------------
            # Automatically generate exact first derivative
            # -------------------------------------------------

            exact_first_expression = sp.diff(
                expression,
                x_symbol
            )

            # -------------------------------------------------
            # Automatically generate exact second derivative
            # -------------------------------------------------

            exact_second_expression = sp.diff(
                expression,
                x_symbol,
                2
            )

            # -------------------------------------------------
            # Convert derivatives into numerical functions
            # -------------------------------------------------

            exact_first_derivative = sp.lambdify(
                x_symbol,
                exact_first_expression,
                "math"
            )

            exact_second_derivative = sp.lambdify(
                x_symbol,
                exact_second_expression,
                "math"
            )

        except Exception as error:

            print(
                "\nError interpreting analytical function:"
            )
            print(error)
            return

    # =========================================================
    # 6. CALCULATE FIRST AND SECOND DERIVATIVES
    # =========================================================

    first_derivatives = []
    second_derivatives = []

    for i in range(len(x_values)):

        # =====================================================
        # FIRST DERIVATIVE
        # =====================================================

        if i == 0:

            # ---------------------------------------------
            # LEFT ENDPOINT
            # Forward difference
            #
            # y'_0 ≈ (y_1 - y_0) / h
            # ---------------------------------------------

            first = (
                y_values[i + 1]
                - y_values[i]
            ) / h

        elif i == len(x_values) - 1:

            # ---------------------------------------------
            # RIGHT ENDPOINT
            # Backward difference
            #
            # y'_n ≈ (y_n - y_(n-1)) / h
            # ---------------------------------------------

            first = (
                y_values[i]
                - y_values[i - 1]
            ) / h

        else:

            # ---------------------------------------------
            # INTERIOR POINT
            # Central difference
            #
            # y'_i ≈ (y_(i+1) - y_(i-1)) / (2h)
            # ---------------------------------------------

            first = (
                y_values[i + 1]
                - y_values[i - 1]
            ) / (2 * h)

        first_derivatives.append(first)

        # =====================================================
        # SECOND DERIVATIVE
        # =====================================================

        if len(y_values) < 3:

            # Not enough data for second derivative
            second_derivatives.append(None)

        elif i == 0:

            # ---------------------------------------------
            # LEFT ENDPOINT
            # Forward second difference
            #
            # y''_0 ≈
            # (y_2 - 2y_1 + y_0) / h²
            # ---------------------------------------------

            second = (
                y_values[i + 2]
                - 2 * y_values[i + 1]
                + y_values[i]
            ) / (h ** 2)

            second_derivatives.append(second)

        elif i == len(x_values) - 1:

            # ---------------------------------------------
            # RIGHT ENDPOINT
            # Backward second difference
            #
            # y''_n ≈
            # (y_n - 2y_(n-1) + y_(n-2)) / h²
            # ---------------------------------------------

            second = (
                y_values[i]
                - 2 * y_values[i - 1]
                + y_values[i - 2]
            ) / (h ** 2)

            second_derivatives.append(second)

        else:

            # ---------------------------------------------
            # INTERIOR POINT
            # Central second difference
            #
            # y''_i ≈
            # (y_(i+1) - 2y_i + y_(i-1)) / h²
            # ---------------------------------------------

            second = (
                y_values[i + 1]
                - 2 * y_values[i]
                + y_values[i - 1]
            ) / (h ** 2)

            second_derivatives.append(second)

    # =========================================================
    # 7. DISPLAY DATA
    # =========================================================

    print("\nDATA")
    print("-" * 80)

    print(
        f"Number of data points = {len(x_values)}"
    )

    print(
        f"Step size h = {h:.10f}"
    )

    # =========================================================
    # 8. DISPLAY EXACT DERIVATIVES
    # =========================================================

    if exact_first_derivative is not None:

        print("\nEXACT DERIVATIVES")
        print("-" * 80)

        print(
            f"f'(x)  = "
            f"{sp.sstr(exact_first_expression)}"
        )

        print(
            f"f''(x) = "
            f"{sp.sstr(exact_second_expression)}"
        )

    # =========================================================
    # 9. DISPLAY NUMERICAL DERIVATIVES
    # =========================================================

    print("\nFINITE DIFFERENCE RESULTS")
    print("-" * 100)

    print(
        f"{'i':>3}"
        f"{'x':>12}"
        f"{'y':>14}"
        f"{'y_prime':>18}"
        f"{'y_double_prime':>20}"
    )

    print("-" * 100)

    for i in range(len(x_values)):

        if second_derivatives[i] is None:

            second_display = "N/A"

        else:

            second_display = (
                f"{second_derivatives[i]:.10f}"
            )

        print(
            f"{i:>3}"
            f"{x_values[i]:>12.6f}"
            f"{y_values[i]:>14.6f}"
            f"{first_derivatives[i]:>18.10f}"
            f"{second_display:>20}"
        )

    # =========================================================
    # 10. ERROR ANALYSIS
    # =========================================================

    if exact_first_derivative is not None:

        print("\nERROR ANALYSIS")
        print("-" * 130)

        print(
            f"{'i':>3}"
            f"{'x':>10}"
            f"{'FD y_prime':>16}"
            f"{'Exact y_prime':>16}"
            f"{'Error %':>14}"
            f"{'FD y_double':>16}"
            f"{'Exact y_double':>16}"
            f"{'Error %':>14}"
        )

        print("-" * 130)

        for i in range(len(x_values)):

            x = x_values[i]

            # -------------------------------------------------
            # Exact derivative values
            # -------------------------------------------------

            exact_first = exact_first_derivative(x)
            exact_second = exact_second_derivative(x)

            # -------------------------------------------------
            # First derivative percentage error
            # -------------------------------------------------

            if exact_first != 0:

                first_error_percent = (
                    abs(
                        first_derivatives[i]
                        - exact_first
                    )
                    / abs(exact_first)
                ) * 100

            else:

                first_error_percent = None

            # -------------------------------------------------
            # Second derivative error
            # -------------------------------------------------

            if second_derivatives[i] is None:

                second_error_percent = None

            elif exact_second != 0:

                second_error_percent = (
                    abs(
                        second_derivatives[i]
                        - exact_second
                    )
                    / abs(exact_second)
                ) * 100

            else:

                second_error_percent = None

            # -------------------------------------------------
            # Display first error
            # -------------------------------------------------

            if first_error_percent is None:

                first_error_display = "N/A"

            else:

                first_error_display = (
                    f"{first_error_percent:.6f}"
                )

            # -------------------------------------------------
            # Display second derivative
            # -------------------------------------------------

            if second_derivatives[i] is None:

                second_fd_display = "N/A"

            else:

                second_fd_display = (
                    f"{second_derivatives[i]:.8f}"
                )

            # -------------------------------------------------
            # Display exact second derivative
            # -------------------------------------------------

            if second_derivatives[i] is None:

                exact_second_display = "N/A"

            else:

                exact_second_display = (
                    f"{exact_second:.8f}"
                )

            # -------------------------------------------------
            # Display second error
            # -------------------------------------------------

            if second_error_percent is None:

                second_error_display = "N/A"

            else:

                second_error_display = (
                    f"{second_error_percent:.6f}"
                )

            print(
                f"{i:>3}"
                f"{x:>10.6f}"
                f"{first_derivatives[i]:>16.8f}"
                f"{exact_first:>16.8f}"
                f"{first_error_display:>14}"
                f"{second_fd_display:>16}"
                f"{exact_second_display:>16}"
                f"{second_error_display:>14}"
            )



# ============================================================
# CENTRAL SECOND DERIVATIVE
# ============================================================


def central_second_derivative_tabulated():
    """
    OPTION 17
    Central Difference Approximation of the Second Derivative

    Formula:
        y''_i ≈ (y_(i-1) - 2y_i + y_(i+1)) / h²

    This formula is applicable only at interior points.
    """

    print("\nCENTRAL DIFFERENCE APPROXIMATION OF SECOND DERIVATIVE")
    print("-" * 75)

    # =========================================================
    # 1. INPUT X-VALUES
    # =========================================================

    x_input = input(
        "Enter x-values separated by commas: "
    )

    try:
        x_values = [
            float(value.strip())
            for value in x_input.split(",")
        ]
    except ValueError:
        print("Error: Please enter valid numerical x-values.")
        return

    # =========================================================
    # 2. INPUT Y-VALUES
    # =========================================================

    y_input = input(
        "Enter y-values separated by commas: "
    )

    try:
        y_values = [
            float(value.strip())
            for value in y_input.split(",")
        ]
    except ValueError:
        print("Error: Please enter valid numerical y-values.")
        return

    # =========================================================
    # 3. CHECK DATA
    # =========================================================

    if len(x_values) != len(y_values):
        print(
            "Error: Number of x-values and y-values "
            "must be the same."
        )
        return

    if len(x_values) < 3:
        print(
            "Error: At least 3 points are required "
            "for central second difference."
        )
        return

    # =========================================================
    # 4. CHECK UNIFORM SPACING
    # =========================================================

    h_values = []

    for i in range(len(x_values) - 1):

        h_values.append(
            x_values[i + 1] - x_values[i]
        )

    h = h_values[0]
    if h == 0:
        print("Error: x-values must be distinct.")
        return

    for value in h_values:

        if abs(value - h) > 1e-10:

            print(
                "Error: x-values must have uniform spacing."
            )
            return

    # =========================================================
    # 5. OPTIONAL ANALYTICAL FUNCTION
    # =========================================================

    function_input = input(
        "\nEnter analytical f(x) "
        "(press ENTER if unavailable): "
    )

    exact_second_derivative = None
    exact_second_expression = None

    if function_input.strip() != "":

        try:

            expression = convert_expression(
                function_input
            )

            # Automatically generate f''(x)
            exact_second_expression = sp.diff(
                expression,
                x_symbol,
                2
            )

            # Convert to numerical function
            exact_second_derivative = sp.lambdify(
                x_symbol,
                exact_second_expression,
                "math"
            )

        except Exception as error:

            print(
                "\nError interpreting analytical function:"
            )
            print(error)
            return

    # =========================================================
    # 6. CENTRAL SECOND-DERIVATIVE CALCULATION
    # =========================================================

    derivative_x_values = []
    second_derivatives = []

    # Central difference is ONLY possible for interior points.
    #
    # i = 1, 2, ..., n-1
    #
    # because we need:
    #
    # y_(i-1), y_i, y_(i+1)

    for i in range(1, len(y_values) - 1):

        second = (
            y_values[i - 1]
            - 2 * y_values[i]
            + y_values[i + 1]
        ) / (h ** 2)

        derivative_x_values.append(
            x_values[i]
        )

        second_derivatives.append(
            second
        )

    # =========================================================
    # 7. DISPLAY DATA
    # =========================================================

    print("\nDATA")
    print("-" * 75)

    print(
        f"Number of data points = {len(x_values)}"
    )

    print(
        f"Step size h = {h:.10f}"
    )

    print(
        "\nCentral second difference is applied "
        "only to interior points."
    )

    # =========================================================
    # 8. DISPLAY EXACT SECOND DERIVATIVE
    # =========================================================

    if exact_second_derivative is not None:

        print("\nEXACT SECOND DERIVATIVE")
        print("-" * 75)

        print(
            f"f''(x) = "
            f"{sp.sstr(exact_second_expression)}"
        )

    # =========================================================
    # 9. DISPLAY APPROXIMATE VALUES
    # =========================================================

    print("\nCENTRAL DIFFERENCE RESULTS")
    print("-" * 75)

    print(
        f"{'i':>5}"
        f"{'x_i':>12}"
        f"{'y_i':>15}"
        f"{'y_i_double_prime':>22}"
    )

    print("-" * 75)

    for j in range(len(second_derivatives)):

        i = j + 1

        print(
            f"{i:>5}"
            f"{x_values[i]:>12.6f}"
            f"{y_values[i]:>15.6f}"
            f"{second_derivatives[j]:>22.8f}"
        )

    # =========================================================
    # 10. ERROR ANALYSIS
    # =========================================================

    if exact_second_derivative is not None:

        error_percentages = []

        print("\nERROR ANALYSIS")
        print("-" * 95)

        print(
            f"{'i':>5}"
            f"{'x_i':>12}"
            f"{'Approx. y_double':>20}"
            f"{'Exact y_double':>20}"
            f"{'Error %':>15}"
        )

        print("-" * 95)

        for j in range(len(second_derivatives)):

            i = j + 1
            x = x_values[i]

            approximate = second_derivatives[j]

            exact = exact_second_derivative(x)

            absolute_error = abs(
                approximate - exact
            )

            if exact != 0:

                error_percent = (
                    absolute_error
                    / abs(exact)
                ) * 100

            else:

                error_percent = None

            error_percentages.append(
                error_percent
            )

            if error_percent is None:

                error_display = "N/A"

            else:

                error_display = (
                    f"{error_percent:.6f}"
                )

            print(
                f"{i:>5}"
                f"{x:>12.6f}"
                f"{approximate:>20.8f}"
                f"{exact:>20.8f}"
                f"{error_display:>15}"
            )

        # =====================================================
        # 11. PLOT ABSOLUTE PERCENTAGE ERRORS
        # =====================================================

        valid_x = []
        valid_errors = []

        for x, error in zip(
            derivative_x_values,
            error_percentages
        ):

            if error is not None:

                valid_x.append(x)
                valid_errors.append(error)

        if len(valid_errors) > 0:

            plt.figure()

            plt.plot(
                valid_x,
                valid_errors,
                marker="o"
            )

            plt.xlabel("x")
            plt.ylabel("Absolute Percentage Error (%)")

            plt.title(
                "Central Difference Approximation of Second Derivative"
            )

            plt.grid(True)

            plt.show()

    else:

        print(
            "\nNo analytical function was supplied."
        )

        print(
            "Percentage errors cannot be calculated."
        )

# ============================================================
# HIGHER ORDER END POINT DIFFERENCE
# ============================================================

def higher_order_endpoint_differences():
    """
    OPTION 18

    Higher-order O(h^2) forward and backward finite-difference
    approximations for the first and second derivatives.

    First derivative:

    Forward O(h^2):
        y'_i = (-3y_i + 4y_(i+1) - y_(i+2)) / (2h)

    Backward O(h^2):
        y'_i = (3y_i - 4y_(i-1) + y_(i-2)) / (2h)

    Second derivative:

    Forward O(h^2):
        y''_i = (2y_i - 5y_(i+1) + 4y_(i+2) - y_(i+3)) / h^2

    Backward O(h^2):
        y''_i = (2y_i - 5y_(i-1) + 4y_(i-2) - y_(i-3)) / h^2

    The analytical function is optional. If supplied,
    SymPy automatically generates the exact first and
    second derivatives.
    """

    print("\nHIGHER-ORDER FORWARD/BACKWARD DIFFERENCES")
    print("-" * 75)

    # =========================================================
    # 1. INPUT X-VALUES
    # =========================================================

    x_input = input(
        "Enter x-values separated by commas: "
    )

    try:

        x_values = [
            float(value.strip())
            for value in x_input.split(",")
        ]

    except ValueError:

        print(
            "Error: Please enter valid numerical x-values."
        )

        return

    # =========================================================
    # 2. INPUT Y-VALUES
    # =========================================================

    y_input = input(
        "Enter y-values separated by commas: "
    )

    try:

        y_values = [
            float(value.strip())
            for value in y_input.split(",")
        ]

    except ValueError:

        print(
            "Error: Please enter valid numerical y-values."
        )

        return

    # =========================================================
    # 3. CHECK DATA
    # =========================================================

    if len(x_values) != len(y_values):

        print(
            "\nError: The number of x-values and y-values "
            "must be the same."
        )

        return

    # Four points are needed for the second derivative
    # higher-order formula.

    if len(x_values) < 4:

        print(
            "\nError: At least 4 data points are required."
        )

        return

    # =========================================================
    # 4. CHECK UNIFORM SPACING
    # =========================================================

    h_values = []

    for i in range(len(x_values) - 1):

        h_values.append(
            x_values[i + 1] - x_values[i]
        )

    h = h_values[0]
    if h == 0:
        print("Error: x-values must be distinct.")
        return

    for value in h_values:

        if abs(value - h) > 1e-10:

            print(
                "\nError: The x-values must have "
                "uniform spacing."
            )

            return

    # =========================================================
    # 5. SELECT FORWARD OR BACKWARD
    # =========================================================

    print("\nChoose endpoint:")
    print("1. Forward difference at x0")
    print("2. Backward difference at xn")

    endpoint_choice = input(
        "Choose an option: "
    )

    # =========================================================
    # 6. OPTIONAL ANALYTICAL FUNCTION
    # =========================================================

    function_input = input(
        "\nEnter analytical f(x) "
        "(press ENTER if unavailable): "
    )

    exact_first_derivative = None
    exact_second_derivative = None

    exact_first_expression = None
    exact_second_expression = None

    if function_input.strip() != "":

        try:

            expression = convert_expression(
                function_input
            )

            # Automatically generate exact first derivative
            exact_first_expression = sp.diff(
                expression,
                x_symbol
            )

            # Automatically generate exact second derivative
            exact_second_expression = sp.diff(
                expression,
                x_symbol,
                2
            )

            exact_first_derivative = sp.lambdify(
                x_symbol,
                exact_first_expression,
                "math"
            )

            exact_second_derivative = sp.lambdify(
                x_symbol,
                exact_second_expression,
                "math"
            )

        except Exception as error:

            print(
                "\nError interpreting analytical function:"
            )

            print(error)

            return

    # =========================================================
    # 7. FORWARD DIFFERENCE AT x0
    # =========================================================

    if endpoint_choice == "1":

        i = 0

        x = x_values[i]

        # -----------------------------------------------------
        # First derivative: O(h^2)
        #
        # y'_0 =
        # (-3y_0 + 4y_1 - y_2) / (2h)
        # -----------------------------------------------------

        higher_first = (
            -3 * y_values[0]
            + 4 * y_values[1]
            - y_values[2]
        ) / (2 * h)

        # -----------------------------------------------------
        # Second derivative: O(h^2)
        #
        # y''_0 =
        # (2y_0 - 5y_1 + 4y_2 - y_3) / h^2
        # -----------------------------------------------------

        higher_second = (
            2 * y_values[0]
            - 5 * y_values[1]
            + 4 * y_values[2]
            - y_values[3]
        ) / (h ** 2)

        # -----------------------------------------------------
        # Basic O(h) formulas for comparison
        # -----------------------------------------------------

        basic_first = (
            y_values[1]
            - y_values[0]
        ) / h

        basic_second = (
            y_values[2]
            - 2 * y_values[1]
            + y_values[0]
        ) / (h ** 2)

        method_name = "Forward"

    # =========================================================
    # 8. BACKWARD DIFFERENCE AT xn
    # =========================================================

    elif endpoint_choice == "2":

        i = len(x_values) - 1

        x = x_values[i]

        # -----------------------------------------------------
        # First derivative: O(h^2)
        #
        # y'_n =
        # (3y_n - 4y_(n-1) + y_(n-2)) / (2h)
        # -----------------------------------------------------

        higher_first = (
            3 * y_values[i]
            - 4 * y_values[i - 1]
            + y_values[i - 2]
        ) / (2 * h)

        # -----------------------------------------------------
        # Second derivative: O(h^2)
        #
        # y''_n =
        # (2y_n - 5y_(n-1)
        #  + 4y_(n-2) - y_(n-3)) / h^2
        # -----------------------------------------------------

        higher_second = (
            2 * y_values[i]
            - 5 * y_values[i - 1]
            + 4 * y_values[i - 2]
            - y_values[i - 3]
        ) / (h ** 2)

        # -----------------------------------------------------
        # Basic O(h) formulas for comparison
        # -----------------------------------------------------

        basic_first = (
            y_values[i]
            - y_values[i - 1]
        ) / h

        basic_second = (
            y_values[i]
            - 2 * y_values[i - 1]
            + y_values[i - 2]
        ) / (h ** 2)

        method_name = "Backward"

    else:

        print("\nInvalid option.")

        return

    # =========================================================
    # 9. DISPLAY RESULTS
    # =========================================================

    print("\nHIGHER-ORDER RESULTS")
    print("-" * 75)

    print(
        f"Method = {method_name} difference"
    )

    print(
        f"Order = O(h^2)"
    )

    print(
        f"x = {x:.10f}"
    )

    print(
        f"h = {h:.10f}"
    )

    # =========================================================
    # 10. DISPLAY FIRST AND SECOND DERIVATIVE APPROXIMATIONS
    # =========================================================

    print("\nAPPROXIMATE DERIVATIVES")
    print("-" * 75)

    print(
        f"First derivative  = "
        f"{higher_first:.10f}"
    )

    print(
        f"Second derivative = "
        f"{higher_second:.10f}"
    )

    # =========================================================
    # 11. DISPLAY EXACT VALUES AND ERRORS
    # =========================================================

    if exact_first_derivative is not None:

        exact_first = exact_first_derivative(x)

        exact_second = exact_second_derivative(x)

        first_absolute_error = abs(
            higher_first - exact_first
        )

        second_absolute_error = abs(
            higher_second - exact_second
        )

        if exact_first != 0:

            first_percent_error = (
                first_absolute_error
                / abs(exact_first)
            ) * 100

        else:

            first_percent_error = None

        if exact_second != 0:

            second_percent_error = (
                second_absolute_error
                / abs(exact_second)
            ) * 100

        else:

            second_percent_error = None

        print("\nEXACT VALUES")
        print("-" * 75)

        print(
            f"Exact first derivative  = "
            f"{exact_first:.10f}"
        )

        print(
            f"Exact second derivative = "
            f"{exact_second:.10f}"
        )

        print("\nHIGHER-ORDER ABSOLUTE ERRORS")
        print("-" * 75)

        print(
            f"First derivative error  = "
            f"{first_absolute_error:.10f}"
        )

        print(
            f"Second derivative error = "
            f"{second_absolute_error:.10f}"
        )

        print("\nHIGHER-ORDER PERCENTAGE ERRORS")
        print("-" * 75)

        if first_percent_error is not None:

            print(
                f"First derivative % error  = "
                f"{first_percent_error:.6f}%"
            )

        else:

            print(
                "First derivative % error  = N/A"
            )

        if second_percent_error is not None:

            print(
                f"Second derivative % error = "
                f"{second_percent_error:.6f}%"
            )

        else:

            print(
                "Second derivative % error = N/A"
            )

        # =====================================================
        # 12. COMPARE O(h) AND O(h^2)
        # =====================================================

        basic_first_error = abs(
            basic_first - exact_first
        )

        basic_second_error = abs(
            basic_second - exact_second
        )

        if exact_first != 0:

            basic_first_percent = (
                basic_first_error
                / abs(exact_first)
            ) * 100

        else:

            basic_first_percent = None

        if exact_second != 0:

            basic_second_percent = (
                basic_second_error
                / abs(exact_second)
            ) * 100

        else:

            basic_second_percent = None

        print("\nCOMPARISON WITH BASIC O(h) FORMULAS")
        print("-" * 85)

        print(
            f"{'Derivative':<20}"
            f"{'Basic O(h)':>18}"
            f"{'Basic Error %':>18}"
            f"{'Higher O(h^2)':>20}"
            f"{'Higher Error %':>18}"
        )

        print("-" * 85)

        basic_first_error_display = (
            f"{basic_first_percent:.6f}%"
            if basic_first_percent is not None
            else "N/A"
        )

        higher_first_error_display = (
            f"{first_percent_error:.6f}%"
            if first_percent_error is not None
            else "N/A"
        )

        basic_second_error_display = (
            f"{basic_second_percent:.6f}%"
            if basic_second_percent is not None
            else "N/A"
        )

        higher_second_error_display = (
            f"{second_percent_error:.6f}%"
            if second_percent_error is not None
            else "N/A"
        )

        print(
            f"{'First derivative':<20}"
            f"{basic_first:>18.8f}"
            f"{basic_first_error_display:>18}"
            f"{higher_first:>20.8f}"
            f"{higher_first_error_display:>18}"
        )

        print(
            f"{'Second derivative':<20}"
            f"{basic_second:>18.8f}"
            f"{basic_second_error_display:>18}"
            f"{higher_second:>20.8f}"
            f"{higher_second_error_display:>18}"
        )

    else:

        print(
            "\nNo analytical function was supplied."
        )

        print(
            "Exact values and errors cannot be calculated."
        )

# ============================================================
# TRAPEZOIDAL RULE
# ============================================================

def trapezoidal_rule():
    """
    OPTION 19
    Composite Trapezoidal Rule

    I ≈ h/2 [y0 + 2(y1 + y2 + ... + y_(n-1)) + yn]
    """

    print("\nTRAPEZOIDAL RULE")
    print("-" * 75)

    print("\nChoose input type:")
    print("1. Function f(x)")
    print("2. Tabulated data")

    input_type = input("Choose an option: ")

    # =========================================================
    # FUNCTION INPUT
    # =========================================================

    if input_type == "1":

        function_input = input(
            "\nEnter f(x): "
        )

        a = get_number(
            "Enter lower limit a: "
        )

        b = get_number(
            "Enter upper limit b: "
        )

        n = int(get_number(
            "Enter number of subintervals n: "
        ))

        if n <= 0:
            print("Error: n must be positive.")
            return

        try:

            expression = convert_expression(
                function_input
            )

            f = sp.lambdify(
                x_symbol,
                expression,
                "math"
            )

        except Exception as error:

            print("Error interpreting f(x):")
            print(error)
            return

        h = (b - a) / n

        x_values = [
            a + i * h
            for i in range(n + 1)
        ]

        y_values = [
            f(x)
            for x in x_values
        ]

    # =========================================================
    # TABULATED DATA
    # =========================================================

    elif input_type == "2":

        x_input = input(
            "\nEnter x-values separated by commas: "
        )

        y_input = input(
            "Enter y-values separated by commas: "
        )

        try:

            x_values = [
                float(v.strip())
                for v in x_input.split(",")
            ]

            y_values = [
                float(v.strip())
                for v in y_input.split(",")
            ]

        except ValueError:

            print("Error: Invalid numerical data.")
            return

        if len(x_values) != len(y_values):

            print(
                "Error: x and y must contain "
                "the same number of values."
            )

            return

        if len(x_values) < 2:

            print(
                "Error: At least two data points are required."
            )

            return

        h_values = [
            x_values[i + 1] - x_values[i]
            for i in range(len(x_values) - 1)
        ]

        h = h_values[0]
        if h == 0:
            print("Error: x-values must be distinct.")
            return

        for value in h_values:

            if abs(value - h) > 1e-10:

                # A-2 also allows unequal intervals for the trapezoidal rule.
                if not (all(v > 0 for v in h_values) or all(v < 0 for v in h_values)):
                    print("Error: x-values must be strictly ordered.")
                    return
                areas = [h_values[i] * (y_values[i] + y_values[i+1]) / 2
                         for i in range(len(h_values))]
                integral = sum(areas)
                print("\nUnequal-spacing trapezoids: A_i = h_i*(y_i + y_(i+1))/2")
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
                return

        n = len(x_values) - 1

    else:

        print("Invalid option.")
        return

    # =========================================================
    # TRAPEZOIDAL CALCULATION
    # =========================================================

    weights = [1] + [2] * (n - 1) + [1]

    weighted_values = [
        weights[i] * y_values[i]
        for i in range(n + 1)
    ]

    weighted_sum = sum(weighted_values)

    integral = (h / 2) * weighted_sum

    # =========================================================
    # RESULTS
    # =========================================================

    print("\nTRAPEZOIDAL RULE RESULTS")
    print("-" * 90)

    print(
        f"{'i':>5}"
        f"{'x_i':>15}"
        f"{'y_i':>18}"
        f"{'w_i':>10}"
        f"{'w_i*y_i':>18}"
    )

    print("-" * 90)

    for i in range(n + 1):

        print(
            f"{i:>5}"
            f"{x_values[i]:>15.6f}"
            f"{y_values[i]:>18.8f}"
            f"{weights[i]:>10}"
            f"{weighted_values[i]:>18.8f}"
        )

    print("-" * 90)

    print(
        f"Sum(w_i*y_i) = {weighted_sum:.10f}"
    )

    print(
        f"\nh = {h:.10f}"
    )

    print(
        f"Integral ≈ {integral:.10f}"
    )

    # =========================================================
    # OPTIONAL EXACT VALUE
    # =========================================================

    exact_input = input(
        "\nEnter exact integral value "
        "(press ENTER if unavailable): "
    )

    if exact_input.strip() != "":

        try:

            exact_value = float(
                exact_input
            )

            error = abs(
                integral - exact_value
            )

            error_percent = (
                error / abs(exact_value)
            ) * 100 if exact_value != 0 else float("nan")

            print(
                f"\nExact integral = "
                f"{exact_value:.10f}"
            )

            print(
                f"Absolute error = "
                f"{error:.10f}"
            )

            print(
                f"Absolute percentage error = "
                f"{error_percent:.6f}%"
            )

        except ValueError:

            print(
                "Invalid exact value."
            )

# ============================================================
# SIMPSON'S 1/3 RULE
# ============================================================

def simpson_one_third_rule():
    """
    OPTION 20
    Composite Simpson's 1/3 Rule

    I ≈ h/3 [
        y0
        + 4y1
        + 2y2
        + 4y3
        + ...
        + 4y_(n-1)
        + yn
    ]

    n must be even.
    """

    print("\nSIMPSON'S 1/3 RULE")
    print("-" * 75)

    print("\nChoose input type:")
    print("1. Function f(x)")
    print("2. Tabulated data")

    input_type = input("Choose an option: ")

    # =========================================================
    # FUNCTION INPUT
    # =========================================================

    if input_type == "1":

        function_input = input(
            "\nEnter f(x): "
        )

        a = get_number(
            "Enter lower limit a: "
        )

        b = get_number(
            "Enter upper limit b: "
        )

        n = int(get_number(
            "Enter number of subintervals n: "
        ))

        if n <= 0:

            print("Error: n must be positive.")
            return

        if n % 2 != 0:

            print(
                "Error: Simpson's 1/3 rule requires "
                "an even value of n."
            )

            return

        try:

            expression = convert_expression(
                function_input
            )

            f = sp.lambdify(
                x_symbol,
                expression,
                "math"
            )

        except Exception as error:

            print("Error interpreting f(x):")
            print(error)
            return

        h = (b - a) / n

        x_values = [
            a + i * h
            for i in range(n + 1)
        ]

        y_values = [
            f(x)
            for x in x_values
        ]

    # =========================================================
    # TABULATED DATA
    # =========================================================

    elif input_type == "2":

        x_input = input(
            "\nEnter x-values separated by commas: "
        )

        y_input = input(
            "Enter y-values separated by commas: "
        )

        try:

            x_values = [
                float(v.strip())
                for v in x_input.split(",")
            ]

            y_values = [
                float(v.strip())
                for v in y_input.split(",")
            ]

        except ValueError:

            print("Error: Invalid numerical data.")
            return

        if len(x_values) != len(y_values):

            print(
                "Error: x and y must contain "
                "the same number of values."
            )

            return

        n = len(x_values) - 1
        if n < 2:
            print("Error: Simpson's rule requires at least three data points.")
            return

        if n % 2 != 0:

            print(
                f"\nError: Simpson's 1/3 rule requires "
                f"an even n."
            )

            print(
                f"Your data gives n = {n}."
            )

            return

        h_values = [
            x_values[i + 1] - x_values[i]
            for i in range(n)
        ]

        h = h_values[0]
        if h == 0:
            print("Error: x-values must be distinct.")
            return

        for value in h_values:

            if abs(value - h) > 1e-10:

                print(
                    "Error: x-values must have uniform spacing."
                )

                return

    else:

        print("Invalid option.")
        return

    # =========================================================
    # SIMPSON WEIGHTS
    # =========================================================

    weights = []

    for i in range(n + 1):

        if i == 0 or i == n:

            weights.append(1)

        elif i % 2 == 1:

            weights.append(4)

        else:

            weights.append(2)

    weighted_values = [
        weights[i] * y_values[i]
        for i in range(n + 1)
    ]

    weighted_sum = sum(weighted_values)

    integral = (h / 3) * weighted_sum

    # =========================================================
    # RESULTS
    # =========================================================

    print("\nSIMPSON'S 1/3 RULE RESULTS")
    print("-" * 90)

    print(
        f"{'i':>5}"
        f"{'x_i':>15}"
        f"{'y_i':>18}"
        f"{'w_i':>10}"
        f"{'w_i*y_i':>18}"
    )

    print("-" * 90)

    for i in range(n + 1):

        print(
            f"{i:>5}"
            f"{x_values[i]:>15.6f}"
            f"{y_values[i]:>18.8f}"
            f"{weights[i]:>10}"
            f"{weighted_values[i]:>18.8f}"
        )

    print("-" * 90)

    print(
        f"Sum(w_i*y_i) = {weighted_sum:.10f}"
    )

    print(
        f"\nh = {h:.10f}"
    )

    print(
        f"n = {n}"
    )

    print(
        f"Integral ≈ {integral:.10f}"
    )

    # =========================================================
    # OPTIONAL EXACT VALUE
    # =========================================================

    exact_input = input(
        "\nEnter exact integral value "
        "(press ENTER if unavailable): "
    )

    if exact_input.strip() != "":

        try:

            exact_value = float(
                exact_input
            )

            error = abs(
                integral - exact_value
            )

            error_percent = (
                error / abs(exact_value)
            ) * 100 if exact_value != 0 else float("nan")

            print(
                f"\nExact integral = "
                f"{exact_value:.10f}"
            )

            print(
                f"Absolute error = "
                f"{error:.10f}"
            )

            print(
                f"Absolute percentage error = "
                f"{error_percent:.6f}%"
            )

        except ValueError:

            print(
                "Invalid exact value."
            )

# ============================================================
# COMPARISON TRAPEZOIDAL AND SIMPSON'S
# ============================================================

def compare_integration_methods():
    """
    OPTION 21

    Compare the composite Trapezoidal Rule and
    Simpson's 1/3 Rule for the same function.

    Simpson's rule requires even n.
    """

    print("\nCOMPARISON OF NUMERICAL INTEGRATION METHODS")
    print("-" * 75)

    # =========================================================
    # INPUT FUNCTION
    # =========================================================

    function_input = input(
        "Enter f(x): "
    )

    a = get_number(
        "Enter lower limit a: "
    )

    b = get_number(
        "Enter upper limit b: "
    )

    n_trap = int(get_number(
        "Enter n for Trapezoidal Rule: "
    ))

    n_simp = int(get_number(
        "Enter n for Simpson's 1/3 Rule "
        "(must be even): "
    ))

    if n_trap <= 0:

        print("Error: Trapezoidal n must be positive.")
        return

    if n_simp <= 0 or n_simp % 2 != 0:

        print(
            "Error: Simpson's n must be positive and even."
        )

        return

    try:

        expression = convert_expression(
            function_input
        )

        f = sp.lambdify(
            x_symbol,
            expression,
            "math"
        )

    except Exception as error:

        print("Error interpreting f(x):")
        print(error)
        return

    # =========================================================
    # TRAPEZOIDAL CALCULATION
    # =========================================================

    h_trap = (b - a) / n_trap

    x_trap = [
        a + i * h_trap
        for i in range(n_trap + 1)
    ]

    y_trap = [
        f(x)
        for x in x_trap
    ]

    trap_sum = (
        y_trap[0]
        + y_trap[-1]
        + 2 * sum(y_trap[1:-1])
    )

    trapezoidal_result = (
        h_trap / 2
    ) * trap_sum

    # =========================================================
    # SIMPSON CALCULATION
    # =========================================================

    h_simp = (b - a) / n_simp

    x_simp = [
        a + i * h_simp
        for i in range(n_simp + 1)
    ]

    y_simp = [
        f(x)
        for x in x_simp
    ]

    simpson_sum = (
        y_simp[0]
        + y_simp[-1]
    )

    for i in range(1, n_simp):

        if i % 2 == 1:

            simpson_sum += (
                4 * y_simp[i]
            )

        else:

            simpson_sum += (
                2 * y_simp[i]
            )

    simpson_result = (
        h_simp / 3
    ) * simpson_sum

    # =========================================================
    # DISPLAY RESULTS
    # =========================================================

    print("\nRESULTS")
    print("-" * 75)

    print(
        f"Trapezoidal Rule:"
    )

    print(
        f"    n = {n_trap}"
    )

    print(
        f"    h = {h_trap:.10f}"
    )

    print(
        f"    Integral = "
        f"{trapezoidal_result:.10f}"
    )

    print(
        f"\nSimpson's 1/3 Rule:"
    )

    print(
        f"    n = {n_simp}"
    )

    print(
        f"    h = {h_simp:.10f}"
    )

    print(
        f"    Integral = "
        f"{simpson_result:.10f}"
    )

    # =========================================================
    # OPTIONAL EXACT VALUE
    # =========================================================

    exact_input = input(
        "\nEnter exact integral value "
        "(press ENTER if unavailable): "
    )

    if exact_input.strip() != "":

        try:

            exact_value = float(
                exact_input
            )

            trap_error = abs(
                trapezoidal_result
                - exact_value
            )

            simp_error = abs(
                simpson_result
                - exact_value
            )

            trap_error_percent = (
                trap_error
                / abs(exact_value)
            ) * 100 if exact_value != 0 else float("nan")

            simp_error_percent = (
                simp_error
                / abs(exact_value)
            ) * 100 if exact_value != 0 else float("nan")

            print("\nERROR COMPARISON")
            print("-" * 75)

            print(
                f"Exact integral = "
                f"{exact_value:.10f}"
            )

            print(
                f"\nTrapezoidal absolute error = "
                f"{trap_error:.10f}"
            )

            print(
                f"Trapezoidal percentage error = "
                f"{trap_error_percent:.6f}%"
            )

            print(
                f"\nSimpson absolute error = "
                f"{simp_error:.10f}"
            )

            print(
                f"Simpson percentage error = "
                f"{simp_error_percent:.6f}%"
            )

            # -------------------------------------------------
            # Accuracy comparison
            # -------------------------------------------------

            print("\nACCURACY COMPARISON")
            print("-" * 75)

            if trap_error < simp_error:

                print(
                    "Trapezoidal Rule produced the "
                    "smaller absolute error."
                )

            elif simp_error < trap_error:

                print(
                    "Simpson's 1/3 Rule produced the "
                    "smaller absolute error."
                )

            else:

                print(
                    "Both methods produced the same "
                    "absolute error."
                )

        except ValueError:

            print(
                "Invalid exact value."
            )


# ============================================================
# 27. VANDERMONDE POLYNOMIAL APPROXIMATION
# ============================================================

def vandermonde_approximation():
    saved_results.pop("vandermonde", None)
    print("\nVANDERMONDE POLYNOMIAL APPROXIMATION")
    print("-" * 70)

    # --------------------------------------------------------
    # User input
    # --------------------------------------------------------

    n = int(input("Enter number of data points: "))

    x_values = []
    y_values = []

    print("\nEnter the data points:")

    for i in range(n):
        x_i = get_number(f"x[{i}] = ")
        y_i = get_number(f"y[{i}] = ")

        x_values.append(x_i)
        y_values.append(y_i)

    if n < 1 or len(set(x_values)) != n:
        print("At least one point and distinct x-values are required.")
        return

    # --------------------------------------------------------
    # Construct Vandermonde matrix
    #
    # [1  x1  x1^2  ... x1^(n-1)]
    # [1  x2  x2^2  ... x2^(n-1)]
    # [        ...              ]
    # [1  xn  xn^2  ... xn^(n-1)]
    # --------------------------------------------------------

    A = np.array([
        [x_values[i] ** j for j in range(n)]
        for i in range(n)
    ], dtype=float)

    b = np.array(y_values, dtype=float)

    # --------------------------------------------------------
    # Solve A a = b
    # --------------------------------------------------------

    try:
        coefficients = np.linalg.solve(A, b)

    except np.linalg.LinAlgError:
        print("\nERROR:")
        print("The Vandermonde matrix is singular.")
        print("Make sure all x-values are different.")
        return

    # --------------------------------------------------------
    # Display data
    # --------------------------------------------------------

    print("\nDATA POINTS")
    print("-" * 70)
    print(f"{'i':>5} {'x_i':>15} {'y_i':>15}")

    for i in range(n):
        print(
            f"{i:>5} "
            f"{x_values[i]:>15.6f} "
            f"{y_values[i]:>15.6f}"
        )

    # --------------------------------------------------------
    # Display Vandermonde matrix
    # --------------------------------------------------------

    print("\nVANDERMONDE MATRIX")
    print("-" * 70)

    for row in A:
        print(" ".join(f"{value:12.6f}" for value in row))

    # --------------------------------------------------------
    # Display coefficients
    # --------------------------------------------------------

    print("\nCOEFFICIENTS")
    print("-" * 70)

    for i, coefficient in enumerate(coefficients):
        print(f"a{i + 1} = {coefficient:.10f}")

    # --------------------------------------------------------
    # Construct polynomial
    #
    # f(x) = a1 + a2*x + a3*x^2 + ... + an*x^(n-1)
    # --------------------------------------------------------

    x = sp.symbols("x")

    polynomial = sum(
        coefficients[j] * x**j
        for j in range(n)
    )

    polynomial = sp.expand(polynomial)

    print("\nAPPROXIMATING POLYNOMIAL")
    print("-" * 70)

    print("f(x) =", polynomial)

    # --------------------------------------------------------
    # Verify that polynomial passes through the data points
    # --------------------------------------------------------

    print("\nVERIFICATION")
    print("-" * 70)
    print(f"{'x':>12} {'given y':>15} {'f(x)':>15} {'error':>15}")

    polynomial_function = sp.lambdify(x, polynomial, "math")

    for i in range(n):
        calculated_y = polynomial_function(x_values[i])
        error = abs(y_values[i] - calculated_y)

        print(
            f"{x_values[i]:>12.6f} "
            f"{y_values[i]:>15.6f} "
            f"{calculated_y:>15.6f} "
            f"{error:>15.6e}"
        )

    # --------------------------------------------------------
    # Optional evaluation
    # --------------------------------------------------------

    evaluate = input(
        "\nEvaluate the approximating polynomial at a value of x? (y/n): "
    ).lower()

    if evaluate == "y":
        x_eval = get_number("Enter x-value: ")
        y_eval = polynomial_function(x_eval)

        print(
            f"\nf({x_eval}) = {y_eval:.10f}"
        )

    # Retain the full-precision polynomial, not the rounded printed formula.
    saved_results["vandermonde"] = {
        "x": list(x_values), "y": list(y_values), "polynomial": polynomial,
        "function": sp.lambdify(x, polynomial, "math")
    }



# ============================================================
# 28. LAGRANGIAN POLYNOMIAL
# ============================================================

def lagrangian_polynomial():
    saved_results.pop("lagrangian", None)
    print("\nLAGRANGIAN POLYNOMIAL")
    print("-" * 70)

    n = int(input("Enter number of data points: "))

    x_values = []
    y_values = []

    print("\nEnter the data points:")

    for i in range(n):
        x_i = get_number(f"x[{i}] = ")
        y_i = get_number(f"y[{i}] = ")

        x_values.append(x_i)
        y_values.append(y_i)

    if n < 1 or len(set(x_values)) != n:
        print("At least one point and distinct x-values are required.")
        return

    x = sp.symbols("x")

    # --------------------------------------------------------
    # Construct Lagrangian basis polynomials
    #
    # phi_i(x) =
    # product[(x-x_j)/(x_i-x_j)], j != i
    # --------------------------------------------------------

    lagrange_basis = []

    for i in range(n):

        phi_i = 1

        for j in range(n):

            if i != j:
                phi_i *= (x - x_values[j]) / (
                    x_values[i] - x_values[j]
                )

        lagrange_basis.append(
            sp.expand(phi_i)
        )

    # --------------------------------------------------------
    # Display Lagrangian basis functions
    # --------------------------------------------------------

    print("\nLAGRANGIAN BASIS POLYNOMIALS")
    print("-" * 70)

    for i, phi_i in enumerate(lagrange_basis):

        print(f"\nphi_{i}(x) =")
        print(sp.N(phi_i, 10))

    # --------------------------------------------------------
    # Construct interpolation polynomial
    #
    # f(x) = sum y_i phi_i(x)
    # --------------------------------------------------------

    polynomial = sum(
        y_values[i] * lagrange_basis[i]
        for i in range(n)
    )

    polynomial = sp.expand(polynomial)

    print("\nAPPROXIMATION FUNCTION")
    print("-" * 70)

    print("f(x) =")
    print(sp.N(polynomial, 10))

    # --------------------------------------------------------
    # First and second derivatives
    # --------------------------------------------------------

    first_derivative = sp.diff(polynomial, x)
    second_derivative = sp.diff(polynomial, x, 2)

    print("\nFIRST DERIVATIVE")
    print("-" * 70)

    print("f'(x) =")
    print(sp.N(first_derivative, 10))

    print("\nSECOND DERIVATIVE")
    print("-" * 70)

    print("f''(x) =")
    print(sp.N(second_derivative, 10))

    # --------------------------------------------------------
    # Optional evaluation
    # --------------------------------------------------------

    evaluate = input(
        "\nEvaluate the Lagrangian polynomial? (y/n): "
    ).lower()

    if evaluate == "y":

        x_eval = get_number("Enter x-value: ")

        polynomial_function = sp.lambdify(
            x,
            polynomial,
            "math"
        )

        y_eval = polynomial_function(x_eval)

        print(
            f"\nf({x_eval}) = {y_eval:.10f}"
        )
    # Retain the full-precision polynomial, not the rounded printed formula.
    saved_results["lagrangian"] = {
        "x": list(x_values), "y": list(y_values), "polynomial": polynomial,
        "function": sp.lambdify(x, polynomial, "math")
    }



# ============================================================
# 29. LAGRANGIAN AND FINITE DIFFERENCE COMPARISON
# ============================================================

def lagrangian_finite_difference():
    saved_results.pop("finite_difference", None)
    print("\nLAGRANGIAN AND FINITE DIFFERENCE COMPARISON")
    print("-" * 70)

    n = int(input("Enter number of data points: "))

    x_values = []
    y_values = []

    print("\nEnter the data points:")

    for i in range(n):

        x_i = get_number(f"x[{i}] = ")
        y_i = get_number(f"y[{i}] = ")

        x_values.append(x_i)
        y_values.append(y_i)

    if n < 1 or len(set(x_values)) != n:
        print("At least one point and distinct x-values are required.")
        return

    # --------------------------------------------------------
    # Check that x-values are equally spaced
    # --------------------------------------------------------

    if n < 3:
        print("\nAt least 3 data points are required.")
        return

    h = x_values[1] - x_values[0]

    equally_spaced = all(
        abs((x_values[i] - x_values[i - 1]) - h) < 1e-10
        for i in range(2, n)
    )

    if not equally_spaced:

        print("\nERROR:")
        print(
            "Finite-difference formulas require equally spaced x-values."
        )
        return

    # --------------------------------------------------------
    # Construct Lagrangian polynomial
    # --------------------------------------------------------

    x = sp.symbols("x")

    lagrange_basis = []

    for i in range(n):

        phi_i = 1

        for j in range(n):

            if i != j:

                phi_i *= (
                    (x - x_values[j])
                    /
                    (x_values[i] - x_values[j])
                )

        lagrange_basis.append(
            sp.expand(phi_i)
        )

    polynomial = sum(
        y_values[i] * lagrange_basis[i]
        for i in range(n)
    )

    polynomial = sp.expand(polynomial)

    # --------------------------------------------------------
    # Analytical derivatives of Lagrangian polynomial
    # --------------------------------------------------------

    lagrange_first = sp.diff(
        polynomial,
        x
    )

    lagrange_second = sp.diff(
        polynomial,
        x,
        2
    )

    fprime = sp.lambdify(
        x,
        lagrange_first,
        "math"
    )

    fdoubleprime = sp.lambdify(
        x,
        lagrange_second,
        "math"
    )

    # --------------------------------------------------------
    # Display Lagrangian polynomial
    # --------------------------------------------------------

    print("\nLAGRANGIAN POLYNOMIAL")
    print("-" * 70)

    print(
        sp.N(polynomial, 10)
    )

    print("\nLAGRANGIAN FIRST DERIVATIVE")
    print("-" * 70)

    print(
        sp.N(lagrange_first, 10)
    )

    print("\nLAGRANGIAN SECOND DERIVATIVE")
    print("-" * 70

    )

    print(
        sp.N(lagrange_second, 10)
    )

    # --------------------------------------------------------
    # Finite-difference approximations
    #
    # First derivative:
    #
    # i = 0:
    # (-3y_i + 4y_i+1 - y_i+2) / (2h)
    #
    # interior:
    # (y_i+1 - y_i-1) / (2h)
    #
    # i = n-1:
    # (3y_i - 4y_i-1 + y_i-2) / (2h)
    #
    # Second derivative:
    #
    # interior:
    # (y_i-1 - 2y_i + y_i+1) / h^2
    # --------------------------------------------------------

    fd_first = []
    fd_second = []

    for i in range(n):

        # First derivative
        if i == 0:

            first = (
                -3 * y_values[i]
                + 4 * y_values[i + 1]
                - y_values[i + 2]
            ) / (2 * h)

        elif i == n - 1:

            first = (
                3 * y_values[i]
                - 4 * y_values[i - 1]
                + y_values[i - 2]
            ) / (2 * h)

        else:

            first = (
                y_values[i + 1]
                - y_values[i - 1]
            ) / (2 * h)

        fd_first.append(first)

        # Second derivative
        if i == 0:

            if n >= 4:

                second = (
                    2 * y_values[i]
                    - 5 * y_values[i + 1]
                    + 4 * y_values[i + 2]
                    - y_values[i + 3]
                ) / (h ** 2)

            else:

                second = None

        elif i == n - 1:

            if n >= 4:

                second = (
                    2 * y_values[i]
                    - 5 * y_values[i - 1]
                    + 4 * y_values[i - 2]
                    - y_values[i - 3]
                ) / (h ** 2)

            else:

                second = None

        else:

            second = (
                y_values[i - 1]
                - 2 * y_values[i]
                + y_values[i + 1]
            ) / (h ** 2)

        fd_second.append(second)

    # --------------------------------------------------------
    # Comparison table
    # --------------------------------------------------------

    print("\nDERIVATIVE COMPARISON")
    print("-" * 100)

    print(
        f"{'i':>3} "
        f"{'x':>12} "
        f"{'Lagrange f_prime':>20} "
        f"{'FD f_prime':>18} "
        f"{'Lagrange f_double_prime':>25} "
        f"{'FD f_double_prime':>20}"
    )

    print("-" * 100)

    for i in range(n):

        lag1 = fprime(x_values[i])
        fd1 = fd_first[i]

        lag2 = fdoubleprime(x_values[i])

        if fd_second[i] is None:
            fd2_text = "N/A"
        else:
            fd2_text = f"{fd_second[i]:.10f}"

        print(
            f"{i:>3} "
            f"{x_values[i]:>12.6f} "
            f"{lag1:>20.10f} "
            f"{fd1:>18.10f} "
            f"{lag2:>25.10f} "
            f"{fd2_text:>20}"
        )

    # --------------------------------------------------------
    # Plot comparison
    # --------------------------------------------------------

    saved_results["finite_difference"] = {
        "x": list(x_values), "y": list(y_values),
        "first": list(fd_first), "second": list(fd_second),
        "lagrange_first": fprime, "lagrange_second": fdoubleprime
    }
    print("Results saved. Use option 36 to plot the derivative comparisons.")


# ============================================================
# 30. LEAST SQUARES
# ============================================================

def least_squares():
    saved_results.pop("least_squares", None)
    print("\nLEAST SQUARES METHOD")
    print("=" * 70)

    # --------------------------------------------------------
    # INPUT DATA
    # --------------------------------------------------------

    n = int(input("Enter number of data points: "))

    x_data = []
    y_data = []

    print("\nEnter the data points:")

    for i in range(n):
        x_i = get_number(f"x[{i}] = ")
        y_i = get_number(f"y[{i}] = ")

        x_data.append(x_i)
        y_data.append(y_i)

    x_data = np.array(x_data, dtype=float)
    y_data = np.array(y_data, dtype=float)

    print("\nDATA")
    print("-" * 70)

    for i in range(n):
        print(f"{i:<5}{x_data[i]:>12.6f}{y_data[i]:>15.6f}")

    # --------------------------------------------------------
    # SELECT CANDIDATE FUNCTIONS
    # --------------------------------------------------------

    print("\nCANDIDATE FUNCTIONS")
    print("-" * 70)
    print("1. Linear        f(x) = a1 + a2*x")
    print("2. Quadratic     f(x) = a1 + a2*x + a3*x^2")
    print("3. Polynomial    f(x) = a1 + a2*x + ... + ap*x^(p-1)")
    print("4. Exponential   f(x) = a*exp(b*x)")
    print("5. All")

    choice = input("\nChoose function(s): ").strip()

    # --------------------------------------------------------
    # STORAGE
    # --------------------------------------------------------

    results = []

    # --------------------------------------------------------
    # POLYNOMIAL LEAST-SQUARES FUNCTION
    # --------------------------------------------------------

    def polynomial_fit(degree):

        p = degree + 1

        if n <= degree:
            print(
                f"\nCannot fit a degree-{degree} polynomial "
                f"with only {n} data points."
            )
            return None

        # Construct S matrix
        #
        # S =
        # [1  x1  x1^2  ... x1^degree]
        # [1  x2  x2^2  ... x2^degree]
        # [ .   .   .          .      ]
        # [1  xn  xn^2  ... xn^degree]

        S = np.column_stack(
            [x_data ** j for j in range(p)]
        )

        # Normal equations:
        #
        # [S]^T [S] {a} = [S]^T {y}

        STS = S.T @ S
        STy = S.T @ y_data

        try:
            coefficients = np.linalg.solve(STS, STy)
        except np.linalg.LinAlgError:
            print("\nUnable to solve the normal equations.")
            print("The matrix may be singular or ill-conditioned.")
            return None

        # Predicted values
        y_fit = S @ coefficients

        # Loss function
        #
        # L = sum [f(x_i) - y_i]^2

        loss = np.sum((y_fit - y_data) ** 2)

        # ----------------------------------------------------
        # DISPLAY
        # ----------------------------------------------------

        print("\n" + "=" * 70)
        print(f"POLYNOMIAL LEAST SQUARES — DEGREE {degree}")
        print("=" * 70)

        print("\nCoefficients:")

        for j, coefficient in enumerate(coefficients):
            print(f"a{j + 1} = {coefficient:.10f}")

        # Construct readable equation
        equation = f"{coefficients[0]:.10g}"

        for j in range(1, len(coefficients)):

            coefficient = coefficients[j]

            if coefficient >= 0:
                sign = " + "
            else:
                sign = " - "

            value = abs(coefficient)

            if j == 1:
                equation += f"{sign}{value:.10g}x"
            else:
                equation += f"{sign}{value:.10g}x^{j}"

        print("\nApproximation:")
        print(f"f(x) = {equation}")

        print(f"\nLoss L = {loss:.10g}")

        # ----------------------------------------------------
        # TABLE
        # ----------------------------------------------------

        print("\nRESULTS")
        print("-" * 70)
        print(
            f"{'i':<5}"
            f"{'x':>12}"
            f"{'y':>15}"
            f"{'f(x)':>15}"
            f"{'Residual':>15}"
        )
        print("-" * 70)

        for i in range(n):

            residual = y_fit[i] - y_data[i]

            print(
                f"{i:<5}"
                f"{x_data[i]:>12.6f}"
                f"{y_data[i]:>15.6f}"
                f"{y_fit[i]:>15.6f}"
                f"{residual:>15.6f}"
            )

        return {
            "name": f"Polynomial degree {degree}",
            "coefficients": coefficients,
            "function": lambda x, c=coefficients: sum(
                c[j] * x ** j
                for j in range(len(c))
            ),
            "loss": loss,
            "x": x_data,
            "y_fit": y_fit
        }

    # --------------------------------------------------------
    # EXPONENTIAL LEAST-SQUARES FUNCTION
    # --------------------------------------------------------

    def exponential_fit():

        # Exponential model:
        #
        # y = a exp(bx)
        #
        # ln(y) = ln(a) + bx
        #
        # Let:
        # A = ln(a)
        #
        # Then:
        # ln(y) = A + bx

        if np.any(y_data <= 0):

            print(
                "\nExponential fitting requires all y-values "
                "to be positive."
            )

            return None

        Y = np.log(y_data)

        # Construct linearized S matrix
        S = np.column_stack(
            [np.ones(n), x_data]
        )

        # Normal equations:
        #
        # [S]^T[S]{A,b} = [S]^T{ln(y)}

        STS = S.T @ S
        STY = S.T @ Y

        try:
            coefficients = np.linalg.solve(STS, STY)
        except np.linalg.LinAlgError:
            print("\nUnable to solve the exponential fit.")
            return None

        A = coefficients[0]
        b = coefficients[1]

        # Recover original parameter:
        #
        # a = exp(A)

        a = np.exp(A)

        # Original exponential function
        y_fit = a * np.exp(b * x_data)

        # Loss must be calculated using the ORIGINAL y-values
        loss = np.sum((y_fit - y_data) ** 2)

        # ----------------------------------------------------
        # DISPLAY
        # ----------------------------------------------------

        print("\n" + "=" * 70)
        print("EXPONENTIAL LEAST SQUARES")
        print("=" * 70)

        print("\nLinearized equation:")
        print(f"ln(y) = {A:.10f} + {b:.10f}x")

        print("\nRecovered parameters:")
        print(f"A = {A:.10f}")
        print(f"a = {a:.10f}")
        print(f"b = {b:.10f}")

        print("\nApproximation:")
        print(f"f(x) = {a:.10g} * exp({b:.10g}x)")

        print(f"\nLoss L = {loss:.10g}")

        # ----------------------------------------------------
        # TABLE
        # ----------------------------------------------------

        print("\nRESULTS")
        print("-" * 70)
        print(
            f"{'i':<5}"
            f"{'x':>12}"
            f"{'y':>15}"
            f"{'f(x)':>15}"
            f"{'Residual':>15}"
        )
        print("-" * 70)

        for i in range(n):

            residual = y_fit[i] - y_data[i]

            print(
                f"{i:<5}"
                f"{x_data[i]:>12.6f}"
                f"{y_data[i]:>15.6f}"
                f"{y_fit[i]:>15.6f}"
                f"{residual:>15.6f}"
            )

        return {
            "name": "Exponential",
            "a": a,
            "b": b,
            "function": lambda x, a=a, b=b: a * np.exp(b * x),
            "loss": loss,
            "x": x_data,
            "y_fit": y_fit
        }

    # --------------------------------------------------------
    # RUN SELECTED MODELS
    # --------------------------------------------------------

    if choice == "1":

        result = polynomial_fit(1)

        if result is not None:
            results.append(result)

    elif choice == "2":

        result = polynomial_fit(2)

        if result is not None:
            results.append(result)

    elif choice == "3":

        degree = int(
            input(
                "\nEnter polynomial degree "
                "(e.g. 3 for cubic): "
            )
        )

        result = polynomial_fit(degree)

        if result is not None:
            results.append(result)

    elif choice == "4":

        result = exponential_fit()

        if result is not None:
            results.append(result)

    elif choice == "5":

        print("\nRunning all candidate functions...")

        linear = polynomial_fit(1)

        if linear is not None:
            results.append(linear)

        quadratic = polynomial_fit(2)

        if quadratic is not None:
            results.append(quadratic)

        degree = int(
            input(
                "\nEnter additional polynomial degree "
                "(or enter 0 to skip): "
            )
        )

        if degree > 2:

            polynomial = polynomial_fit(degree)

            if polynomial is not None:
                results.append(polynomial)

        exponential = exponential_fit()

        if exponential is not None:
            results.append(exponential)

    else:

        print("\nInvalid choice.")
        return

    # --------------------------------------------------------
    # COMPARE LOSS VALUES
    # --------------------------------------------------------

    if len(results) > 1:

        print("\n" + "=" * 70)
        print("COMPARISON OF LEAST-SQUARES MODELS")
        print("=" * 70)

        print(
            f"{'Model':<30}"
            f"{'Loss L':>20}"
        )

        print("-" * 70)

        for result in results:

            print(
                f"{result['name']:<30}"
                f"{result['loss']:>20.10g}"
            )

        # Lowest loss = best fit

        best = min(
            results,
            key=lambda result: result["loss"]
        )

        print("\nBEST-FIT FUNCTION")
        print("-" * 70)
        print(f"Model: {best['name']}")
        print(f"Lowest loss: {best['loss']:.10g}")

        # ----------------------------------------------------
        # PLOT DATA AND FITTED FUNCTIONS
        # ----------------------------------------------------

    if results:
        saved_results["least_squares"] = {
            "x": x_data.copy(), "y": y_data.copy(), "models": results
        }
        print("Results saved. Use option 37 to plot the fitted models.")

# ============================================================
# OPTION 31 — BISECTION METHOD
# ============================================================

def bisection_method():
    saved_results.pop("root", None)
    saved_results.pop("bisection", None)
    print("\n" + "=" * 70)
    print("BISECTION METHOD")
    print("=" * 70)

    # --------------------------------------------------------
    # 1. INPUT FUNCTION
    # --------------------------------------------------------
    f_input = input("\nEnter f(x): ")

    try:
        f_expr = convert_expression(f_input)

        # Create numerical function of x
        f = sp.lambdify(x_symbol, f_expr, "math")

    except Exception as e:
        print("\nInvalid function.")
        print("Error:", e)
        return

    print("\nFunction:")
    print(f"f(x) = {formula_string(f_expr)}")

    # --------------------------------------------------------
    # 2. INPUT STARTING INTERVAL
    # --------------------------------------------------------
    a = get_number("\nEnter lower bound a: ")
    b = get_number("Enter upper bound b: ")

    # Make sure a < b
    if a > b:
        a, b = b, a

    # --------------------------------------------------------
    # 3. CHECK WHETHER THE INTERVAL BRACKETS A ROOT
    # --------------------------------------------------------
    try:
        fa = f(a)
        fb = f(b)
    except Exception as e:
        print("\nCould not evaluate the function at the interval endpoints.")
        print("Error:", e)
        return

    print("\nInitial interval:")
    print(f"[a, b] = [{a:.10f}, {b:.10f}]")
    print(f"f(a) = {fa:.10f}")
    print(f"f(b) = {fb:.10f}")

    if fa == 0:
        print(f"\nExact root found at a = {a:.10f}")
        result = {"method": "Bisection", "function": f, "root": a,
                  "iterations": [a], "converged": True}
        saved_results["bisection"] = saved_results["root"] = result
        return

    if fb == 0:
        print(f"\nExact root found at b = {b:.10f}")
        result = {"method": "Bisection", "function": f, "root": b,
                  "iterations": [b], "converged": True}
        saved_results["bisection"] = saved_results["root"] = result
        return

    if fa * fb > 0:
        print("\nERROR:")
        print("f(a) and f(b) have the same sign.")
        print("The interval does not bracket a root.")
        print("Bisection requires f(a) * f(b) < 0.")

        return

    # --------------------------------------------------------
    # 4. INPUT STOPPING CRITERION
    # --------------------------------------------------------
    tolerance = get_number(
        "\nEnter tolerance epsilon [default = 1e-4]: ",
        1e-4
    )

    max_iterations = int(
        get_number(
            "Enter maximum number of iterations [default = 100]: ",
            100
        )
    )

    if tolerance <= 0 or max_iterations < 1:
        print("Tolerance and maximum iterations must be positive.")
        return

    # --------------------------------------------------------
    # 5. BISECTION ITERATION
    # --------------------------------------------------------
    print("\n" + "=" * 110)
    print("BISECTION ITERATIONS")
    print("=" * 110)

    print(
        f"{'k':>3} "
        f"{'a^k':>14} "
        f"{'b^k':>14} "
        f"{'c^(k+1)':>14} "
        f"{'f(c)':>16} "
        f"{'|f(c)|':>16} "
        f"{'|f(c)| < eps':>14}"
    )

    print("-" * 110)

    rows = []

    root = None

    for k in range(max_iterations):

        # Midpoint
        c = (a + b) / 2

        fc = f(c)
        abs_fc = abs(fc)

        # Check stopping criterion
        converged = abs_fc < tolerance

        rows.append({
            "k": k,
            "a": a,
            "b": b,
            "c": c,
            "f_c": fc,
            "abs_f_c": abs_fc,
            "converged": converged
        })

        print(
            f"{k:3d} "
            f"{a:14.8f} "
            f"{b:14.8f} "
            f"{c:14.8f} "
            f"{fc:16.8f} "
            f"{abs_fc:16.8f} "
            f"{'YES' if converged else 'No':>14}"
        )

        # Stopping criterion
        if converged:
            root = c
            break

        # ----------------------------------------------------
        # Decide which half contains the root
        # ----------------------------------------------------
        if fa * fc < 0:
            # Root is in [a, c]
            b = c
            fb = fc

        else:
            # Root is in [c, b]
            a = c
            fa = fc

    # --------------------------------------------------------
    # 6. FINAL RESULT
    # --------------------------------------------------------
    print("\n" + "=" * 70)
    print("BISECTION RESULT")
    print("=" * 70)

    if root is not None:

        print(f"\nApproximate root = {root:.10f}")
        print(f"f(root)          = {f(root):.10f}")
        print(f"|f(root)|        = {abs(f(root)):.10f}")
        print(f"Tolerance        = {tolerance:.10f}")
        print(f"Iterations       = {len(rows)}")

    else:

        # If tolerance was not reached
        root = rows[-1]["c"]

        print("\nMaximum number of iterations reached.")
        print("The requested tolerance was not achieved.")

        print(f"\nLast approximation = {root:.10f}")
        print(f"f(c)               = {f(root):.10f}")
        print(f"|f(c)|             = {abs(f(root)):.10f}")
        print(f"Iterations         = {len(rows)}")


    result = {
        "method": "Bisection", "function": f,
        "root": root,
        "iterations": [row["c"] for row in rows], "converged": rows[-1]["converged"]
    }
    saved_results["bisection"] = saved_results["root"] = result
    print("Results saved. Options 38-39 use this root-finding run.")

# ============================================================
# OPTION 32 — NEWTON-RAPHSON METHOD
# ============================================================

def newton_raphson_method():
    saved_results.pop("root", None)
    saved_results.pop("newton", None)
    print("\n" + "=" * 75)
    print("NEWTON-RAPHSON METHOD")
    print("=" * 75)

    # --------------------------------------------------------
    # 1. INPUT FUNCTION
    # --------------------------------------------------------
    f_input = input("\nEnter f(x): ")

    try:
        # Convert typed formula into a SymPy expression
        f_expr = convert_expression(f_input)

        # Automatically calculate the first derivative
        df_expr = sp.diff(f_expr, x_symbol)

        # Convert symbolic expressions into numerical functions
        f = sp.lambdify(x_symbol, f_expr, "math")
        df = sp.lambdify(x_symbol, df_expr, "math")

    except Exception as e:
        print("\nInvalid function.")
        print("Error:", e)
        return

    # --------------------------------------------------------
    # 2. DISPLAY FUNCTION AND AUTOMATIC DERIVATIVE
    # --------------------------------------------------------
    print("\nFunction:")
    print(f"f(x)  = {formula_string(f_expr)}")

    print("\nAutomatically calculated derivative:")
    print(f"f'(x) = {formula_string(df_expr)}")

    # --------------------------------------------------------
    # 3. INPUT STARTING POINT
    # --------------------------------------------------------
    x_current = get_number(
        "\nEnter starting point x^0: "
    )

    # --------------------------------------------------------
    # 4. INPUT STOPPING CRITERION
    # --------------------------------------------------------
    tolerance = get_number(
        "\nEnter tolerance epsilon (%) [default = 1e-4]: ",
        1e-4
    )

    max_iterations = int(
        get_number(
            "Enter maximum number of iterations [default = 100]: ",
            100
        )
    )

    if tolerance <= 0 or max_iterations < 1:
        print("Tolerance and maximum iterations must be positive.")
        return

    # --------------------------------------------------------
    # 5. NEWTON-RAPHSON ITERATION
    # --------------------------------------------------------
    print("\n" + "=" * 115)
    print("NEWTON-RAPHSON ITERATIONS")
    print("=" * 115)

    print(
        f"{'k':>3} "
        f"{'x^k':>16} "
        f"{'f(x^k)':>18} "
        + "f'(x^k)".rjust(18) + " "
        + f"{'x^(k+1)':>16} "
        f"{'Approx. Error %':>20} "
        f"{'Converged?':>13}"
    )

    print("-" * 115)

    initial_x = x_current
    rows = []
    root = None

    for k in range(max_iterations):

        # Evaluate function and derivative at current x
        fx = f(x_current)
        if fx == 0:
            root = x_current
            rows.append({"error_percent": 0.0})
            break
        dfx = df(x_current)

        # ----------------------------------------------------
        # Check for zero derivative
        # ----------------------------------------------------
        if abs(dfx) < 1e-14:
            print(
                f"{k:3d} "
                f"{x_current:16.10f} "
                f"{fx:18.10f} "
                f"{dfx:18.10f}"
            )

            print("\nERROR:")
            print("The derivative is zero or extremely close to zero.")
            print("Newton-Raphson cannot continue from this point.")

            return

        # ----------------------------------------------------
        # Newton-Raphson equation
        #
        # x^(k+1) = x^k - f(x^k) / f'(x^k)
        # ----------------------------------------------------
        x_new = x_current - fx / dfx

        # ----------------------------------------------------
        # Lecture stopping criterion:
        #
        # |(x^(k+1) - x^k) / x^(k+1)| * 100%
        # ----------------------------------------------------
        if abs(x_new) > 1e-14:
            error_percent = (
                abs((x_new - x_current) / x_new) * 100
            )
        else:
            error_percent = 0.0 if f(x_new) == 0 else float("inf")

        converged = error_percent < tolerance

        rows.append({
            "k": k,
            "x_old": x_current,
            "f_x": fx,
            "df_x": dfx,
            "x_new": x_new,
            "error_percent": error_percent,
            "converged": converged
        })

        print(
            f"{k:3d} "
            f"{x_current:16.10f} "
            f"{fx:18.10f} "
            f"{dfx:18.10f} "
            f"{x_new:16.10f} "
            f"{error_percent:20.10f} "
            f"{'YES' if converged else 'No':>13}"
        )

        # ----------------------------------------------------
        # Check stopping criterion
        # ----------------------------------------------------
        if converged:
            root = x_new
            break

        # Move to next iteration
        x_current = x_new

    # --------------------------------------------------------
    # 6. FINAL RESULT
    # --------------------------------------------------------
    print("\n" + "=" * 75)
    print("NEWTON-RAPHSON RESULT")
    print("=" * 75)

    if root is not None:

        print(f"\nApproximate root = {root:.10f}")
        print(f"f(root)          = {f(root):.10e}")
        print(f"|f(root)|        = {abs(f(root)):.10e}")
        print(f"Approx. error    = {rows[-1]['error_percent']:.10e}%")
        print(f"Tolerance        = {tolerance:.10e}%")
        print(f"Iterations       = {len(rows)}")

    else:

        root = rows[-1]["x_new"]

        print("\nMaximum number of iterations reached.")
        print("The requested tolerance was not achieved.")

        print(f"\nLast approximation = {root:.10f}")
        print(f"f(x)               = {f(root):.10e}")
        print(f"Approx. error      = {rows[-1]['error_percent']:.10e}%")
        print(f"Iterations         = {len(rows)}")


    result = {
        "method": "Newton-Raphson", "function": f,
        "root": root if root is not None else x_new,
        "iterations": [initial_x] + [row["x_new"] for row in rows if "x_new" in row], "converged": rows[-1]["error_percent"] < tolerance
    }
    saved_results["newton"] = saved_results["root"] = result
    print("Results saved. Options 38-39 use this root-finding run.")


# ============================================================
# OPTION 33 — SECANT METHOD
# ============================================================

def secant_method():
    saved_results.pop("root", None)
    saved_results.pop("secant", None)
    print("\n" + "=" * 75)
    print("SECANT METHOD")
    print("=" * 75)

    # --------------------------------------------------------
    # 1. INPUT FUNCTION
    # --------------------------------------------------------
    f_input = input("\nEnter f(x): ")

    try:
        # Convert typed formula into a SymPy expression
        f_expr = convert_expression(f_input)

        # Convert symbolic expression into numerical function
        f = sp.lambdify(x_symbol, f_expr, "math")

    except Exception as e:
        print("\nInvalid function.")
        print("Error:", e)
        return

    # --------------------------------------------------------
    # 2. DISPLAY FUNCTION
    # --------------------------------------------------------
    print("\nFunction:")
    print(f"f(x) = {formula_string(f_expr)}")

    # --------------------------------------------------------
    # 3. INPUT TWO STARTING POINTS
    # --------------------------------------------------------
    x_previous = get_number(
        "\nEnter first starting point x^0: "
    )

    x_current = get_number(
        "Enter second starting point x^1: "
    )

    # --------------------------------------------------------
    # 4. INPUT STOPPING CRITERION
    # --------------------------------------------------------
    tolerance = get_number(
        "\nEnter tolerance epsilon [default = 1e-5]: ",
        1e-5
    )

    max_iterations = int(
        get_number(
            "Enter maximum number of iterations [default = 100]: ",
            100
        )
    )

    if tolerance <= 0 or max_iterations < 1:
        print("Tolerance and maximum iterations must be positive.")
        return

    # --------------------------------------------------------
    # 5. SECANT ITERATIONS
    # --------------------------------------------------------
    print("\n" + "=" * 110)
    print("SECANT METHOD ITERATIONS")
    print("=" * 110)

    print(
        f"{'k':>3} "
        f"{'x^k':>16} "
        f"{'f(x^k)':>18} "
        f"{'x^(k+1)':>16} "
        f"{'|f(x^k)|':>18} "
        f"{'Converged?':>13}"
    )

    print("-" * 110)

    initial_points = [x_previous, x_current]
    root = None
    rows = []

    # Evaluate the two initial points
    f_previous = f(x_previous)
    f_current = f(x_current)
    x_new, f_new = x_current, f_current

    # --------------------------------------------------------
    # Initial values
    # --------------------------------------------------------
    print(
        f"{0:3d} "
        f"{x_previous:16.10f} "
        f"{f_previous:18.10f} "
        f"{'-':>16} "
        f"{abs(f_previous):18.10e} "
        f"{'':>13}"
    )

    print(
        f"{1:3d} "
        f"{x_current:16.10f} "
        f"{f_current:18.10f} "
        f"{'-':>16} "
        f"{abs(f_current):18.10e} "
        f"{'':>13}"
    )

    # --------------------------------------------------------
    # Secant iteration
    # --------------------------------------------------------
    for k in range(1, max_iterations + 1):
        if abs(f_current) < tolerance:
            root = x_current
            break
        if abs(f_previous) < tolerance:
            root = x_previous
            break


        denominator = f_current - f_previous

        # ----------------------------------------------------
        # Prevent division by zero
        # ----------------------------------------------------
        if abs(denominator) < 1e-14:
            print("\nERROR:")
            print("The denominator became zero or extremely small.")
            print("The Secant Method cannot continue.")
            return

        # ----------------------------------------------------
        # Secant equation
        #
        # x^(k+1) =
        # x^k - f(x^k)(x^k-x^(k-1))
        #       -----------------------
        #       f(x^k)-f(x^(k-1))
        # ----------------------------------------------------
        x_new = (
            x_current
            - f_current * (x_current - x_previous)
            / denominator
        )

        f_new = f(x_new)

        # ----------------------------------------------------
        # Stopping criterion from the lecture:
        #
        # |f(x^k)| < epsilon
        # ----------------------------------------------------
        converged = abs(f_new) < tolerance

        rows.append({
            "k": k + 1,
            "x": x_new,
            "f_x": f_new,
            "converged": converged
        })

        print(
            f"{k + 1:3d} "
            f"{x_current:16.10f} "
            f"{f_current:18.10f} "
            f"{x_new:16.10f} "
            f"{abs(f_new):18.10e} "
            f"{'YES' if converged else 'No':>13}"
        )

        # ----------------------------------------------------
        # Check stopping criterion
        # ----------------------------------------------------
        if converged:
            root = x_new
            break

        # ----------------------------------------------------
        # Shift the points:
        #
        # x^(k-1) <- x^k
        # x^k     <- x^(k+1)
        # ----------------------------------------------------
        x_previous = x_current
        f_previous = f_current

        x_current = x_new
        f_current = f_new

    # --------------------------------------------------------
    # 6. FINAL RESULT
    # --------------------------------------------------------
    print("\n" + "=" * 75)
    print("SECANT METHOD RESULT")
    print("=" * 75)

    if root is not None:

        print(f"\nApproximate root = {root:.10f}")
        print(f"f(root)          = {f(root):.10e}")
        print(f"|f(root)|        = {abs(f(root)):.10e}")
        print(f"Tolerance        = {tolerance:.10e}")
        print(f"Iterations       = {len(rows)}")

    else:

        print("\nMaximum number of iterations reached.")
        print("The requested tolerance was not achieved.")

        print(f"\nLast approximation = {x_new:.10f}")
        print(f"f(x)               = {f_new:.10e}")
        print(f"|f(x)|             = {abs(f_new):.10e}")


    result = {
        "method": "Secant", "function": f,
        "root": root if root is not None else x_new,
        "iterations": initial_points + [row["x"] for row in rows], "converged": root is not None
    }
    saved_results["secant"] = saved_results["root"] = result
    print("Results saved. Options 38-39 use this root-finding run.")

# ============================================================
# GRAPHICAL REPRESENTATION — GENERAL PLOTTING FUNCTIONS
# ============================================================



# ------------------------------------------------------------
# 1. Plot one function
# ------------------------------------------------------------

def plot_function(f, x_start, x_end, title="Function"):
    """
    Plot a single mathematical function over a specified domain.
    """

    x_values = np.linspace(x_start, x_end, 500)
    y_values = [f(x) for x in x_values]

    plt.figure(figsize=(8, 5))

    plt.plot(x_values, y_values)

    plt.axhline(0, linewidth=0.8)
    plt.axvline(0, linewidth=0.8)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)

    plt.grid(True, alpha=0.3)
    plt.show()


# ------------------------------------------------------------
# 2. Compare two functions
# ------------------------------------------------------------

def plot_function_comparison(
    f1,
    f2,
    x_start,
    x_end,
    label1="Function 1",
    label2="Function 2",
    title="Function Comparison"
):
    """
    Plot two functions over the same domain.
    Useful for comparing:
        numerical vs exact
        Lagrangian vs exact
        Vandermonde vs exact
        fitted function vs exact
    """

    x_values = np.linspace(x_start, x_end, 500)

    y1_values = [f1(x) for x in x_values]
    y2_values = [f2(x) for x in x_values]

    plt.figure(figsize=(8, 5))

    plt.plot(
        x_values,
        y1_values,
        label=label1
    )

    plt.plot(
        x_values,
        y2_values,
        label=label2
    )

    plt.axhline(0, linewidth=0.8)
    plt.axvline(0, linewidth=0.8)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 3. Plot numerical data against an exact function
# ------------------------------------------------------------

def plot_numerical_vs_exact(
    x_values,
    y_values,
    exact_function,
    x_start=None,
    x_end=None,
    title="Numerical Solution vs Exact Solution",
    numerical_label="Numerical",
    exact_label="Exact"
):
    """
    Plot discrete numerical results together with
    the exact continuous solution.
    """

    x_values = np.array(x_values)
    y_values = np.array(y_values)

    if x_start is None:
        x_start = min(x_values)

    if x_end is None:
        x_end = max(x_values)

    x_exact = np.linspace(x_start, x_end, 500)
    y_exact = [exact_function(x) for x in x_exact]

    plt.figure(figsize=(8, 5))

    # Numerical points
    plt.plot(
        x_values,
        y_values,
        marker="o",
        label=numerical_label
    )

    # Exact function
    plt.plot(
        x_exact,
        y_exact,
        label=exact_label
    )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 4. Plot several numerical methods together
# ------------------------------------------------------------

def plot_method_comparison(
    results,
    exact_function=None,
    x_start=None,
    x_end=None,
    title="Comparison of Numerical Methods"
):
    """
    results should be a dictionary of the form:

    {
        "Euler": (x_values, y_values),
        "RK2.1": (x_values, y_values),
        "RK4": (x_values, y_values)
    }

    """

    plt.figure(figsize=(9, 6))

    for method_name, data in results.items():

        x_values, y_values = data

        plt.plot(
            x_values,
            y_values,
            marker="o",
            markersize=3,
            label=method_name
        )

    # Plot exact solution if available
    if exact_function is not None:

        if x_start is None:
            x_start = min(
                min(data[0]) for data in results.values()
            )

        if x_end is None:
            x_end = max(
                max(data[0]) for data in results.values()
            )

        x_exact = np.linspace(x_start, x_end, 500)
        y_exact = [exact_function(x) for x in x_exact]

        plt.plot(
            x_exact,
            y_exact,
            linewidth=2,
            label="Exact"
        )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 5. Plot error comparison
# ------------------------------------------------------------

def plot_error_data(
    error_results,
    title="Error Comparison"
):
    """
    error_results should be:

    {
        "Euler": (x_values, errors),
        "RK2.1": (x_values, errors),
        "RK4": (x_values, errors)
    }
    """

    plt.figure(figsize=(9, 6))

    for method_name, data in error_results.items():

        x_values, errors = data

        plt.plot(
            x_values,
            errors,
            marker="o",
            markersize=3,
            label=method_name
        )

    plt.xlabel("x")
    plt.ylabel("Absolute Error")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 6. Plot data points and approximation
# ------------------------------------------------------------

def plot_data_and_fit(
    x_data,
    y_data,
    approximation,
    x_start=None,
    x_end=None,
    title="Data and Approximation",
    data_label="Data",
    fit_label="Approximation"
):
    """
    Used for:
        Vandermonde
        Lagrangian
        Least Squares
    """

    x_data = np.array(x_data)
    y_data = np.array(y_data)

    if x_start is None:
        x_start = min(x_data)

    if x_end is None:
        x_end = max(x_data)

    x_fit = np.linspace(x_start, x_end, 500)
    y_fit = [approximation(x) for x in x_fit]

    plt.figure(figsize=(8, 5))

    # Original data points
    plt.scatter(
        x_data,
        y_data,
        label=data_label
    )

    # Approximation
    plt.plot(
        x_fit,
        y_fit,
        label=fit_label
    )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 7. Plot approximation vs exact function
# ------------------------------------------------------------

def plot_approximation_vs_exact(
    x_data,
    approximation,
    exact_function,
    x_start=None,
    x_end=None,
    title="Approximation vs Exact Function",
    approximation_label="Approximation",
    exact_label="Exact"
):
    """
    Used for:
        Vandermonde vs exact
        Lagrangian vs exact
        Least-squares vs exact
    """

    x_data = np.array(x_data)

    if x_start is None:
        x_start = min(x_data)

    if x_end is None:
        x_end = max(x_data)

    x_plot = np.linspace(x_start, x_end, 500)

    y_approx = [
        approximation(x)
        for x in x_plot
    ]

    y_exact = [
        exact_function(x)
        for x in x_plot
    ]

    plt.figure(figsize=(8, 5))

    plt.plot(
        x_plot,
        y_approx,
        label=approximation_label
    )

    plt.plot(
        x_plot,
        y_exact,
        label=exact_label
    )

    plt.scatter(
        x_data,
        [approximation(x) for x in x_data],
        s=20
    )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 8. Plot derivative approximation vs exact derivative
# ------------------------------------------------------------

def plot_derivative_comparison(
    x_values,
    derivative_values,
    exact_derivative,
    title="Derivative Approximation vs Exact",
    numerical_label="Numerical",
    exact_label="Exact"
):
    """
    Used for finite-difference and Lagrangian derivative
    comparisons.
    """

    x_values = np.array(x_values)
    derivative_values = np.array(derivative_values)

    x_exact = np.linspace(
        min(x_values),
        max(x_values),
        500
    )

    y_exact = [
        exact_derivative(x)
        for x in x_exact
    ]

    plt.figure(figsize=(8, 5))

    plt.plot(
        x_values,
        derivative_values,
        marker="o",
        label=numerical_label
    )

    plt.plot(
        x_exact,
        y_exact,
        label=exact_label
    )

    plt.xlabel("x")
    plt.ylabel("Derivative")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 9. Plot second derivative comparison
# ------------------------------------------------------------

def plot_second_derivative_comparison(
    x_values,
    derivative_values,
    exact_second_derivative,
    title="Second Derivative Approximation vs Exact",
    numerical_label="Numerical",
    exact_label="Exact"
):
    """
    Used for finite-difference and Lagrangian
    second-derivative comparisons.
    """

    x_values = np.array(x_values)
    derivative_values = np.array(derivative_values)

    x_exact = np.linspace(
        min(x_values),
        max(x_values),
        500
    )

    y_exact = [
        exact_second_derivative(x)
        for x in x_exact
    ]

    plt.figure(figsize=(8, 5))

    plt.plot(
        x_values,
        derivative_values,
        marker="o",
        label=numerical_label
    )

    plt.plot(
        x_exact,
        y_exact,
        label=exact_label
    )

    plt.xlabel("x")
    plt.ylabel("Second Derivative")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 10. Plot root-finding function and root
# ------------------------------------------------------------

def plot_root(
    f,
    root,
    x_start,
    x_end,
    title="Root-Finding"
):
    """
    Plot f(x) and mark the computed root.
    """

    x_values = np.linspace(
        x_start,
        x_end,
        500
    )

    y_values = [
        f(x)
        for x in x_values
    ]

    plt.figure(figsize=(8, 5))

    plt.plot(
        x_values,
        y_values,
        label="f(x)"
    )

    plt.axhline(
        0,
        linewidth=0.8
    )

    plt.scatter(
        [root],
        [f(root)],
        s=60,
        label=f"Estimate = {root:.6f}"
    )

    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ------------------------------------------------------------
# 11. Plot root-finding iterations
# ------------------------------------------------------------

def plot_root_iterations(
    iteration_values,
    true_root=None,
    title="Root-Finding Convergence"
):
    """
    iteration_values:
        list of successive root estimates

    Example:
        [1.5, 1.152057, 1.091851, 1.090278, ...]
    """

    iterations = range(
        len(iteration_values)
    )

    plt.figure(figsize=(8, 5))

    plt.plot(
        iterations,
        iteration_values,
        marker="o",
        label="Root estimate"
    )

    if true_root is not None:

        plt.axhline(
            true_root,
            linestyle="--",
            label="True root"
        )

    plt.xlabel("Iteration")
    plt.ylabel("Root estimate")
    plt.title(title)

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


# ============================================================
# MENU
# ============================================================


def plot_vandermonde():
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


def plot_root_finding():
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


def plot_root_convergence():
    result = saved_results.get("root")
    if result is None:
        print("Run option 31, 32, or 33 first to calculate root iterations.")
        return
    plot_root_iterations(result["iterations"], title=result["method"] + " iterations")


# ============================================================
# PROBLEM SET 2: INDEPENDENT WORKED-PROBLEM BLOCKS (OPTIONS 40-48)
# ============================================================
# These settings belong to the tutorial and do not change the main ODE.
# Keep full precision internally; round only when displaying results.
# Each calculation saves its results. Option 48 only plots saved results.

problem_set_results = {}


def ps2_1a():
    """Problem 1a: RK4 on two grids, h_m = h0 / 2**m."""
    # y' + y**2 + exp(-x) - x = 0, y(0) = 1.
    def f1(x, y):
        return x - y**2 - math.exp(-x)

    initial_h = get_number("Problem 1 h0 (given 0.4; ENTER uses 0.4): ", 0.4)
    if initial_h <= 0 or not math.isfinite(initial_h):
        print("h0 must be finite and positive.")
        return
    n0 = round(2 / initial_h)
    if n0 < 1 or not math.isclose(n0 * initial_h, 2, rel_tol=1e-9):
        print("h0 must divide [0, 2] into an integer number of intervals.")
        return
    if not math.isclose(initial_h, 0.4):
        print("NOTE: this uses a changed h0, not the printed value 0.4.")

    # A new grid invalidates the calculations based on the previous grid.
    for key in ["1a", "1b", "1c", "1d"]:
        problem_set_results.pop(key, None)

    runs = {}
    for m in [0, 1]:
        h = initial_h / (2 ** m)
        xs, ys, actual_h = rk4(f1, 0.0, 1.0, 2.0, h)
        runs[m] = {"x": xs, "y": ys, "h": actual_h}
        print(f"\nProblem 1a: m = {m}, h = {actual_h:.10f}")
        print(f"{'i':>3} {'x_i':>14} {'RK4 y_i':>20}")
        for i, (xi, yi) in enumerate(zip(xs, ys)):
            print(f"{i:3d} {xi:14.8f} {yi:20.12f}")

    problem_set_results["1a"] = {"runs": runs, "h0": initial_h}
    print("Results saved. Use option 48 to plot both grids on one graph.")
    return problem_set_results["1a"]


def ps2_1b():
    """Problem 1b: compare coarse/fine values at the same x."""
    data = problem_set_results.get("1a")
    if data is None:
        print("Run option 40 (Problem 1a) first.")
        return
    coarse, fine = data["runs"][0], data["runs"][1]
    requested_x = [0.0, 0.5, 1.0, 1.5, 2.0]

    print("\nProblem 1b has a grid inconsistency when h0 = 0.4.")
    print("The requested 0, 0.5, 1, 1.5, 2 are not all shared grid nodes.")
    print("1. Requested x-values, with explicitly labelled Lagrange interpolation")
    print("2. Actual common nodes of the saved grids (no interpolation)")
    mode = input("Choose 1 or 2 (ENTER uses 2): ").strip() or "2"
    if mode not in ["1", "2"]:
        print("Invalid choice.")
        return

    if mode == "1":
        targets = requested_x
        note = "Requested points; off-grid values use Lagrange interpolation."
        print("ASSUMPTION: interpolation is used at off-grid points.")
        print("Confirm this interpretation with the instructor; it is not specified in 1b.")
    else:
        targets = list(coarse["x"])
        note = "Actual shared grid nodes; this differs from the printed x-list if h0=0.4."
    print(note)

    # Values at existing nodes are used directly. At other requested points,
    # P(x) = sum_i y_i * product_(j != i) [(x-x_j)/(x_i-x_j)].
    values = []
    for grid in [coarse, fine]:
        interpolated_values = []
        for target in targets:
            node = next((i for i, xi in enumerate(grid["x"])
                         if math.isclose(xi, target, abs_tol=1e-10)), None)
            if node is not None:
                interpolated_values.append((grid["y"][node], "node"))
                continue
            total = 0.0
            for i in range(len(grid["x"])):
                phi = 1.0
                for j in range(len(grid["x"])):
                    if i != j:
                        phi *= (target-grid["x"][j]) / (grid["x"][i]-grid["x"][j])
                total += grid["y"][i] * phi
            interpolated_values.append((total, "interpolated"))
        values.append(interpolated_values)

    print("\nAbsolute approximate relative error: |(y_fine-y_coarse)/y_fine|*100")
    print("This compares the two approximations; it is not a true-error estimate.")
    print(f"{'x':>9} {'y (m=0)':>18} {'y (m=1)':>18} {'epsilon_a %':>16} {'source m0 / m1':>27}")
    rows = []
    for i, target in enumerate(targets):
        yc, source_c = values[0][i]
        yf, source_f = values[1][i]
        error = abs((yf-yc)/yf)*100 if yf != 0 else float("nan")
        rows.append({"x": target, "coarse": yc, "fine": yf, "error": error,
                     "coarse_source": source_c, "fine_source": source_f})
        print(f"{target:9.5f} {yc:18.12f} {yf:18.12f} {error:16.8f} {source_c+' / '+source_f:>27}")
    problem_set_results["1b"] = {"rows": rows, "note": note}
    return problem_set_results["1b"]


def ps2_1c():
    """Problem 1c: exact-fit Lagrangian polynomial of the fine RK4 data."""
    data = problem_set_results.get("1a")
    if data is None:
        print("Run option 40 (Problem 1a) first.")
        return
    problem_set_results.pop("1d", None)
    grid = data["runs"][1]
    x = x_symbol
    # Rational arithmetic avoids cancellation while expanding the polynomial.
    # It preserves the stored numerical data; it does not make RK4 exact.
    xs = [sp.Rational(str(xi)) for xi in grid["x"]]
    ys = [sp.Rational(str(yi)) for yi in grid["y"]]
    polynomial = 0
    for i in range(len(xs)):
        phi = 1
        for j in range(len(xs)):
            if i != j:
                phi *= (x-xs[j]) / (xs[i]-xs[j])
        polynomial += ys[i] * phi
    polynomial = sp.expand(polynomial)
    first = sp.diff(polynomial, x)
    second = sp.diff(polynomial, x, 2)

    print(f"\nProblem 1c: degree {len(xs)-1} exact fit of the m=1 RK4 data")
    print("Exact fit means passing through the data, not the exact IVP solution.")
    print("f(x) =", sp.N(polynomial, 16))
    print("f'(x) =", sp.N(first, 16))
    print("f''(x) =", sp.N(second, 16))
    residuals = [abs(float(polynomial.subs(x, xi)-yi)) for xi, yi in zip(xs, ys)]
    print(f"Maximum interpolation residual at data nodes = {max(residuals):.3e}")
    result = {"x": list(grid["x"]), "y": list(grid["y"]),
              "polynomial": polynomial, "first": first, "second": second,
              "function": sp.lambdify(x, polynomial, "math")}
    problem_set_results["1c"] = result
    print("Results saved. Use option 48 to plot f(x).")
    return result


def ps2_1d():
    """Problem 1d: secant solution of g(x)=f'(x)=0; classify the point."""
    data = problem_set_results.get("1c")
    if data is None:
        print("Run option 42 (Problem 1c) first.")
        return
    problem_set_results.pop("1d", None)
    # The assignment calls the root-finding function y(x); call it g here
    # to distinguish it from the IVP solution y and its polynomial fit f.
    g = sp.lambdify(x_symbol, data["first"], "math")
    second = sp.lambdify(x_symbol, data["second"], "math")
    print("\nRoot-finding problem: g(x) = f'(x) = 0")
    print("g(x) =", sp.N(data["first"], 16))
    previous = get_number("First secant point (default 0.4): ", 0.4)
    current = get_number("Second secant point (default 0.8): ", 0.8)
    tolerance = 1e-5
    history = [previous, current]
    converged = False
    print(f"{'k':>4} {'x_k':>20} {'g(x_k)':>20}")
    print(f"{0:4d} {previous:20.12f} {g(previous):20.12e}")
    print(f"{1:4d} {current:20.12f} {g(current):20.12e}")
    for k in range(2, 102):
        if abs(g(current)) < tolerance:
            converged = True
            break
        denominator = g(current)-g(previous)
        if abs(denominator) < 1e-14:
            print("Secant cannot continue: equal or nearly equal function values.")
            break
        next_x = current - g(current)*(current-previous)/denominator
        if not math.isfinite(next_x):
            print("Secant produced a nonfinite estimate; choose new starting points.")
            break
        previous, current = current, next_x
        history.append(current)
        print(f"{k:4d} {current:20.12f} {g(current):20.12e}")
    converged = abs(g(current)) < tolerance
    in_domain = 0 <= current <= 2
    curvature = second(current)
    minimum = converged and in_domain and curvature > 0
    print(f"\nStationary-point estimate x = {current:.12f}")
    print(f"f(x) = {data['function'](current):.12f}")
    print(f"|g(x)| = {abs(g(current)):.12e}; required < 1e-5")
    print(f"f''(x) = {curvature:.12f}")
    print("Verified local minimum in [0,2]." if minimum else
          "A local minimum in [0,2] was not verified; inspect the graph/start values.")
    result = {"x": current, "y": data["function"](current), "history": history,
              "residual": abs(g(current)), "curvature": curvature, "minimum": minimum}
    problem_set_results["1d"] = result
    return result


def ps2_2a():
    """Problem 2a: one-pass PC, and true errors at the provided exact nodes."""
    # y' + 2y - sqrt(x)*sin(2x) = 0, y(0)=0, h=0.15.
    def f2(x, y):
        return math.sqrt(x)*math.sin(2*x) - 2*y

    for key in ["2a", "2b", "2c"]:
        problem_set_results.pop(key, None)
    xs, ys, h = predictor_corrector(f2, 0.0, 0.0, 2.4, 0.15)
    print("\nProblem 2a: PC, h = 0.15, one corrector pass per step")
    print(f"{'i':>3} {'x_i':>14} {'PC y_i':>20}")
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        print(f"{i:3d} {xi:14.8f} {yi:20.12f}")

    exact_x = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4]
    exact_y = [0.0, 0.0323187, 0.140207, 0.270998, 0.342188,
               0.290355, 0.104387, -0.164493, -0.419418]
    errors = []
    print("\nTrue absolute percentage errors: |(y_exact-y_PC)/y_exact|*100")
    print(f"{'x':>9} {'PC y':>18} {'given exact y':>18} {'abs error':>18} {'error %':>16}")
    for j, (xi, reference) in enumerate(zip(exact_x, exact_y)):
        numerical = ys[2*j]  # 0.3 = 2h; these are actual PC nodes.
        absolute = abs(reference-numerical)
        percent = absolute/abs(reference)*100 if reference != 0 else float("nan")
        errors.append(percent)
        label = f"{percent:.8f}" if reference != 0 else "N/A (exact=0)"
        print(f"{xi:9.5f} {numerical:18.12f} {reference:18.12f} {absolute:18.12f} {label:>16}")
    result = {"x": xs, "y": ys, "h": h, "exact_x": exact_x,
              "exact_y": exact_y, "errors": errors}
    problem_set_results["2a"] = result
    print("Results saved. Use option 48 for the solution and error graphs.")
    return result


def ps2_2b():
    """Problem 2b: O(h**2) endpoint and central finite differences."""
    data = problem_set_results.get("2a")
    if data is None:
        print("Run option 44 (Problem 2a) first.")
        return
    xs, ys, h = data["x"], data["y"], data["h"]
    n = len(ys)-1
    first, second = [], []
    for i in range(n+1):
        if i == 0:
            d1 = (-3*ys[0] + 4*ys[1] - ys[2])/(2*h)
            d2 = (2*ys[0] - 5*ys[1] + 4*ys[2] - ys[3])/h**2
        elif i == n:
            d1 = (3*ys[n] - 4*ys[n-1] + ys[n-2])/(2*h)
            d2 = (2*ys[n] - 5*ys[n-1] + 4*ys[n-2] - ys[n-3])/h**2
        else:
            d1 = (ys[i+1]-ys[i-1])/(2*h)
            d2 = (ys[i-1]-2*ys[i]+ys[i+1])/h**2
        first.append(d1)
        second.append(d2)
    print("\nProblem 2b: central interior / higher-order endpoint differences")
    print(f"{'i':>3} {'x_i':>12} {'y_i':>18} {'first derivative':>20} {'second derivative':>20}")
    for i in range(n+1):
        print(f"{i:3d} {xs[i]:12.8f} {ys[i]:18.12f} {first[i]:20.12f} {second[i]:20.12f}")
    result = {"x": list(xs), "first": first, "second": second}
    problem_set_results["2b"] = result
    print("Results saved. Option 48 plots both derivatives on one graph.")
    return result


def ps2_2c():
    """Problem 2c: trapezoidal and Simpson integration of all PC nodes."""
    data = problem_set_results.get("2a")
    if data is None:
        print("Run option 44 (Problem 2a) first.")
        return
    xs, ys, h = data["x"], data["y"], data["h"]
    n = len(ys)-1
    trap_weights = [1] + [2]*(n-1) + [1]
    trap_sum = sum(trap_weights[i]*ys[i] for i in range(n+1))
    trapezoidal = h/2 * trap_sum

    # Simpson: 1,4,2,4,...,2,4,1; n=16 is even.
    simpson_weights = []
    for i in range(n+1):
        if i == 0 or i == n:
            simpson_weights.append(1)
        elif i % 2 == 1:
            simpson_weights.append(4)
        else:
            simpson_weights.append(2)
    simpson_sum = sum(simpson_weights[i]*ys[i] for i in range(n+1))
    simpson = h/3 * simpson_sum if n % 2 == 0 else None
    print(f"\nProblem 2c: all {n+1} PC nodes, n={n}, h={h:.8f}")
    print(f"{'i':>3} {'x':>10} {'y':>18} {'w_T':>5} {'w_T*y':>18} {'w_S':>5} {'w_S*y':>18}")
    for i in range(n+1):
        print(f"{i:3d} {xs[i]:10.5f} {ys[i]:18.12f} {trap_weights[i]:5d} "
              f"{trap_weights[i]*ys[i]:18.12f} {simpson_weights[i]:5d} {simpson_weights[i]*ys[i]:18.12f}")
    print(f"Trapezoidal: (h/2)*sum(w_i*y_i) = {trapezoidal:.12f}")
    if simpson is not None:
        print(f"Simpson 1/3: (h/3)*sum(w_i*y_i) = {simpson:.12f}")
    else:
        print("Simpson 1/3 is not applicable: n must be even.")
    result = {"trapezoidal": trapezoidal, "simpson": simpson,
              "trap_sum": trap_sum, "simpson_sum": simpson_sum, "n": n}
    problem_set_results["2c"] = result
    return result


def ps2_2d():
    """Problem 2d: degree-six least squares of GIVEN exact data; integrate f."""
    # These are the supplied exact values, not the PC approximations.
    xs = [sp.Rational(i*3, 10) for i in range(9)]
    ys = [sp.Rational(v) for v in ["0", "0.0323187", "0.140207", "0.270998",
           "0.342188", "0.290355", "0.104387", "-0.164493", "-0.419418"]]
    degree = 6
    # f(x) = a0 + a1*x + ... + a6*x**6.
    # S_ij = x_i**j; normal equations: (S.T*S)*a = S.T*y.
    # Rational arithmetic keeps these lesson equations without floating-point
    # round-off in the matrix solve. Input accuracy is still that of the table.
    S = sp.Matrix([[xi**j for j in range(degree+1)] for xi in xs])
    y_vector = sp.Matrix(ys)
    STS = S.T*S
    STy = S.T*y_vector
    coefficients = STS.LUsolve(STy)
    x = x_symbol
    polynomial = sum(coefficients[j]*x**j for j in range(degree+1))
    print("\nProblem 2d: degree-six least squares of the nine given exact values")
    print("S.T*S =")
    print(sp.N(STS, 10))
    print("S.T*y =")
    print(sp.N(STy, 10))
    print("\nCoefficients (constant term first):")
    for j, aj in enumerate(coefficients):
        print(f"a{j} = {float(aj):.16e}")
    print("f(x) =", sp.N(polynomial, 16))
    fitted = S*coefficients
    loss = sum((fitted[i]-ys[i])**2 for i in range(len(xs)))
    print(f"{'x':>10} {'given exact y':>18} {'fitted f(x)':>18} {'residual':>18}")
    for i in range(len(xs)):
        print(f"{float(xs[i]):10.5f} {float(ys[i]):18.12f} {float(fitted[i]):18.12f} "
              f"{float(fitted[i]-ys[i]):18.12f}")
    print(f"Loss sum((f(x_i)-y_i)**2) = {float(loss):.12e}")

    # Integrate the polynomial term by term, not by a quadrature approximation.
    # F(x) = sum_j a_j*x**(j+1)/(j+1); integral = F(b)-F(a).
    antiderivative = sum(coefficients[j]*x**(j+1)/(j+1) for j in range(degree+1))
    integral = antiderivative.subs(x, sp.Rational(24,10))-antiderivative.subs(x,0)
    print("F(x) =", sp.N(antiderivative, 16))
    print(f"Integral of fitted f(x), 0 to 2.4 = {float(integral):.12f}")
    result = {"x": list(map(float,xs)), "y": list(map(float,ys)),
              "coefficients": coefficients, "polynomial": polynomial,
              "function": sp.lambdify(x, polynomial, "math"),
              "antiderivative": antiderivative, "integral": float(integral), "loss": float(loss)}
    problem_set_results["2d"] = result
    print("Results saved. Option 48 can plot the fit and exact data.")
    return result


def ps2_plots():
    """Plot saved tutorial results only; never rerun any numerical method."""
    print("\nSaved-result plots: 1a, 1c, 1d, 2a, 2b, 2d, or all")
    choice = input("Plot which part? (ENTER uses all): ").strip().lower() or "all"
    parts = ["1a", "1c", "1d", "2a", "2b", "2d"] if choice == "all" else [choice]
    created_plot = False
    for part in parts:
        if part not in ["1a", "1c", "1d", "2a", "2b", "2d"]:
            print("Choose a listed part.")
            return
        data = problem_set_results.get(part)
        if data is None:
            print(f"No saved result for {part}; run its calculation first.")
            continue
        created_plot = True
        print(f"Preparing plot for {part}...", flush=True)
        if part == "1a":
            plt.figure()
            for m, grid in data["runs"].items():
                plt.plot(grid["x"], grid["y"], "o-", label=f"m={m}, h={grid['h']:g}")
            plt.title("Problem 1a: RK4 on two grids")
            plt.xlabel("x")
            plt.ylabel("y")
        elif part in ["1c", "2d"]:
            dense = np.linspace(min(data["x"]), max(data["x"]), 500)
            plt.figure()
            plt.plot(dense, [data["function"](x) for x in dense], label="f(x)")
            plt.scatter(data["x"], data["y"], label="RK4 data" if part=="1c" else "Given exact data")
            plt.title(f"Problem {part}: polynomial fit")
            plt.xlabel("x")
            plt.ylabel("y")
        elif part == "1d":
            fit = problem_set_results["1c"]
            dense = np.linspace(0, 2, 500)
            plt.figure()
            plt.plot(dense, [fit["function"](x) for x in dense], label="f(x)")
            plt.scatter([data["x"]], [data["y"]], label="Verified minimum" if data["minimum"] else "Stationary-point estimate")
            plt.title("Problem 1d: local minimum")
            plt.xlabel("x")
            plt.ylabel("f(x)")
        elif part == "2a":
            plt.figure()
            plt.plot(data["x"], data["y"], "o-", label="PC, h=0.15")
            plt.scatter(data["exact_x"], data["exact_y"], marker="x", label="Given exact values")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.title("Problem 2a: PC solution")
            plt.legend()
            plt.grid(True)
            plt.figure()
            valid = [(x,e) for x,e in zip(data["exact_x"],data["errors"]) if math.isfinite(e)]
            plt.plot([v[0] for v in valid], [v[1] for v in valid], "o-", label="Absolute true error %")
            plt.title("Problem 2a: true percentage errors (x=0 undefined)")
            plt.xlabel("x")
            plt.ylabel("Error (%)")
        else:  # Problem 2b: both requested derivatives on ONE graph.
            plt.figure()
            plt.plot(data["x"], data["first"], "o-", label="First derivative")
            plt.plot(data["x"], data["second"], "s--", label="Second derivative")
            plt.title("Problem 2b: finite-difference derivatives")
            plt.xlabel("x")
            plt.ylabel("Derivative")
        plt.legend()
        plt.grid(True)

    if created_plot:
        print("\nThe graphs open in separate Figure windows.", flush=True)
        print("If they are behind VS Code, use Alt+Tab to find them.", flush=True)
        print("Close ALL Figure windows to return to the menu.", flush=True)
        plt.show()
        print("Plot windows closed. Returning to the menu.", flush=True)




def menu():

    # Show problem when program starts
    print()

    print(
        "########## SHOWING PROBLEM ##########"
    )

    show_problem()

    input(
        "\nPress ENTER to continue to the menu..."
    )

    # --------------------------------------------------------
    # Main menu loop
    # --------------------------------------------------------

    while True:

        print()

        print("=" * 70)

        print(
            "                  NUMERICAL METHODS PROGRAM"
        )

        print("=" * 70)

        print()

        # ----------------------------------------------------
        # Numerical methods
        # ----------------------------------------------------

        print(
            "ODE NUMERICAL METHODS"
        )

        print(
            "    1. Euler"
        )

        print(
            "    2. Predictor-Corrector"
        )

        print(
            "    3. RK2.1"
        )

        print(
            "    4. RK2.2"
        )

        print(
            "    5. RK3"
        )

        print(
            "    6. RK4"
        )

        print()

        # ----------------------------------------------------
        # Analysis
        # ----------------------------------------------------

        print(
            "ANALYSIS"
        )

        print(
            "    7. Compare all methods"
        )

        print(
            "    8. Calculate final-point error"
        )

        print(
            "    9. Find minimum n"
        )

        print(
            "   10. Euler m-refinement"
        )

        print(
            "   11. RK4 m-refinement"
        )

        print(
            "   12. Show errors at all points"
        )

        print(
            "   13. Lagrange Interpolation - Estimate h_max"
        )

        print(
            "   14. RK4 Successive m - Approximate Error Criterion"
        )

        print()

        print(
            "NUMERICAL DIFFERENTIATION"
        )

        print(
            "   15. Finite Difference Approximations"
        )

        print(
            "   16. Finite Difference Approximations Tabulated"
        )

        print(
            "   17. Central Difference - Second Derivative"
        )

        print(
            "   18. Higher Order End Point Difference"
        )

        print()

        print(
            "NUMERICAL INTEGRATION"
        )

        print(
            "   19. Trapezoidal Rule"
        )

        print(
            "   20. Simpson's 1/3 Rule"
        )

        print(
            "   21. Compare Trapezoidal and Simpson's 1/3"
        )


        print()

        # ----------------------------------------------------
        # Plots
        # ----------------------------------------------------

        print(
            "PLOTS"
        )

        print(
            "   22. Plot numerical solutions"
        )

        print(
            "   23. Plot Euler maximum error vs m"
        )

        print(
            "   24. Plot final-point error vs n"
        )

        print(
            "   25. Plot error comparison (RK2.1 vs RK4)"
        )

        print()

        # ----------------------------------------------------
        # Problem information
        # ----------------------------------------------------

        print(
            "PROBLEM INFORMATION"
        )

        print(
            "   26. Show problem definition"
        )

        # ----------------------------------------------------
        # Function Approximation
        # ----------------------------------------------------
        
        print("\nFUNCTION APPROXIMATION")

        print("    27. Vandermonde Polynomial")


        print("    28. Lagrangian Polynomial")

        print("    29. Lagrangian and Finite Difference")

        print("    30. Least Squares")

   

        # ----------------------------------------------------
        # Root Finding
        # ----------------------------------------------------

        print("\nROOT-FINDING")

        print("    31. Bisection Method")

        print("    32. Newton-Raphson Method")

        print("    33. Secant Method")

        print()


        # ----------------------------------------------------
        # Function Approximation Plots
        # ----------------------------------------------------

        print("\nFUNCTION APPROXIMATION PLOTS")

        print("    34. Plot Vandermonde approximation")

        print("    35. Plot Lagrangian approximation")

        print("    36. Plot finite-difference derivatives")

        print("    37. Plot least-squares approximation")


        # ----------------------------------------------------
        # ROOT FINDING PLOTS
        # ----------------------------------------------------
        
        print("\nROOT-FINDING PLOTS")

        print("    38. Plot root and function")

        print("    39. Plot root-finding convergence")

        print("\nPROBLEM SET 2 (TUTORIALS)")
        print("   40. Problem 1a - RK4 on two grids")
        print("   41. Problem 1b - Pointwise approximate errors")
        print("   42. Problem 1c - Exact-fit polynomial")
        print("   43. Problem 1d - Local minimum by secant")
        print("   44. Problem 2a - PC and tabulated exact errors")
        print("   45. Problem 2b - First and second derivatives")
        print("   46. Problem 2c - Trapezoidal and Simpson integrals")
        print("   47. Problem 2d - Degree-6 least squares and integral")
        print("   48. Plot saved problem-set results")




        # ----------------------------------------------------
        # Exit
        # ----------------------------------------------------

        print(
            "    0. Exit"
        )

        print()

        choice = input(
            "Choose an option: "
        ).strip()

        # ====================================================
        # 1. EULER
        # ====================================================

        if choice == "1":

            h = get_number(
                f"Enter h "
                f"(press Enter for {h0}): ",
                h0
            )

            display_solution(
                "Euler",
                h
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 2. PREDICTOR-CORRECTOR
        # ====================================================

        elif choice == "2":

            h = get_number(
                f"Enter h "
                f"(press Enter for {h0}): ",
                h0
            )

            iterations = int(get_number("Corrector iterations per step (default 1): ", 1))
            display_solution("PC", h, corrector_iterations=iterations)

            if not proceed_to_menu():

                break

        # ====================================================
        # 3. RK2.1
        # ====================================================

        elif choice == "3":

            h = get_number(
                f"Enter h "
                f"(press Enter for {h0}): ",
                h0
            )

            display_solution(
                "RK2.1",
                h
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 4. RK2.2
        # ====================================================

        elif choice == "4":

            h = get_number(
                f"Enter h "
                f"(press Enter for {h0}): ",
                h0
            )

            display_solution(
                "RK2.2",
                h
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 5. RK3
        # ====================================================

        elif choice == "5":

            h = get_number(
                f"Enter h "
                f"(press Enter for {h0}): ",
                h0
            )

            display_solution(
                "RK3",
                h
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 6. RK4
        # ====================================================

        elif choice == "6":

            h = get_number(
                f"Enter h "
                f"(press Enter for {h0}): ",
                h0
            )

            display_solution(
                "RK4",
                h
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 7. COMPARE ALL METHODS
        # ====================================================

        elif choice == "7":

            h = get_number(
                f"Enter h "
                f"(press Enter for {h0}): ",
                h0
            )

            compare_all_methods(
                h
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 8. FINAL-POINT ERROR
        # ====================================================

        elif choice == "8":

            if exact_at_xf is None:

                print()

                print(
                    "No exact value at xf is available."
                )

                if not proceed_to_menu():

                    break

                continue

            print()

            print(
                "Available methods:"
            )

            for name in methods:

                print(
                    f"    {name}"
                )

            print()

            method_name = input(
                "Enter method name exactly: "
            ).strip()

            if method_name not in methods:

                print(
                    "Invalid method."
                )

                if not proceed_to_menu():

                    break

                continue

            print()

            print(
                "Enter n values separated by commas."
            )

            print(
                "Example: 10, 20, 40, 80"
            )

            print()

            user_input = input(
                "n values: "
            )

            try:

                n_values = [

                    int(
                        value.strip()
                    )

                    for value in
                    user_input.split(",")

                ]

                final_point_error_table(
                    method_name,
                    n_values
                )

            except Exception:

                print(
                    "Invalid n values."
                )

            if not proceed_to_menu():

                break

        # ====================================================
        # 9. MINIMUM n
        # ====================================================

        elif choice == "9":

            if exact_at_xf is None:

                print()

                print(
                    "No exact value at xf is available."
                )

                if not proceed_to_menu():

                    break

                continue

            print()

            print(
                "Available methods:"
            )

            for name in methods:

                print(
                    f"    {name}"
                )

            print()

            method_name = input(
                "Enter method name exactly: "
            ).strip()

            if method_name not in methods:

                print(
                    "Invalid method."
                )

                if not proceed_to_menu():

                    break

                continue

            print()

            tolerance = get_number(
                "Enter maximum final-point "
                "error % "
                "(default 0.001): ",
                0.001
            )

            find_minimum_n(
                method_name,
                tolerance=tolerance
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 10. EULER m-REFINEMENT
        # ====================================================

        elif choice == "10":

            tolerance = get_number(
                "Enter maximum error % "
                "(default 1.0): ",
                1.0
            )

            euler_m_refinement(
                tolerance=tolerance
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 11. RK4 m-REFINEMENT
        # ====================================================

        elif choice == "11":

            tolerance = get_number(
                "Enter maximum % change "
                "(default 0.5): ",
                0.5
            )

            rk4_m_refinement(
                tolerance=tolerance
            )

            if not proceed_to_menu():

                break

        # ====================================================
        # 12. ALL POINT ERRORS
        # ====================================================

        elif choice == "12":

            show_all_point_errors()

            if not proceed_to_menu():

                break

        # ====================================================
        # 22. PLOT SOLUTIONS
        # ====================================================

        elif choice == "22":

            plot_solutions()

            if not proceed_to_menu():

                break

        # ====================================================
        # 23. PLOT EULER ERROR VS m
        # ====================================================

        elif choice == "23":

            plot_euler_error_vs_m()

            if not proceed_to_menu():

                break

        # ====================================================
        # 24. PLOT FINAL ERROR VS n
        # ====================================================

        elif choice == "24":

            plot_final_error_vs_n()

            if not proceed_to_menu():

                break

        # ====================================================
        # 26. SHOW PROBLEM
        # ====================================================

        elif choice == "26":

            show_problem()

            if not proceed_to_menu():

                break

        # ============================================================
        # 27. VANDERMONDE POLYNOMIAL
        # ============================================================

        elif choice == "27":

            vandermonde_approximation()

            if not proceed_to_menu():
                break


        # ============================================================
        # 28. LAGRANGIAN
        # ============================================================

        elif choice == "28":

            lagrangian_polynomial()

            if not proceed_to_menu():
                break


        # ============================================================
        # 29. LAGRANGIAN + FINITE DIFFERENCE
        # ============================================================

        elif choice == "29":

            lagrangian_finite_difference()

            if not proceed_to_menu():
                break

        # ========================================================
        # 30. LEAST SQUARES
        # ========================================================

        elif choice == "30":

            least_squares()

            if not proceed_to_menu():
                break

        # ========================================================
        # 31. BISECTION METHOD
        # ========================================================          

        elif choice == "31":

            bisection_method()

            if not proceed_to_menu():
                            break

        # ========================================================
        # 32. NEWTON RAPHSON METHOD
        # ======================================================== 

        elif choice == "32":

            newton_raphson_method()


            if not proceed_to_menu():
                            break


        # ========================================================
        # 33. SECANT METHOD
        # ======================================================== 

        elif choice == "33":

            secant_method()

            if not proceed_to_menu():
                            break

        elif choice == "34":
            plot_vandermonde()
            if not proceed_to_menu():
                break

        elif choice == "35":
            plot_lagrangian()
            if not proceed_to_menu():
                break

        elif choice == "36":
            plot_finite_difference()
            if not proceed_to_menu():
                break

        elif choice == "37":
            plot_least_squares()
            if not proceed_to_menu():
                break

        elif choice == "38":
            plot_root_finding()
            if not proceed_to_menu():
                break

        elif choice == "39":
            plot_root_convergence()
            if not proceed_to_menu():
                break



        # ====================================================
        # 13. LANGRAGIAN
        # ====================================================

        elif choice == "13":

            method_name = input("Enter method to use (Euler/RK2.1/RK2.2/RK3/RK4): ")

            if method_name in methods:

                lagrange_hmax(
                    methods[method_name],
                    f,
                    exact,
                    x0,
                    y0,
                    xf
                )
         
            else:
                print("Invalid method.")

            if not proceed_to_menu():
            
                break

        # ==========================================================
        # 14. RK4 SUCCESSIVE m - APPROXIMATE ERROR CRITERION
        # ==========================================================

        elif choice == "14":

            rk4_successive_m_error(
                f,
                x0,
                y0,
                xf,
                h0
            )

            if not proceed_to_menu():
                break


        # ==========================================================
        # 25. ERROR COMPARISON
        # ==========================================================

        elif choice == "25":

            plot_error_comparison(
                f,
                exact,
                x0,
                y0,
                xf,
                h0,
                methods
            )

            if not proceed_to_menu():
                break

        # ==========================================================
        # 15. FINITE DIFFERENCE APPROXIMATIONS
        # ==========================================================

        elif choice == "15":

            finite_difference_approximations()

            if not proceed_to_menu():
                break               

        # ==========================================================
        # 16. FINITE DIFFERENCE FROM TABULATED DATA
        # ==========================================================

        elif choice == "16":

            finite_difference_tabulated()

            if not proceed_to_menu():
                break


        # ==========================================================
        # 17. CENTRAL SECOND DERIVATIVE
        # ==========================================================

        elif choice == "17":

            central_second_derivative_tabulated()

            if not proceed_to_menu():
                break


        # ==========================================================
        # 18. HIGHER-ORDER FORWARD/BACKWARD DIFFERENCES
        # ==========================================================

        elif choice == "18":

            higher_order_endpoint_differences()

            if not proceed_to_menu():
                break

        elif choice == "19":

            trapezoidal_rule()

            if not proceed_to_menu():
                break


        elif choice == "20":

            simpson_one_third_rule()

            if not proceed_to_menu():
                break


        elif choice == "21":

            compare_integration_methods()

            if not proceed_to_menu():
                break



        # ====================================================
        # 0. EXIT
        # ====================================================

        # 40. Problem 1a - RK4 on two grids
        elif choice == "40":
            ps2_1a()
            if not proceed_to_menu():
                break

        # 41. Problem 1b - Pointwise approximate errors
        elif choice == "41":
            ps2_1b()
            if not proceed_to_menu():
                break

        # 42. Problem 1c - Exact-fit polynomial
        elif choice == "42":
            ps2_1c()
            if not proceed_to_menu():
                break

        # 43. Problem 1d - Local minimum by secant
        elif choice == "43":
            ps2_1d()
            if not proceed_to_menu():
                break

        # 44. Problem 2a - PC and tabulated exact errors
        elif choice == "44":
            ps2_2a()
            if not proceed_to_menu():
                break

        # 45. Problem 2b - First and second derivatives
        elif choice == "45":
            ps2_2b()
            if not proceed_to_menu():
                break

        # 46. Problem 2c - Trapezoidal and Simpson integrals
        elif choice == "46":
            ps2_2c()
            if not proceed_to_menu():
                break

        # 47. Problem 2d - Degree-6 least squares and integral
        elif choice == "47":
            ps2_2d()
            if not proceed_to_menu():
                break

        # 48. Plot saved problem-set results
        elif choice == "48":
            ps2_plots()
            if not proceed_to_menu():
                break

        elif choice == "0":

            print()

            print(
                "Exiting program."
            )

            break

        # ====================================================
        # INVALID
        # ====================================================

        else:

            print()

            print(
                "Invalid option. "
                "Please choose again."
            )


# ============================================================
# PROGRAM START
# ============================================================


if __name__ == "__main__":

    menu()