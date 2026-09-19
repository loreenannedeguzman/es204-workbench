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
    for part in parts:
        if part not in ["1a", "1c", "1d", "2a", "2b", "2d"]:
            print("Choose a listed part.")
            return
        data = problem_set_results.get(part)
        if data is None:
            print(f"No saved result for {part}; run its calculation first.")
            continue
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
            plt.show()
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
        plt.show()
