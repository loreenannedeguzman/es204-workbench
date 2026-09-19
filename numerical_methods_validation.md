# ES204 numerical methods validation

Reviewed September 7, 2026. The current program contains options **1-39**, including the function-approximation and root-finding additions. Those additions were preserved.

**Conclusion:** the original file was not ready to use as a reliable study program: it could not start because of a syntax error, and two calculations did not implement the lesson methods. The revised file fixes those issues, preserves separate formula blocks, and passes the focused numerical and menu checks described below. This is a validation of the stated formulas and tested examples, not a guarantee for every function, input, or polynomial degree.

## Sources and scope

All pages of these three scanned PDFs were visually reviewed. Page references below mean PDF pages; slide numbers are included where helpful.

- [Lesson A-1: first-order IVPs](F:/ES204/notes/ES204-Unit%20I.A-1_2026notes.pdf): 15 PDF pages, including Euler, Predictor-Corrector, RK2.1, RK2.2, RK3, and RK4.
- [Lesson A-2: post-processing](F:/ES204/notes/ES204Unit%20II-A-2_2026notes.pdf): 28 PDF pages, covering errors, refinement, differentiation, and integration.
- [Lesson A-3: approximation and roots](F:/ES204/notes/ES204Unit%20II-A-3_2026notes.pdf): 45 PDF pages, covering polynomial approximation, least squares, root-finding, plots, and repeated corrector iteration.

Assignments and administrative directions inside the PDFs were treated as source material, not instructions to submit work or contact anyone.

## Coverage and accuracy

| Lesson topic | Menu options | Validation outcome |
|---|---|---|
| Euler, PC, RK2.1, RK3, RK4 | 1-3, 5-6 | Update formulas match the lessons. Worked examples and expected convergence orders pass. |
| RK2.2 | 4 | Corrected from midpoint to the lesson's Ralston version. |
| Error tables and refinement | 7-12, 14 | Percentage and successive-grid calculations checked; zero-reference handling corrected. Minimum-n search retains the qualification below. |
| Quadratic estimate of h_max | 13 | Corrected to h_m = h0/2^m, maximum percentage error, and solving p(m) = 1%, using m = 6, 7, 8. The equal-grid candidate is checked directly. |
| First/second differences and higher-order endpoints | 15-18 | Stencil coefficients and endpoint indexing match the lessons; verified with tabulated data. |
| Trapezoidal and Simpson 1/3 integration | 19-21 | Composite weights verified. Simpson enforces even n; option 19 now also handles unequal intervals as taught in A-2. |
| ODE plots and problem information | 22-26 | Existing behavior retained; repaired the overwritten error-plot function. |
| Vandermonde and Lagrangian polynomials | 27-28 | Matrix powers and basis products are correct. Interpolation at the nodes and agreement between the two methods checked. |
| Lagrangian versus finite differences | 29 | First- and second-derivative formulas match the lessons, including higher-order endpoint differences. |
| Least squares | 30 | Polynomial normal equations, log-linearized exponential fit, and comparison using loss in the original y-values verified. |
| Bisection, Newton-Raphson, secant | 31-33 | Update formulas and lesson stopping criteria verified. Repaired runtime failures and several zero/degenerate cases. |
| Approximation and root plots | 34-39 | Now consume results saved by options 27-33. They do not refit a model, rerun an iteration, or request the data again. |
| Repeated PC correction, A-3 Example 10 | 2 | Optional corrector count added; default one pass preserves the original behavior. The previous node's y remains fixed within the corrector loop. |

### Important formula corrections

**RK2.2, A-1 PDF page 10, slide 18.** With the code's convention that k values are slopes:

```python
k1 = f(x, y)
k2 = f(x + 3*h/4, y + 3*h*k1/4)
y = y + h*(k1/3 + 2*k2/3)
```

The notes instead include h in each k. These conventions are mathematically equivalent when used consistently. The original midpoint implementation was a valid second-order method, but it was not the RK2.2 specified in these notes.

**Lagrange h_max, A-2 PDF pages 10-11, slides 19-20.** The old block treated m as a number of subintervals, interpolated absolute error, and used the quadratic vertex. The lesson requires refinement level m, percentage error, and the decreasing threshold crossing. The corrected calculation estimates h_max near **0.005648** for Example 2. It reports and checks an integer-n grid rather than asserting that interpolation alone guarantees the tolerance.

### Lesson arithmetic and rounded references

- **A-2 PDF page 16, slide 31:** the displayed forward first difference is arithmetically inconsistent. From the printed values, `(-0.56966 - (-1))/0.5 = 0.86068`, not `1.15076`. Using the six-decimal data elsewhere in the lessons gives `0.860678`. The program correctly uses the formula; the associated percentage error is about **13.932%**, not the later table's 15.076%.
- The default exact function uses the lesson's rounded `2.16667`. It gives `y(0) = 6.50001`, not exactly 6.5. It remains unchanged for reproducing the printed reference tables. For mathematical error studies of that IVP, use `exact_input = "exp(x)*(6.5 - 3*x)"`; this alternative is now noted beside the setting.
- Small differences from printed tables can result from rounded input values and intermediate arithmetic. The validation compares both mathematical formulas and the precision shown in the examples.

## Runtime repairs and study structure

Fixed the malformed Newton header, calls to nonexistent `format_formula` and `proceeded_to_menu`, the duplicate `plot_error_comparison` name, and undefined plotting menu functions. Removed the repeated top-level plotting import. Updated stale option numbers inside docstrings.

The main numerical methods remain separate. PC and RK2.1 intentionally retain their own blocks, because their equivalence with one correction is part of the lesson. Finite-difference and integration routines also retain their visible formula calculations. No class hierarchy or general solver framework was introduced.

Mesh checks prevent zero division and silent changes to the requested h. Use **h = (xf - x0)/n** with a positive integer n and a forward interval. Relative percentage error is undefined at a zero reference: integration displays `nan` for that percentage while preserving the absolute error; ODE comparison routines treat an exactly matched zero as zero for comparison and a mismatched zero as infinity so it cannot falsely pass a tolerance test. This zero-match convention is not a mathematical definition of relative error.

## Saved-result plotting

| Calculate first | Plot afterward |
|---|---|
| 27: Vandermonde | 34 |
| 28: Lagrangian polynomial | 35 |
| 29: Lagrangian/finite differences | 36 |
| 30: Least squares | 37 |
| 31, 32, or 33: root-finding | 38: root/function; 39: iteration history |

Results are held in `saved_results` until the program exits. Each new calculation replaces that method's saved result; starting an unsuccessful calculation clears its prior result rather than plotting stale data. Root plots use the most recent root-finding run and distinguish convergence from a last approximation. They have no dependency on the current ODE problem definition.

Plotting still evaluates the **already calculated** polynomial or saved function on a plotting grid. That is necessary to draw a curve; it does not repeat the numerical solve. Existing ODE plots 22-25 retain their prior calculation behavior; the requested cache change targets 27-39.

## Verification

The accompanying `validate_numerical_methods.py` runs with the same Python dependencies as the program, using Matplotlib's noninteractive Agg backend. Checks include:

- A-1 Euler decay/instability and PC examples; Example 4 gives Euler 3.7024408821, PC/RK2.1 3.7535337830, and RK4 3.7515324640.
- A distinct nonlinear, nonautonomous Ralston step and convergence-order checks for all six solvers.
- A-2 Example 2 Euler refinement m = 7, RK4 refinement m = 1, and quadratic h_max interpolation.
- Finite-difference endpoint/interior values, trapezoidal/Simpson weights, unequal-spacing integration, and the integration exercise with n = 7 and n = 6.
- Interpolation consistency; A-3 Example 4's exponential model has the lowest original-scale loss and parameters approximately a = 1.999, b = 0.250.
- A-3 Example 6 bisection and Example 8 Newton (root approximately 1.090277), secant residual stopping, and PC corrector counts 1-4.
- Invalid mesh and zero-reference cases; menu startup/exit; one unique handler for every option 0-39; defined menu callables.
- Plotting with saved results, no calculation prompts, and linear solving disabled; missing-result messages; an actual calculate-then-plot menu sequence.

**Result: 259 checks passed**, including repeated static menu-callable checks. This count is not 259 independent mathematical examples. Charts were rendered through Agg for execution checks; every chart was not individually reviewed for visual styling.

## Remaining limits and study guidance

1. **Minimum n:** the fast doubling/binary search assumes error decreases with n. It checks n and n-1, but that does not prove a global minimum for nonmonotonic error. The program now states that assumption. Its maximum-n cutoff and finite refinement limits are search limits, not mathematical impossibility proofs.
2. **High-degree polynomials:** Vandermonde systems, expanded Lagrange polynomials, and normal equations can lose substantial accuracy. The lesson itself discusses this. Their simple forms were retained for learning; low-degree checks passing does not certify the 12th- or 24th-degree cases. Inspect residuals and compare against known values before trusting high-degree results.
3. **Lagrange h_max:** samples m = 6, 7, 8 belong to the lesson example. For another problem whose threshold is outside that window, use refinement first. The block now reports that case instead of extrapolating a misleading maximum.
4. **Input assumptions:** these are educational routines, not a fully defensive numerical library. Arbitrary discontinuities, overflow, invalid symbolic domains, or impossible tolerances are not exhaustively handled. Bisection's sign-change test presumes continuity. Newton and secant need suitable starting points.
5. **Manual study steps:** the menu does not classify ODEs, derive analytic IVP solutions, automatically perform all lesson exercises, or provide a dedicated exact integral of an interpolating polynomial. Approximation plots 34-37 use saved data/models; they do not automatically request an analytic overlay. Existing generic plotting helpers can support those overlays in code. Mentioned but untaught alternatives such as splines, Gaussian quadrature, and fixed-point iteration were not added.
6. **PC iteration exercise:** option 2 lets you repeat the calculation with 1-4 corrector passes. Side-by-side tables, analytic reference values, and integrating those results remain explicit study steps rather than an automated exercise runner.

Use the program to follow each recursion or stencil, then compare selected rows with a hand calculation. Avoid copying a rounded displayed polynomial back into another method when high accuracy is required; the saved plotting results keep full precision.
