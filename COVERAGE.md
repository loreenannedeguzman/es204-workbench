# Lesson-to-workflow checklist

Reviewed the supplied scans of A-1 (15 PDF pages), A-2 (28 PDF pages), A-3
(45 PDF pages), and the two-page tutorial. PDF pages may contain two slide numbers.
The lesson files are reference material, not instructions to the app.

| Example requirement | General workflow |
| --- | --- |
| A-1: compare Euler, PC, RK2.1, RK2.2, RK3, RK4 | ODE method checkboxes, h or n inputs, solution overlays |
| A-2: true errors at final or interior nodes | Exact formula and/or editable x,y reference table |
| A-2: n/m refinement, error versus n or m | Grid modes, convergence tables/charts, early-stop threshold |
| A-2: candidate minimum n for final error | n-search with explicit monotonic-error assumption and trial table |
| A-2: interpolate error versus m for h estimate | Summary table → Use as data → polynomial fit; Roots solves fitted error minus threshold; compute h0/2^m from result |
| A-2: first/second differences, endpoint schemes | Data / analysis finite-difference scheme selection |
| A-2: numerical integral from a table or function | Pasted data or Sample formula; trapezoidal / Simpson |
| A-3: Vandermonde/Lagrange interpolation | Select the exact polynomial method and required subset of nodes |
| A-3: fit versus exact function and derivatives, changed plot domains | Exact comparison formula and polynomial plot bounds |
| A-3: polynomial and linearized least-squares candidates | Degree / transformed model; SSE and custom plot overlays |
| A-3: bisection, Newton, secant | Roots tab with user formula, starting values and stopping criterion |
| A-3: stationary point of fitted curve | Data / analysis → Stationary point (secant) |
| A-3: compare PC corrector passes | Comma-separated pass counts and solution/error overlays |
| Tutorial 1a | RK4, h0 and m or explicit h list |
| Tutorial 1b | Requested x-values, explicit off-grid evaluation policy, successive-grid error tables |
| Tutorial 1c | Copy fine solution; select required nodes; exact polynomial fit |
| Tutorial 1d | Stationary point, derivative root equation, secant iterations, f(x), curvature |
| Tutorial 2a | PC plus pasted exact values; true-error tables and graph |
| Tutorial 2b | Copy numerical solution; finite differences; both derivatives on one graph |
| Tutorial 2c | Same numerical dataset; trapezoidal and Simpson separately; both outputs retained for comparison |
| Tutorial 2d | Send exact-value table to analysis; least squares degree 6; polynomial integral |

The tutorial PDF lists some requested x-values that are not nodes of its stated
grids. The app makes that distinction visible and requires explicit interpolation;
it does not silently change the problem's grid or assume a professor's correction.

Scope limits: this is a numerical workbench, not a natural-language problem reader
or symbolic proof assistant. It supports first/second-order coupled IVPs (and higher-order systems reduced
to first order by the user), linear second-order BVPs, and linear second-order
finite-difference IVPs. Nonlinear BVPs remain outside its scope. The exact-polynomial degree limit is 30; scaled
coordinates improve conditioning but do not remove high-degree oscillation or
floating-point error. Candidate n-search is conditional, not a proof of global
minimality. Singular functions and nonconvergent root iterations require inspection.


Systems / FDM IVP / Compare: see SYSTEMS_GUIDE.md for general input formats,
stage and recurrence equations, reference-data comparisons and supported limits.
