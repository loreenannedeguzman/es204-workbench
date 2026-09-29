# Systems, initial-value finite differences, and comparisons

These tools accept user formulas and data. There is no assignment selector or
built-in homework solution. The existing ODE, analysis, roots, BVP and matrix
workflows retain their inputs and behavior.

## Systems

Choose **First order** for `u'=F(x,u,v,...)`, or **Second order** for
`u''=F(x,u,du,v,dv,...)`. Enter names such as `u,v`, one right-hand-side formula
per named variable, and initial values in the displayed order.

For second-order equations, the workbench introduces `du=u'`, `dv=v'`, etc.
Initial values must follow `u,du,v,dv,...`. For first-order equations, use
`u,v,...`. Names must not collide with x, mathematical constants/functions,
or generated derivative names. All components are evaluated simultaneously
at each stage, including coupled equations. For mixed/higher-order systems,
reduce them to first order yourself and enter the resulting state equations.

Choose Euler, predictor–corrector, either existing RK2 variant, RK3 or RK4.
PC uses an Euler predictor followed by the specified number of trapezoidal
corrector passes; it does not silently iterate to convergence. The grid is
explicit (h or intervals); h must divide the interval.

Outputs include:

- The reduced differential equations, initial conditions, grid, problem-specific
  stage formulas and update equation, in Report / formulas and LaTeX export.
- A solution table containing every state, and a stage-evaluation table containing
  the stage x, trial states and slopes. Slopes are derivatives, not h times slopes.
- Solution overlays; separate first-derivative overlays for second-order input.
- An optional versus plot: enter the horizontal state and vertical state explicitly.
  For `u versus v`, horizontal is `v` and vertical is `u`. Node order is retained.
- Every state as a dataset for subsequent analysis. Arbitrary additional overlays
  are available through Build custom versus plot.

Limits: 1–6 named variables (up to 12 reduced states), 2–10000 intervals,
1–100 PC passes, and at most 500000 scalar stage slopes per run. These are
full-table teaching tools, not adaptive or stiff-system solvers.

## FDM IVP

Enter `a(x)y'' + b(x)y' + c(x)y = d(x)`, both initial conditions at the left
endpoint, and the grid. The coefficient of y'' must be nonzero. This is distinct
from BVP finite differences, which impose conditions at opposite endpoints.

The central recurrence is

```text
(a_i + h*b_i/2)*y_(i+1)
  = (2*a_i-h^2*c_i)*y_i - (a_i-h*b_i/2)*y_(i-1) + h^2*d_i.
```

Choose the starting formula explicitly:

1. **Central ghost / Taylor (second order)** eliminates a ghost point using the
   centered initial derivative and the differential equation at the initial node.
   This gives `y1=y0+h*v0+h^2*(d0-b0*v0-c0*y0)/(2*a0)`.
2. **Forward derivative (first order)** uses `y1=y0+h*v0`. This generally reduces
   global accuracy to first order.
3. **Forward derivative (second order)** combines
   `(-3*y0+4*y1-y2)/(2*h)=v0` with the first interior recurrence to determine y1.

Use the scheme required by your lesson. The app does not infer that requirement.
Second-order stencils alone do not guarantee a small error: stability, starting
errors, initial-data precision and propagation through the ODE also matter.
Zero/near-zero recurrence denominators are rejected.

The report includes the start, every substituted recurrence and a numerical
recurrence table. The initial derivative is the supplied value; later derivative
columns use second-order finite-difference estimates. No right boundary value
is imposed. Optional reference pairs are matched at grid nodes only.

## Compare

Paste a table with a header. The first column is x, the second is the supplied
reference, and the remaining columns are numerical series:

```text
x, reference, method A, method B
0, 1, 1, 1
0.5, 1.64872, 1.645, 1.649
1, 2.71828, 2.7, 2.719
```

Alternatively, begin with only x and reference columns, then click **Add column
from result history**. Choose the numerical solution table, its x and value
columns, and a distinct series label. Repeat for other methods, then Compare.
Only matching nodes are accepted. The reference strings are retained during
column insertion, and no interpolation or replacement is performed.

Outputs: numerical/reference overlay, signed-error and percentage-error plots,
per-node signed and absolute errors, maximum absolute error, RMSE and maximum
absolute percentage error. Percentage error uses `(reference-numerical)/reference`
times 100. Zero references yield N/A percentages and are excluded from percentage
plots/statistics; their absolute errors are retained. The existing Signed/Absolute
plot control applies to percentage-error charts. Rounded references limit the
meaningful precision of comparisons. Plots use distinct markers and line styles.

## BVP additions

The BVP form now accepts a supplied reference table alongside the optional
analytical formula; the comparisons are separate and neither reference is
silently replaced. BVP results also include a combined y/y' chart.

Shooting already supports PC and the strict boundary-residual stopping test.
A small boundary residual does not establish mesh accuracy. If no step size is
specified in a problem, select and document one or confirm the expected grid.

## Saved work and updates

Save session preserves all new inputs, results and history. Older sessions keep
their original tab indices and load with defaults for new fields. An already-open
window is not closed by installation; save its session, reopen the application
when ready, then load that session to access the new tabs.

## Verification

Tests check analytic coupled solutions, RK4 convergence, PC passes, finite-
difference convergence and independent matrix solutions, zero/off-grid references,
comparison statistics, worker execution, plots, equation export and session
compatibility. The runtime has no SciPy dependency; independent validation can
use SciPy without adding it to the application.
