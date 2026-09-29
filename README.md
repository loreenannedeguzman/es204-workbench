# ES204 Numerical Workbench

Open **ES204 Numerical Methods** from the desktop. This version accepts your own
inputs and contains no built-in homework answers. Numerical routines are in
`numerics.py`; the Tkinter interface is in `ES204_App.pyw`.

## ODE workflow

1. Enter the right side of `y' = f(x,y)`, the initial x and y, and final x.
2. Select a grid mode:
   - **Step sizes h:** enter one or more comma-separated step sizes.
   - **Interval counts n:** enter one or more positive integer interval counts.
   - **Refinement h0 / 2^m:** enter a single initial h0 and the last m.
   - **Search n (final error):** enter a starting n, maximum n, a final reference,
     and a percentage tolerance. This brackets and searches under a monotone-error
     assumption; it reports a candidate minimum, not a guaranteed minimum for
     nonmonotonic errors.
3. Check one or more methods. PC can use a single corrector-pass count or a list
   to compare different counts. Each pass keeps the old y fixed.
4. Optionally enter an exact solution formula, an exact-value table, or both.
   Tables contain one `x, y` pair per line with no header. A final reference is a
   one-row table. Never put a final reference into the exact-formula field.
5. Optionally enter requested x-values. By default, missing grid nodes are N/A.
   Choose linear or polynomial interpolation explicitly if that is appropriate
   for the problem. Interpolated values are numerical estimates, not exact data.
6. Click **Solve ODE**. Select a table in the results dropdown or open **Plots**.

All ODE methods require a positive h dividing the forward interval. Input numbers
may be expressions such as `pi/4`. Formulas use explicit multiplication, for example
`2*x`, and support `^` or `**` for powers. This version solves scalar first-order
IVPs, not ODE systems or boundary-value problems.

True error is measured against your supplied references. Percentage error is N/A
when the exact value is zero. Successive-grid error is `abs((fine-coarse)/fine)*100`,
not a Richardson error estimate. Maximum errors exclude undefined percentages.
With a reference table, the maximum is over matching/evaluated reference points,
not an assertion about the whole continuous domain.

## Post-processing workflow

Open **Data / analysis**. Select an ODE result and click **Use selected solution as
data**, or paste/import points. The optional node-selection field lets you copy
only specific x nodes. You can also send the exact-value table directly to analysis.
Always select the data source required by the problem: numerical results and exact
reference values are different datasets.

Available operations:

- Exact Lagrange or Vandermonde interpolation; least-squares polynomial of a chosen
  degree. Coefficients, fit residuals, SSE, evaluations, derivatives, integral,
  and comparison plots are produced. Calculations use scaled coordinates for
  numerical stability. Floating-point interpolation is not symbolic exact arithmetic.
- Linearized least squares for exponential, power-law, reciprocal, and inverse-x
  models. SSE in the original y coordinates allows comparisons between models.
- Forward, backward, central, or higher-order finite differences. Unsupported
  endpoints are N/A rather than silently switching schemes. The combined default
  uses central interiors and higher-order endpoints. An exact function enables
  derivative-error tables and plots.
- Trapezoidal integration; Simpson's 1/3 integration on uniform grids with an even
  number of intervals. These integrate over the first and last data points. To
  integrate a subinterval, supply its points. A supplied exact integral enables
  integral-error calculation.
- Secant solution of the fitted polynomial's derivative for a stationary point,
  with residual tolerance, iteration table and curvature classification.

**Sample formula** can generate data from a formula over chosen bounds and n,
so differentiation and integration do not require an ODE. Polynomial plot/integral
bounds can differ from the data domain; extrapolation is explicitly reported.

For a given formula `g(x)=0`, use **Roots**: bisection, Newton-Raphson or secant,
with residual or relative-step stopping criteria. Inspect both convergence status
and the residual; a small step alone does not establish an accurate root.

## Versus plots and exports

- Built-in result charts include solution overlays, reference errors, convergence
  against n/m, derivative overlays, polynomial comparisons and root convergence.
- **Plots → Build custom versus plot** selects numeric X and Y columns from any
  table generated during this session. Add multiple series to compare results
  from different runs. Enter labels and choose linear/log axes. Non-numeric and
  N/A rows are omitted; log axes require positive data.
- **Tables → Use as data** copies any two numeric result columns into analysis.
  For example, fit error versus refinement level, then solve a threshold equation.
- Use the embedded plot toolbar to zoom and save PNG/PDF/SVG figures.
- Export tables as CSV or copy them into a spreadsheet. Save report writes all
  current tables and notes to a text file.
- Save/load inputs preserves your form setup in JSON. It does not store computed
  history. Loading a setup clears the old results; rerun to compute them.

Scroll the input panel to reach all fields. Inputs take effect only when you run.
Cancel terminates the calculation worker and keeps the previous successful results.
The table viewer shows up to 5,000 rows; CSV export retains every row.

## AI231 study guide

Read `numerics.py` first: `ode_step` shows the individual recurrences; `solve`
handles grids and reference comparisons; `analyze` handles supplied data; `roots`
handles root iterations. `execute` selects a general operation.

In `ES204_App.pyw`, start at the bottom (`Tk`, `App`, `mainloop`). Then read the
form builders, `run`, `poll`, and `render`. Calculations run in a subprocess and
exchange a JSON request/response. A reader thread puts results in a queue, and
Tk's `after` callback updates the UI on its own thread. Matplotlib figures are
embedded in the window. No formula-specific answers are built into either file.

Practice: add an output column, add a plotting style, or implement another method
with an independent recurrence and an analytic regression test.

This app assumes trusted local formulas and input files. SymPy expression parsing
is not a sandbox for untrusted executable input. The app uses the laptop's existing
Python 3.13 installation with NumPy, SymPy, Tkinter and Matplotlib.

## Unit II-A update

The BVP and Matrix tabs now support linear second-order boundary-value problems, shooting, finite differences, and Direct/Jacobi/Gauss-Seidel/SOR solvers. Show equations generates problem-specific formulas; Export equations writes LaTeX. Save session and Load session preserve inputs, current results, datasets, and table history.

See [BVP_GUIDE.md](BVP_GUIDE.md) for boundary inputs, lecture examples, stopping criteria, and limitations.


## General systems and initial-value finite differences

The **Systems**, **FDM IVP**, and **Compare** tabs add coupled first/second-order
ODEs, finite-difference initial-value recurrences, supplied-reference comparisons,
and detailed equations/stage tables. BVP also accepts reference tables and plots
y with its first derivative. See [SYSTEMS_GUIDE.md](SYSTEMS_GUIDE.md) for input
formats, starting schemes, limits and session-preservation instructions.
