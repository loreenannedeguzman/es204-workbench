# General 1D finite-element workflow

This edition adds an **FEA** tab to the existing workbench. It contains no
assignment-specific inputs or answers. The previous installed files are preserved in the timestamped backup folder.

## Scope

Scalar, linear, regular, one-dimensional boundary-value equations; two-node
linear Galerkin elements; nonuniform meshes; variable and piecewise coefficients;
piecewise and concentrated nodal loads; Dirichlet, Neumann and Robin conditions.
It does not solve nonlinear equations, 2D/3D meshes, beams requiring fourth-order
elements, or singular leading coefficients.

## Enter the equation

Choose **Self-adjoint** for `-(p*y')' + r*y = f`. Signs matter: a lesson equation
`(p*y')' + c*y + q = 0` uses `r=-c`, `f=q` in this form.

Choose **Standard** for `a*y'' + b*y' + c*y = f`. The program derives
`mu=exp(integral((b-a')/a))`, then uses `p=mu*a`, `r=-mu*c`, and the transformed
load `-mu*f`. The derivation is included in the rendered equations. If an explicit
integrating factor cannot be obtained, supply a self-adjoint equation yourself.

Enter increasing node coordinates separated by commas. Each adjacent pair is an
element, and global node IDs start at 1. Place every material or load interface at
a node. Connectivity is displayed explicitly; arbitrary disconnected meshes are
not supported.

Optional self-adjoint region rows use semicolons:

`left; right; p(x); r(x); f(x)`

Alternatively, forcing-only regions use `left; right; f(x)` and work with either
equation form. Tables must cover the whole domain without gaps or overlaps.
Do not combine the two region tables. In Standard mode the entered forcing is
the original right side, before multiplication by `-mu`.

## Boundary conditions

Enter `alpha*y + beta*y' = gamma` at each endpoint. The derivative is in the
increasing-x direction at both ends. For `y=value`, use `(alpha,beta,gamma) =
(1,0,value)`. For `y'=value`, use `(0,1,value)`.

The workbench derives the outward-flux contribution with normal -1 on the left
and +1 on the right. Robin terms modify both the diagonal and right side. At
internal interfaces the continuous solution and weak flux balance are implicit
in shared-node assembly. Optional concentrated loads add directly to the weak
right side at the chosen node.

## Assemble, then manually partition

1. Leave Workflow at **Assemble only** and click **Run FEA**.
2. Inspect **Tables**, **Plots**, and **Rendered equations**. The equation viewer
   has Previous/Next controls, 24 equations per page. Export equations includes
   every page, not just the visible page.
3. Find **Assembled system for manual partition**. The report identifies free
   and prescribed node IDs. Form `K_FF` and `F_F - K_FD*y_D` yourself.
4. Enter the free-node IDs in your chosen order, your reduced matrix (one row
   per line), and your reduced right side. Choose **Check partition and solve**.
5. The program checks both arrays against the assembly and your node ordering.
   A mismatch stops the solve. On acceptance it solves **your entered arrays**,
   so rounding in your inputs can affect the solution. The default absolute
   entry tolerance is `1e-6`; retain sufficient digits when copying values.
6. Read **FEA nodal solution**, **Element slopes**, the residual report and plots.
   Prescribed-node residuals are reactions. A small free-node residual verifies
   the algebraic solve, not mesh convergence. Refine the mesh to assess accuracy.

No numerical reduced system is supplied by Assemble only. The generic partition
identity is shown to guide manual working. A mesh with all nodes prescribed has
no reduced equations to solve.

## Integration and working

**Gauss 2** uses exactly two Gauss points per element for every integral. Tables
give xi, x, weights, Jacobian, integrand values and weighted contributions.
**Auto** uses symbolic antiderivatives for a bounded supported family (polynomials
and suitable linear-argument sine, cosine and exponential expressions), and
falls back to Gauss 2. Fallback does not assert that an antiderivative is
mathematically impossible. **Analytical only** rejects unsupported integrals.

Rendered working includes the transformation, weighted residual, integration by
parts, local shapes and gradients, all independent local matrix/load integrals,
local-to-global connectivity, matrices, assembly, endpoint contributions, manual
partition identity and reconstructed approximation. Matrices are rendered as
numeric grids; LaTeX export uses matrix environments. Full-precision arrays are
available in the tables and saved session. Rendered values are abbreviated.

Use Save session to preserve inputs, results and table history; existing sessions
without FEA fields remain loadable. Existing CSV, plot and equation exports are
available. Use the solution dataset in Data / analysis for additional processing.

## Run and verify this installation

Run `python ES204_App.pyw` from this directory using the existing Python
installation with NumPy, SymPy, Matplotlib and Pillow (and Tk support).
Use your existing workbench shortcut to open this installation.

Run `python -B -m unittest discover -p "test_*.py"` for numerical regression tests.
GUI, rendering, session and export checks were completed on the staged build before installation.

Installed from the tested workspace build. GitHub was not changed by this installation.

## Termwise Auto integration

Auto expands each matrix/load integrand into additive terms, integrates terms
within its supported symbolic family analytically, and applies two-point Gauss
quadrature only to the remaining terms. Reports show each contribution and their
sum. This applies generally to coefficients and loads, not just one assignment.
Explicit Gauss 2 still integrates the entire integrand with two points.
Analytical only still rejects expressions outside the supported family.
Fallback does not prove that no antiderivative exists.

Existing session inputs, saved results, and open displayed results are not
migrated or recalculated. A new run uses the new integration policy. Previously
entered manual partitions remain unchanged and are validated against the new
assembly; if they no longer match, assemble and update your partition manually.