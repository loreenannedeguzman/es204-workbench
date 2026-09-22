# Unit II-A: boundary-value problems and matrix equations

## BVP inputs

Enter a linear second-order equation in coefficient form:

`a(x)*y'' + b(x)*y' + c(x)*y = d(x)`.

The four coefficient boxes accept formulas of x, including `sin`, `cos`, `exp`, `ln`, `pi` and powers. Do not enter an equals sign. All boundary derivatives mean differentiation toward increasing x.

At each end enter `alpha*y + beta*y' = gamma`:

- Dirichlet y=value: alpha=1, beta=0, gamma=value.
- Neumann y'=value: alpha=0, beta=1, gamma=value.
- Robin: enter all three coefficients as given.
- Mixed: select the appropriate coefficients independently at each end.

Grid choices explicitly distinguish **points n**, **intervals**, and **step h**. The lecture counts points: n points means n-1 intervals. Step h must divide the domain; the program never changes it silently. Current study/report limit: 3–251 points.

## Finite differences

Interior first and second derivatives use central second-order formulas. Boundary derivatives use three-point second-order one-sided formulas. Endpoints are eliminated into an interior matrix when possible, matching the lesson's mixed-boundary example. Otherwise the full matrix is retained.

Choose Direct, Jacobi, Gauss-Seidel, or SOR. Initial vector values follow the unknown names printed in the matrix; a blank vector means zeros. Set tolerance, < or <=, maximum iterations, and one of:

- Residual: max |A*u-b| (not percent).
- Absolute update: max |u_new-u_old| (not percent).
- Relative update %: max |(u_new-u_old)/u_new|*100. A changed component ending at zero is undefined and cannot pass this criterion.

Always inspect the residual even when stopping on update size. A singular/incompatible matrix is reported rather than presented as a unique solution. Iteration failure is explicitly reported. SOR requires 0<omega<2; convergence is not guaranteed just by this range.

## Shooting

The app reduces the ODE to y'=z and z'=(d-b*z-c*y)/a. It applies the selected vector Euler, PC, RK2.1, lesson RK2.2, RK3, or RK4 method. PC passes specifies the number of corrections per step.

For a Dirichlet left boundary, the parameter is the initial slope z(left). For any left boundary containing z, the parameter is y(left), with z(left) calculated from that boundary. Supply two guesses. Secant updates target the right boundary residual alpha*y(right)+beta*z(right)-gamma.

Choose boundary residual or absolute parameter-change stopping; tolerances are not percentages. Parameter convergence alone does not establish a small boundary residual. Tables show both, plus every trajectory and its y/z values. Negligible shooting sensitivity is rejected because uniqueness cannot be established.

## Equations and reports

Click **Show equations** in the BVP or Matrix tab to generate equations without computing a solution. Open **Report / formulas** for the problem-specific derivation; choose the augmented matrix under Tables. After solving, the same formulas remain alongside results.

Reports include governing/boundary equations, grid, derivative stencils, discrete interior relation, substituted row equations, boundary elimination, reduced matrix, chosen iterative updates, shooting IVP conversion, problem-specific RK stage expressions, and secant relation as applicable. They describe the selected method, not every possible algebraic derivation or every method simultaneously.

**Export equations** saves a standalone LaTeX file with the current generated equations. **Save report** saves equations and numerical tables as text. Large equations are fitted to page width in the LaTeX export; use smaller matrices for readable classroom presentation.

## Matrix tab

Enter square A one row per line, and right-side b as comma/space/newline separated values. Choose Direct/Jacobi/Gauss-Seidel/SOR. The tab accepts any square real nonsingular system up to 250 unknowns, not just BVP matrices. Zero diagonal entries prevent the iterative formulas; reorder the system yourself or choose Direct. No silent equation reordering occurs.

## Lecture entry examples (validation examples, not built-in solutions)

1. Example 1: a=1, b=-2, c=1, d=-8+7*x-x^2; interval [0,2]; left (1,0,0), right (1,0,2); Step h=.5; Shooting/RK4. Optional exact: 3*x-x^2.
2. Example 2: a=1, b=0, c=1, d=1; interval [0,pi/2]; left (0,1,-3), right (0,1,-2); Step h=pi/10; Shooting/RK4; guesses 2 and 4. Optional exact: 1+2*cos(x)-3*sin(x). Use the governing equation's domain; a requirement line in the slide has a conflicting endpoint.
3. Example 3: a=5+x, b=-2, c=1, d=8*x; interval [0,1.5]; left (1,0,1), right (1,0,10); Step h=.25; Finite differences.
4. Example 4: a=1, b=0, c=-2, d=cos(x); interval [0,pi/2]; left (1,0,3), right (1,1,0); Step h=pi/8; Finite differences.
5. Matrix example: A rows `-2.30843 1 0`, `1 -2.30843 1`, `0 .73583 -1.25174`; b=`-2.85753 .10905 .05901`. Enter these rounded numbers to solve that printed matrix literally. The unrounded BVP-generated matrix gives slightly different final digits. Even the printed matrix's exact solution differs slightly from some printed slide values; numerical residuals are the check.

## Preserve work

Save inputs keeps only settings. **Save session** keeps inputs, current report/results, datasets and table history; **Load session** restores them. Old version-2 input files still load. Running equation preview replaces the current displayed report, while prior tables and computed datasets remain available. Save a session before replacing any report you want to retain in full.

New tabs require restarting the GUI once after installation. The updater does not close existing windows. In the old window, save inputs and reports before reopening; the newly added full-session button becomes available after reopening.

## Scope

This release covers the lesson's linear second-order scalar BVPs and real matrix iteration methods. It does not add nonlinear BVPs, arbitrary ODE systems, finite elements, complex roots, or unrelated roadmap items. Optional exact solutions are comparisons only; no exact reference is inferred or substituted.
