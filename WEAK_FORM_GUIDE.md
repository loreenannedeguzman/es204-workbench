# Weak forms

Open **Weak form**, enter the equation coefficients, interval and boundary
conditions, then click **Derive weak form**. No grid or numerical initial
guess is required. This is a symbolic derivation, not an equation solver.

## Supported input

- **Standard:** `a(x)y'' + b(x)y' + c(x)y = f(x)`.
- **Divergence:** `-(a(x)y')' + b(x)y' + c(x)y = f(x)`.
  Here a is the flux coefficient often called p. Include the displayed leading
  minus sign when translating your original equation to this input form.

Coefficients depend on x only; use the same sin, cos, exp, ln, sqrt and power
syntax as the other workbench tools. Endpoints and boundary coefficients must
be finite real constants. Boundary values and nonzero Dirichlet data are allowed.
The workflow assumes regular smooth real coefficients on the closed interval
and a nonvanishing leading coefficient. Known singular/degenerate and piecewise
inputs are rejected; inconclusive symbolic domain checks produce warnings.

At each endpoint use `alpha*y + beta*y' = gamma`:

- beta=0: Dirichlet, an essential condition. The trial function takes the given
  value; the test function is zero at this endpoint.
- alpha=0: Neumann, a natural condition. The test value is free here.
- Both nonzero: Robin, a natural condition. The y-dependent endpoint term
  remains in the bilinear form, and its known part moves to the load functional.

All derivatives use increasing x. For an outward-normal derivative, the normal
is -1 at the left and +1 at the right; incorporate that sign in beta. Convert
prescribed fluxes to this boundary convention using the endpoint coefficient.

## Derivation

The steps show the original strong equation, boundary conditions, trial/test
spaces, weighted residual, integration by parts, coefficient differentiation,
endpoint substitutions, bilinear form, load functional and final weak identity.

For standard form, integration by parts produces the volume integrand

```text
-a*y'*v' + (b-a')*y'*v + c*y*v
```

and the boundary bracket `[a*y'*v]_left^right`. The derivative of a is included.
For divergence form, the volume integrand is

```text
a*y'*v' + b*y'*v + c*y*v
```

and the bracket is `-[a*y'*v]_left^right`.

Both forms retain the signs of the entered equation. Multiplying the entire
weak identity by -1 gives an equivalent form; a different overall sign from
course notes is not by itself an error. No integrating factor or other scaling
is silently introduced.

The boundary table separately lists coefficients of endpoint y*v in A and
endpoint v in L. Natural derivative data are substituted before stating the
weak problem in H^1; derivatives of arbitrary H^1 trial functions are not
imposed as classical boundary traces.

## Visually rendered equations

The **Rendered equations** results tab opens automatically after derivation.
It displays each formula as typeset mathematical notation with its explanation,
using local Matplotlib math rendering. It needs neither internet nor a LaTeX
installation. Scroll vertically for the steps and horizontally for long formulas.
Maximize the window when reading a long derivation.

**Report / formulas** retains searchable/copyable text and LaTeX. **Export
equations** produces a standalone LaTeX document with both explanations and
equations for this workflow. **Save session** preserves the inputs and symbolic
results; the visual formulas are recreated on load.

## Limits

This release covers scalar, linear, one-dimensional second-order equations.
It does not derive nonlinear, multidimensional, fourth-order or interface weak
forms; it does not build a finite-element mesh or stiffness matrix. A formal
weak derivation does not prove existence, uniqueness or coercivity. In particular,
pure Neumann problems may require compatibility and normalization conditions.

The displayed spaces encode essential boundary conditions and homogeneous test
traces. Nonzero essential values make the trial set affine, not a vector space.

## Validation

The tests substitute manufactured exact solutions into the generated weak
identities and check them by independent symbolic integration for combinations
of Dirichlet, Neumann and Robin data. They also check equivalence of divergence
and expanded standard formulations, coefficient product rules, invalid domains,
decimal-zero boundary entries, rendering, worker execution, session restoration
and explanatory equation export. Existing numerical-method tests are retained.

Installations do not close running windows. Save the current session before
reopening the workbench to access the new tab; previous tab indices are unchanged.
