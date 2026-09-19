# Large runs and final-value error tables

Use this workflow when a problem asks only for n, the final estimate, and its error.

1. Select **Output mode: Final values only**.
2. Select **Grid mode: Powers of ten n=10^p**.
3. Set **First p** and **Last p**, for example 1 and 6. This tests 10 through
   1,000,000 steps. Each selected method stops independently.
4. Enter the ODE, initial condition, final x and selected methods as usual.
5. Choose **Final true error %**, enter the tolerance in percent, and select
   **<** for a strict condition or **<=** for an inclusive condition.
6. Enter your reference explicitly: an exact formula, or an xf,y row in the
   given-value table. **Final-values reference source** determines which is used
   when both are present. A formula is used only if you enter and select it.
7. Solve and choose **Final-value trials**. Columns include n, final y, signed
   true error, absolute true error, successive final error and stopping status.

The h/n input and Last m are not used in powers-of-ten mode. Maximum n for
candidate n-search is also unrelated to this mode. Use Last p to bound the run.
Powers-of-ten mode automatically uses compact final-value output.

The signed error is `(reference-final y)/reference*100`. Stopping uses the absolute
value. A tolerance of 10^-5 percent must be entered as `1e-5`, not as a fraction
of a percent converted again. `No early stop` computes every requested grid.
`Final successive error %` compares the current final estimate with the previous
tested grid; it is not the true error.

## Reference precision

A decimal reference can be an exact decimal or a rounded measurement/value; the
app cannot decide which. If rounding to its last displayed digit could exceed
the requested tolerance, a warning describes that uncertainty. The value is
still used unchanged. The app does not infer an analytical solution or add digits.
If a more precise reference is available, explicitly enter that value or formula
and rerun. If not, report that the requested accuracy cannot be certified from
the supplied reference precision alone.

## Limits and existing workflows

- Final-values mode uses constant storage per grid and supports up to 10 million
  steps per grid, with 50 million total steps per run. Cancel remains available.
- Full tables retain the previous node-storage limit and all previous
  differentiation, interpolation and plotting workflows.
- Final-only runs provide convergence plots but no y-versus-x node arrays for
  post-processing. Run a chosen grid in Full tables mode when those are needed.
- No convergence is claimed if the last tested grid fails the chosen condition.
- Previously saved input files still load; missing new fields get their defaults.
- Candidate n-search now also honors the selected strict/inclusive comparison;
  its monotonic-error assumption and smaller n limit are unchanged.

## AI231 source

`final_values.py` contains the compact loop, reference-precision check and stopping
comparison. It calls the existing independent recurrences in `numerics.py`.
`ES204_App.pyw` adds the controls and displays reference warnings.
`test_final_values.py` checks large runs, strict equality, rounding warnings,
reference selection and agreement with the original full-table calculations.
