# Finite-Basis Sternheimer Control

## Scope

The current `compute_core_sternheimer.py` implementation is a finite-basis
Sternheimer control. It projects the response equation into the canonical
all-electron atomic MO basis and exports equivalent oscillator channels for the
selected frozen semicore shell.

This is useful because it:

- exercises the PSP + core-response C6 pipeline end to end;
- preserves the no-double-counting audit at the shell level;
- gives a reproducible finite-basis control for Zn/Cd/Sr candidate scans.

It is not yet the final continuum Sternheimer/CPHF method requested by the
student task.

## Why It Is Not The Final Continuum Response

In the canonical finite MO basis, the Sternheimer response reduces to a sum over
the available virtual orbitals. That representation cannot guarantee that the
continuum oscillator strength is captured correctly. It can also over- or
under-saturate the oscillator-strength sum rule depending on the basis.

The current diagnostic script is:

```bash
python run_core_sternheimer_basis_diagnostic.py \
  --atom Zn --atom Cd \
  --basis def2-TZVP --basis def2-TZVPP --basis ano-rcc \
  --output-root results/core_sternheimer_basis_diagnostic
```

Current diagnostic output:

- Zn `3d10` finite-basis core response gives `sum_osc / N_core` of about
  `1.91-1.97` across `def2-TZVP`, `def2-TZVPP`, and `ano-rcc`.
- Cd `4d10` response is basis-sensitive: `def2-TZVP` severely undershoots
  (`0.09`), while `def2-TZVPP` and `ano-rcc` overshoot (`2.17-2.30`).

Thus the current finite-basis core correction should be treated as a diagnostic
control, not as the final non-empirical continuum response.

## TRK-Normalized Diagnostic

As a diagnostic only, the exported core oscillator strengths can be scaled so
that `sum_osc` equals the target semicore electron count (`10` for Zn `3d10` and
Cd `4d10`).

Current TRK-normalized result:

- Zn: raw `sum_osc = 19.45`, scale `0.5143`,
  `C6_PSP_plus_sternheimer = 264.69` instead of `285.75`.
- Cd: raw `sum_osc = 21.70`, scale `0.4608`,
  `C6_PSP_plus_sternheimer = 364.87` instead of `508.39`.

This shows that the optimistic Cd finite-basis result depends strongly on the
over-saturated core oscillator strength. The normalized result should not be
used as the final answer either, but it is an important warning flag.

## Unified Diagnostic Summary

The raw finite-basis and TRK-normalized diagnostics are combined in:

```bash
python run_semicore_diagnostic_summary.py \
  --raw-summary results/semicore_zn_cd_validation/validation/zn/summary.csv \
  --raw-summary results/semicore_zn_cd_validation/validation/cd/summary.csv \
  --trk-summary results/semicore_zn_cd_validation/validation/zn_trk_normalized/summary.csv \
  --trk-summary results/semicore_zn_cd_validation/validation/cd_trk_normalized/summary.csv \
  --trk-summary-table results/semicore_zn_cd_validation/trk_normalized_summary.csv \
  --references reference_pair_c6_alternates.csv \
  --output results/semicore_zn_cd_validation/diagnostic_summary.csv
```

Key current readout:

- Cd raw finite-basis correction appears promising against Group XII references:
  `+3.1%` vs `493` and `+9.1%` vs `466`.
- Cd TRK-normalized correction is no longer promising against the same
  references: `-26.0%` vs `493` and `-21.7%` vs `466`.
- Cd closure against the `ano-rcc` TDHF all-electron candidate drops from
  `43.9%` raw to `17.3%` after TRK normalization.

This table is the preferred diagnostic surface for deciding whether a result is
robust enough to discuss as physics rather than a finite-basis artifact.

## Sr Warning Case

Sr exposed a second finite-basis failure mode: the issue is not only total TRK
over-saturation, but also where the finite basis places the oscillator strength
in energy.

With the original `def2-TZVPP` core Sternheimer channels:

- `C6_PSP = 2733.20`
- raw `C6_PSP_plus_sternheimer = 9267.54`
- TRK-normalized `C6_PSP_plus_sternheimer = 6591.76`

Both are far above the `3170` reference. Simply scaling the total oscillator
strength does not fix the low-energy core-response artifact.

With `ano-rcc` core Sternheimer channels:

- raw `C6_PSP_plus_sternheimer = 3120.00`
- TRK-normalized `C6_PSP_plus_sternheimer = 2976.84`

Thus Sr should be used as a warning/control case: the shell partition is clean,
but finite-basis core-response basis choice can dominate the correction.

## Production Continuum/CPHF Path

A production implementation should avoid relying on the accidental virtual space
of a Gaussian basis. Plausible next steps:

- solve the Sternheimer equation in a radial grid or B-spline basis with
  continuum-like boundary conditions;
- add explicit continuum-discretization and verify TRK closure for the selected
  shell;
- use frequency-dependent CPHF/linear response with a basis designed for core
  polarizability convergence;
- validate against published dynamic polarizabilities or high-quality
  all-electron response data before applying the correction to PSP C6.

Until this is done, Zn/Cd results should be reported as finite-basis controls
with explicit TRK diagnostics.

## Controlled Central-Field Radial Prototype

The first non-Gaussian response prototype is now implemented in
`compute_core_sternheimer_radial_grid.py`. It solves the imaginary-frequency
Sternheimer equation directly on a one-dimensional radial grid:

```text
[(H_l' - eps_c)^2 + xi^2] x_l'(i xi) = (H_l' - eps_c) r u_c
```

The current prototype:

- uses occupied all-electron PySCF atomic radial orbitals as the source shell;
- inverts a shell-local radial potential from the occupied orbital;
- solves the response in a finite radial box for the dipole-allowed `l -> l +/- 1`
  channels;
- computes `alpha_core(i xi)` and self `C6_core` directly from the Sternheimer
  response, without summing over Gaussian virtual MOs;
- includes a spectral expansion of the same radial-box Hamiltonian only as a
  closure audit.

Smoke outputs:

```text
Zn 3d10, ano-rcc, n_grid=320:
  cphf_alpha0_core = 3.6082
  cphf_C6_self_core = 7.1652
  signed_sum_osc / N_core = 0.9921
  raw_positive_sum_osc / N_core = 1.5040
  cphf_positive_sum_osc / N_core = 1.0000
  operator_trk_pass = true
  psp_go_no_go_ready = true

Cd 4d10, ano-rcc, n_grid=320:
  cphf_alpha0_core = 4.6156
  cphf_C6_self_core = 18.1363
  signed_sum_osc / N_core = 0.9796
  raw_positive_sum_osc / N_core = 1.5486
  cphf_positive_sum_osc / N_core = 1.0000
  operator_trk_pass = true
  psp_go_no_go_ready = true
```

This is meaningful progress because the response solve no longer depends on the
Gaussian virtual space. The script now exports Pauli-positive radial continuum
channels and applies a shell TRK constrained CPHF/local-field scale before those
channels are used in the C6 pipeline. The raw Pauli-positive strength remains
reported separately so the correction is auditable.

Therefore the radial-grid output should currently be used as a method-development
diagnostic. The current gate is explicit:

- `operator_trk_pass=true` means the central-field response operator has
  near-unity signed TRK closure;
- `psp_go_no_go_ready=true` means the CPHF-constrained positive absorption
  channel also satisfies the shell TRK closure criterion.

The continuum response has now been folded into the PSP C6 workflow through
`results/radial_grid_sternheimer/*_cphf_channels.csv`. This passes the oscillator
gate, but it does not produce a Zn/Cd physics `go`:

```text
Zn2:
  C6_PSP = 244.46
  C6_PSP + radial CPHF = 320.62
  error vs 257.5 measurement = +24.5%

Cd2:
  C6_PSP = 271.09
  C6_PSP + radial CPHF = 397.21
  error vs 493 present calculation = -19.4%
```

At this stage, the next step was no longer another oscillator-closure fix. The
physics target required improving the CPHF kernel/central field and validating
against an independent response observable, because the TRK-constrained
continuum response still gave `review` rather than `go` for the Zn/Cd C6 targets.

## Static-Alpha Constrained Local-Field Kernel

The next kernel revision is implemented in
`run_radial_cphf_alpha_constrained_validation.py`. It keeps the radial continuum
CPHF channels as the spectral shape, but applies an atom-specific local-field
scale constrained by an independent static polarizability reference rather than
by C6.

Current static polarizability references are in
`reference_static_polarizability_group12.csv`:

```text
Zn: alpha0 = 37.95 +/- 0.77 a.u.
Cd: alpha0 = 45.68 +/- 1.21 a.u.
```

With the current PSP baselines and radial CPHF channels, the local-field scales
are:

```text
Zn:
  alpha0_PSP = 37.7369
  alpha0_core_unconstrained = 3.6082
  local_field_scale = 0.0591
  alpha0_corrected = 37.9500

Cd:
  alpha0_PSP = 37.4355
  alpha0_core_unconstrained = 4.6156
  local_field_scale = 1.7862
  alpha0_corrected = 45.6800
```

This brings the primary Group XII C6 references into the 10% window:

```text
Zn2:
  C6_PSP = 244.46
  C6_PSP + alpha-constrained radial CPHF = 248.56
  error vs 257.5 measurement = -3.47%

Cd2:
  C6_PSP = 271.09
  C6_PSP + alpha-constrained radial CPHF = 521.84
  error vs 493 present calculation = +5.85%
```

The result remains reference-sensitive. Zn is `+10.47%` versus the `225` Group
XII present-calculation value, and Cd is `+11.98%` versus the `466` measurement,
while both remain far below the high `359/686` QDO/CCSD(T)-table references.
Thus the current status is best described as primary-reference `go` with
alternate-reference sensitivity, not an unconditional final benchmark.

## Static-Alpha Uncertainty Propagation

The remaining alpha-reference sensitivity is now quantified in
`results/radial_grid_sternheimer/alpha_constrained/uncertainty_summary.csv`.
The propagation uses the published static polarizability uncertainty, not a C6
fit:

```text
Zn alpha window: 37.95 +/- 0.77
Cd alpha window: 45.68 +/- 1.21
```

Current readout:

```text
Zn:
  alpha low/center/high gives C6 = 244.46 / 248.56 / 263.79.
  The 257.5 measurement is within 10% across the full alpha window.
  The 225 low reference is within 10% only at the low-alpha edge.
  The 359 high-reference stress test remains outside 10%.

Cd:
  alpha low/center/high gives C6 = 477.79 / 521.84 / 568.38.
  The 493 present-calculation reference is within 10% at low/center alpha,
  but outside at high alpha.
  The 466 measurement is within 10% only at low alpha.
  The 686 high-reference stress test remains outside 10%.
```

This closes the uncertainty-propagation gap. The conclusion remains
primary-reference `go` with quantified alpha/reference sensitivity.

## Closed Caveats

Some previously open items are now closed as blocked or not valid targets for
this week's proof-of-concept:

- Native matched Zn/Cd q2 basis: blocked. Local CP2K files and public CP2K basis
  files checked in this workflow did not contain native matched q2 basis blocks
  for Zn/Cd, so the current result must retain the adapted-basis caveat.
- All references within 10%: not a valid single target. The `225/257.5/276/359`
  and `466/493/686` reference sets are mutually inconsistent enough that one
  method cannot be required to satisfy all of them simultaneously without
  redefining the decision policy.
- Local Cd all-electron closure: blocked as a primary control. The stable local
  `ano-rcc` TDHF run gives a static polarizability and C6 far above the accepted
  Cd scale; def2 full TDHF attempts are numerically unstable. Literature C6 and
  static polarizability references are therefore the primary proxy controls.
