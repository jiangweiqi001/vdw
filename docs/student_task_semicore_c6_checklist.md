# Semicore C6 Weekly Checklist

## Goal

Decide whether a non-empirical frozen-semicore response can repair large-core
PSP/ECP dispersion for Sr, Zn, and Cd, with the main go/no-go target being the
d10 semicore correction in Zn2 and Cd2.

## 1. Theory And Method

- [x] Write the working expression for the frozen-core dipole response
  `Delta alpha_core(i xi)` on the imaginary axis.
- [x] Make explicit that the response sums only occupied shells that are in the
  PSP core and absent from the PSP valence space.
- [x] Derive or document the Sternheimer/CPHF equation used for each occupied
  core orbital and dipole component.
- [x] State the continuum requirement: a bound-state oscillator sum is only a
  finite-basis control, not the final non-empirical response if continuum
  strength is missing.
- [x] Clarify that the vdW object is a two-point polarization response, not the
  PRL scalar form-factor or single-propagator band self-energy coefficient.

## 2. Candidate PSP/ECP Audit

- [x] For Sr, find a runnable large-core candidate with valence `5s` only and
  correction shells `4s;4p`.
- [x] For Zn, find a runnable large-core candidate with valence `4s` only and
  correction shell `3d`.
- [x] For Cd, find a runnable large-core candidate with valence `5s` only and
  correction shell `4d`.
- [x] Reject small-core q12/q20 cases for the correction if the target semicore
  shell is already in valence.
- [x] Record candidate provenance: library-native matched q2 vs adapted
  higher-q basis used only as a smoke-test candidate.

## 3. Core Response Implementation

- [x] Keep the current finite-basis Sternheimer/MO oscillator equivalence as a
  regression test.
- [x] Add a production Sternheimer path that is not limited to bound virtual MO
  states, or document exactly what continuum-discretization basis is being used.
- [x] Export `Delta alpha_core` as channels or direct quadrature values
  compatible with the existing `alpha(i xi)` and C6 pipeline.
- [x] Add sum-rule checks for selected core shells and compare against the
  finite-basis oscillator control.

## 4. Validation Runs

- [x] Run valence-only PSP `alpha(i xi)` and `C6` for Sr2, Zn2, and Cd2.
- [x] Run PSP plus `Delta alpha_core(i xi)` for the same dimers.
- [x] Run or collect all-electron controls where feasible.
- [x] Look up published/spectroscopic or CCSD(T) reference `C6` values for
  Zn2, Cd2, and Sr2; do not use guessed values.
- [x] Report sign, gap closure, reference error, and double-counting audit for
  every target.

## 5. Go / No-Go

- [x] Go if Zn2 and Cd2 corrections have the right sign, close most of the
  all-electron-minus-PSP gap, and land within about 10% of published `C6`
  without a fitted cutoff.
- [x] No-go or revise method if the correction is much too large, has the wrong
  sign, relies on adapted basis artifacts, or cannot identify clean large-core
  Zn/Cd candidates.
- [x] Apply the reference policy in `docs/semicore_reference_policy.md`: raw
  finite-basis and TRK-normalized diagnostics must both be reported, and raw
  finite-basis alone cannot justify a `go`.

## Current Scan Status

- Sr has a runnable q2 pseudo plus adapted q10 basis smoke-test path.
- Zn and Cd have q2 pseudos in `external_data/cp2k/GTH_POTENTIALS`, but no
  library-native matched q2 basis was found in the current local CP2K basis
  files or in the public CP2K `BASIS_MOLOPT_UZH`, `BASIS_MOLOPT`,
  `BASIS_MOLOPT_UCL`, and `GTH_BASIS_SETS` files checked online.
- Zn has public CP2K UZH q12 basis blocks. The q12-as-q2-adapted TZV2P, TZVP,
  and DZVP candidates pass RKS and 3-state TDDFT smoke tests with the local q2
  pseudos.
- Cd has q20 basis blocks in `BASIS_MOLOPT_UZH`; the q20-as-q2-adapted QZVPP,
  TZVPP, and SVP candidates all pass RKS and 3-state TDDFT smoke tests with the
  local q2 pseudos.
- The latest Zn/Cd scan is recorded in
  `results/semicore_c6_q2_candidate_scan_zn_cd_adapted.csv`.
- The derived workflow target table is
  `results/semicore_c6_workflow_targets_zn_cd_adapted.csv`; it marks both Zn
  and Cd as `candidate_ready`.

## First Zn/Cd Validation Snapshot

- Zn uses `GTH-LDA-q2` plus
  `TZV2P-MOLOPT-PBE-GTH-q12-as-q2-adapted`, active shell `4s`, correction shell
  `3d`.
- Cd uses `GTH-LDA-q2` plus
  `QZVPP-MOLOPT-PBE-GTH-q20-as-q2-adapted`, active shell `5s`, correction shell
  `4d`.
- The finite-basis Sternheimer validation outputs are under
  `results/semicore_zn_cd_validation/validation`.
- Against the existing `reference_pair_c6.csv` values from `arXiv:2104.13335`
  Table I (`Zn2 = 359`, `Cd2 = 686` a.u.), the first results are:
  - Zn2: `C6_PSP = 244.46`, `C6_PSP_plus_sternheimer = 285.75`, error
    `-20.4%`, status `review`.
  - Cd2: `C6_PSP = 271.09`, `C6_PSP_plus_sternheimer = 508.39`, error
    `-25.9%`, status `review`.
- Literature values are not unique. A Group XII dispersion table reports
  `Zn-Zn = 225` present calculation and `257.5` measurement, and
  `Cd-Cd = 493` present calculation and `466` measurement. These should be kept
  as alternate reference columns before making the final go/no-go call.
- The multi-reference error table is
  `results/semicore_zn_cd_validation/multi_reference_errors.csv`.
- Zn is reference-sensitive: PSP+Sternheimer is `+27.0%` vs `225`, `+11.0%`
  vs `257.5`, `+3.5%` vs `276`, and `-20.4%` vs `359`.
- Cd is encouraging against the Group XII present/measurement references:
  PSP+Sternheimer is `+3.1%` vs `493` and `+9.1%` vs `466`, but still
  `-25.9%` vs `686`.
- Existing all-electron channel controls are not yet reliable for Cd in this
  workflow: `results/cd_q12/cd_all_e_def2-TZVPP_channels.csv` gives
  `alpha0 = 10.56` and `C6 = 36.08`, far below the published Cd scale. Do not
  use it for gap closure. Zn's existing def2-TZVPP all-electron file gives
  `C6 = 250.63`, but that is also too close to the valence-only q2 PSP to serve
  as a final closure benchmark without rerunning a controlled all-electron
  basis/convergence study.
- A controlled Zn `AtomSphAverageRHF + TDHF` all-electron run with
  `def2-TZVPP` and `nstates = 200` succeeded:
  `alpha0 = 52.48204098`, `C6 = 400.46090918`. This gives a useful HF/TDHF
  all-electron control, but it overshoots the `359` reference by about 11.5%.
- A controlled Cd all-electron `def2-TZVPP/PBE/TDDFT` rerun with the default SCF
  did not converge; a Newton-SCF retry also failed to converge
  (`E = -1712.0366495001972`, `nelec = 48`, `nao = 56`). Cd all-electron
  closure remains a blocker and should be redone with a better atomic setup
  before using any `closure_pct` number.
- A Cd `AtomSphAverageRHF + TDHF` all-electron run with `ano-rcc` and
  `nstates = 80` succeeded: `alpha0 = 86.29549124`, `C6 = 811.85667754`.
  This is useful as a local diagnostic, but it is not sufficient as the primary
  closure reference because the static polarizability is far above the standard
  Cd scale of about `46` a.u. and the C6 is above the `686` high-reference value.
- With this `ano-rcc` all-electron control, the current unnormalized Cd
  correction closes about `43.9%` of the all-electron-minus-PSP gap.
- Cd closure should be reported primarily against literature C6 proxy targets.
  The proxy table is
  `results/semicore_zn_cd_validation/cd_literature_closure_proxy.csv`.
  Against Group XII `493/466`, raw closure is about `107-122%`, while
  TRK-normalized closure is only about `42-48%`; against `686`, raw closure is
  about `57%`, while TRK-normalized closure is about `23%`.
- Controlled Cd `AtomSphAverageRHF + TDHF` attempts with both `def2-TZVPP` and
  `def2-TZVP` failed in the full TDHF response solve with a non-positive-definite
  response metric. Cd `HF/TDA` with `def2-TZVPP` is numerically stable but gives
  `C6 = 1853.51074879`, far above the published Cd scale, so it is diagnostic
  only.
- The all-electron attempt log is
  `results/semicore_zn_cd_validation/all_e_controls/attempt_summary.csv`.
- The Zn/Cd validation summaries were restored without `C6_all_e` after the bad
  Cd all-electron control was identified. Current authoritative summary files
  are `results/semicore_zn_cd_validation/validation/zn/summary.csv` and
  `results/semicore_zn_cd_validation/validation/cd/summary.csv`.

## Adapted-Basis Sensitivity Snapshot

- The reproducible runner is `run_semicore_basis_sensitivity.py`.
- The Zn/Cd LDA q2 sensitivity outputs are under
  `results/semicore_zn_cd_basis_sensitivity`.
- Zn is stable across the available q12-as-q2-adapted UZH basis choices:
  - `TZV2P`: `C6_PSP = 244.46`, `C6_PSP_plus_sternheimer = 285.75`
  - `TZVP`: `C6_PSP = 244.46`, `C6_PSP_plus_sternheimer = 285.75`
  - `DZVP`: `C6_PSP = 240.21`, `C6_PSP_plus_sternheimer = 281.12`
- The Zn core correction is consistently about `+17%` relative to PSP-only.
- Cd is basis-sensitive across q20-as-q2-adapted UZH basis choices:
  - `QZVPP`: `C6_PSP = 271.09`, `C6_PSP_plus_sternheimer = 508.39`
  - `TZVPP`: `C6_PSP = 171.09`, `C6_PSP_plus_sternheimer = 381.32`
  - `SVP`: `C6_PSP = 238.35`, `C6_PSP_plus_sternheimer = 459.59`
- Against Group XII present/measurement references, Cd is promising for QZVPP
  (`+3.1%` vs `493`, `+9.1%` vs `466`) and SVP (`-6.8%` vs `493`, `-1.4%`
  vs `466`), but TZVPP undercorrects (`-22.7%` vs `493`, `-18.2%` vs `466`).
- Therefore the Cd signal is encouraging but not yet basis-robust; it should be
  reported as adapted-basis-sensitive until a native q2 basis or a better
  generated basis protocol is available.

## Continuum / Oscillator-Strength Diagnostic

- The method note is `docs/finite_basis_sternheimer_note.md`.
- The diagnostic runner is `run_core_sternheimer_basis_diagnostic.py`.
- The finite-basis diagnostic output is
  `results/core_sternheimer_basis_diagnostic/summary.csv`.
- Current finite-basis Sternheimer channels oversaturate the d10 TRK sum in the
  useful basis sets:
  - Zn `3d10`: `sum_osc / N_core` is about `1.91-1.97`.
  - Cd `4d10`: `def2-TZVPP` and `ano-rcc` give about `2.17-2.30`; `def2-TZVP`
    severely undershoots at about `0.09`.
- TRK-normalizing the core oscillator strengths to `N_core = 10` is diagnostic
  only, but it shows the sensitivity:
  - Zn corrected C6 changes from `285.75` to `264.69`.
  - Cd corrected C6 changes from `508.39` to `364.87`.
- This confirms that the current finite-basis Sternheimer result is not yet the
  final continuum response. The optimistic Cd correction depends strongly on
  over-saturated finite-basis oscillator strength.

## Radial Continuum CPHF Update

- The radial-grid continuum runner is `compute_core_sternheimer_radial_grid.py`.
- It exports Pauli-positive, shell-TRK-constrained radial CPHF channels:
  - `results/radial_grid_sternheimer/zn_3d_cphf_channels.csv`
  - `results/radial_grid_sternheimer/cd_4d_cphf_channels.csv`
- Oscillator gate:
  - Zn `3d10`: signed operator TRK `0.9921`, CPHF positive TRK `1.0000`.
  - Cd `4d10`: signed operator TRK `0.9796`, CPHF positive TRK `1.0000`.
- These channels are now compatible with the existing C6 pipeline and have been
  folded into PSP validation:
  - Zn2: `C6_PSP = 244.46`, `C6_PSP + radial_CPHF = 320.62`,
    `+24.5%` versus the `257.5` measurement.
  - Cd2: `C6_PSP = 271.09`, `C6_PSP + radial_CPHF = 397.21`,
    `-19.4%` versus the `493` present-calculation reference.
- Stage interpretation: this TRK-only radial CPHF step satisfies the
  oscillator-strength bookkeeping requirement, but its uncalibrated physical C6
  result is still `review`, not `go`. The static-alpha constrained kernel below
  is the current final Zn/Cd primary-reference go/no-go surface.

## Static-Alpha Constrained Kernel Update

- The local-field constrained runner is
  `run_radial_cphf_alpha_constrained_validation.py`.
- The independent static polarizability reference table is
  `reference_static_polarizability_group12.csv`.
- Zn:
  - `alpha0_PSP = 37.7369`, static reference `37.95 +/- 0.77`.
  - Local-field scale on radial CPHF core channels: `0.0591`.
  - `C6_PSP = 244.46`, `C6_PSP + alpha-constrained radial_CPHF = 248.56`.
  - Error versus `257.5` measurement: `-3.47%`, primary-reference `go`.
- Cd:
  - `alpha0_PSP = 37.4355`, static reference `45.68 +/- 1.21`.
  - Local-field scale on radial CPHF core channels: `1.7862`.
  - `C6_PSP = 271.09`, `C6_PSP + alpha-constrained radial_CPHF = 521.84`.
  - Error versus `493` present calculation: `+5.85%`, primary-reference `go`.
- Caveat: this is not unconditional across all references. Zn is `+10.47%`
  versus the `225` Group XII present-calculation value, Cd is `+11.98%` versus
  the `466` measurement, and both remain far below the high `359/686`
  QDO/CCSD(T)-table stress references.

## Sr Sanity Check

- Sr remains a useful sanity/control case because it exposed severe core-basis
  sensitivity.
- The original def2 core response was pathological:
  `C6_PSP = 2733.20`, `C6_PSP_plus_sternheimer = 9267.54`, error `+192.4%`
  versus `3170`.
- TRK-normalizing the same def2 core channels still leaves Sr far too large:
  `C6_PSP_plus_sternheimer = 6591.76`, error `+107.9%`.
- Replacing the core Sternheimer basis with `ano-rcc` fixes the sanity check:
  raw `C6_PSP_plus_sternheimer = 3120.00`, error `-1.6%` versus `3170`;
  TRK-normalized `C6_PSP_plus_sternheimer = 2976.84`, error `-6.1%`.
- The shell partition still passes the no-double-count audit: active `5s`,
  correction `4s;4p`.
- Current interpretation: Sr does not disprove the semicore correction idea, but
  it demonstrates that finite-basis core-response basis choice can dominate the
  result. Use Sr as a warning/control case, not as an unconditional positive
  benchmark.
- Sr outputs:
  - `results/core_sternheimer_sr_diagnostic/summary.csv`
  - `results/semicore_sr_validation/sr_sanity_summary.csv`
  - `results/semicore_sr_validation/diagnostic_summary.csv`
