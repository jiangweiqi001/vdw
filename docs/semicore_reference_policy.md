# Semicore C6 Reference Policy

## Purpose

This policy defines how Zn2/Cd2 semicore-correction results should be judged.
It exists because the numerical conclusion depends on both:

- which published C6 reference is used; and
- whether the current finite-basis Sternheimer result is judged raw or with a
  TRK-normalized diagnostic.

No Zn/Cd result should be called `go` from raw finite-basis data alone.

## Reference Sets

The multi-reference audit table is `reference_pair_c6_alternates.csv`.

| System | Reference Set | C6 (a.u.) | Role |
|---|---:|---:|---|
| Zn2 | Group XII present calculation | 225 | Alternate, low reference |
| Zn2 | Group XII measurement | 257.5 | Primary experimental/spectroscopic-scale reference |
| Zn2 | Benchmark atom data | 276 | Alternate atom-data benchmark |
| Zn2 | QDO/CCSD(T) table | 359 | High-reference stress test |
| Cd2 | Group XII present calculation | 493 | Primary theory reference |
| Cd2 | Group XII measurement | 466 | Primary experimental/spectroscopic-scale reference |
| Cd2 | QDO/CCSD(T) table | 686 | High-reference stress test |

## Primary Decision Target

For this week's proof-of-concept, the primary decision target is the Group XII
present/measurement scale:

- Zn2 primary window: compare against `257.5` first, with `225` and `276` as
  alternate checks.
- Cd2 primary window: compare against both `493` and `466`.

Reason:

- The task is about a missing semicore dispersion correction for Group XII
  weakly-bound dimers.
- The Group XII table directly reports Zn-Zn and Cd-Cd dispersion coefficients
  and measurement comparisons.
- The QDO/CCSD(T) table values (`359/686`) are retained as high-reference stress
  tests because they imply a different, stricter go/no-go conclusion.

## Required Diagnostic Variants

Every Zn/Cd go/no-go statement must report both variants from
`results/semicore_zn_cd_validation/diagnostic_summary.csv`:

1. `raw_finite_basis`
   - Current finite-basis Sternheimer channels as exported.
   - Useful for pipeline and adapted-basis sensitivity.
   - Not sufficient by itself for a go decision.

2. `trk_normalized`
   - Diagnostic scaling of selected d-shell core oscillator strengths to the
     d10 TRK count (`sum_osc = 10`).
   - Not the final continuum response either.
   - Required as a guard against finite-basis oscillator-strength
     over-saturation.

## Go / Review / No-Go Rules

Use these rules until a production continuum/CPHF implementation exists.

`go_candidate`:

- raw finite-basis corrected C6 is within about `10%` of the primary reference;
- TRK-normalized corrected C6 is also within about `10%` of the same primary
  reference;
- correction has the right sign and improves PSP-only;
- no-double-count audit passes;
- basis sensitivity does not change the conclusion.

`review`:

- raw finite-basis looks good, but TRK-normalized does not;
- result depends strongly on adapted basis choice;
- all-electron closure is only a candidate or proxy;
- high-reference stress test gives a conflicting conclusion.

`no_go_or_rework`:

- correction worsens PSP-only against the primary reference in both raw and
  TRK-normalized variants;
- double-counting audit fails;
- core oscillator diagnostic is pathological and no continuum fix is available.

## Current Classification

Zn2:

- raw finite-basis is reference-sensitive.
- TRK-normalized is close to `257.5` and `276`, but not to `225` or `359`.
- Current classification: `review`.

Cd2:

- raw finite-basis is promising against Group XII present/measurement:
  `+3.1%` vs `493`, `+9.1%` vs `466`.
- TRK-normalized is not promising:
  `-26.0%` vs `493`, `-21.7%` vs `466`.
- Current classification: `review`, not `go_candidate`.

## Cd Closure Proxy Decision

The local Cd `ano-rcc` all-electron TDHF calculation is useful as a diagnostic
control, but it is not sufficient as the primary closure reference:

- `ano-rcc` TDHF gives `alpha0 = 86.29549124`, much larger than the standard
  Cd static polarizability scale of about `46` a.u.
- `ano-rcc` TDHF gives `C6 = 811.85667754`, above the `686` high-reference value
  and far above the Group XII present/measurement values (`493/466`).

Therefore, Cd closure should be reported using literature C6 values as proxy
all-electron targets, while keeping `ano-rcc` TDHF as a local diagnostic.

The closure proxy table is:

`results/semicore_zn_cd_validation/cd_literature_closure_proxy.csv`

Current readout:

- Against Group XII `493/466`, raw finite-basis correction closes about
  `107-122%` of the PSP-to-reference gap, while TRK-normalized correction closes
  only about `42-48%`.
- Against the `686` high-reference value, raw finite-basis correction closes
  about `57%`, while TRK-normalized correction closes about `23%`.
- Against local `ano-rcc` TDHF, raw finite-basis correction closes about `44%`,
  while TRK-normalized correction closes about `17%`.

This reinforces the current `review` classification: the raw correction is large
enough to repair the Group XII gap, but the TRK-normalized diagnostic is not.

## Consequence

The next methodological priority is not to optimize the raw finite-basis result.
It is to replace or validate the finite-basis Sternheimer response with a
continuum-aware core response whose oscillator strength passes the selected-shell
sum-rule checks without ad hoc scaling.
