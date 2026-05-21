# Open Issues And Action Plan

This document tracks the parts of issue #1 that are not fully solved yet. It is
intended to drive the next implementation steps in order.

## Current Resolution Status

Resolved for the current prototype:

- All-electron KS-LDA/PBE TDDFT reference leg.
- True PSP-valence TDDFT/RPA baseline, not all-electron masking.
- PSP vs all-electron missing-C6 summary.
- Double-counting guard for explicit valence shells.
- Clean Mg q2 benchmark.
- Diagnostic Ca q2 adapted benchmark.
- Clear labeling that TDHF projections and radial shell-average routes are
  diagnostics, not final EFT-core corrections.

Not fully resolved:

- Rigorous derivation of `Delta alpha_core^EFT(i xi)` from PRL Wilson
  coefficients.
- Final `l=1` dipole Wilson coefficient, derived from EFT rather than from
  all-electron MO transition dipoles.
- Matched Ca/Sr large-core q2 benchmark.
- Screened valence interaction `W_v`.
- Full screened EFT-vdW energy or periodic implementation.

## Priority 1: Close The EFT-Core Derivation Gap

### Problem

The current code has two diagnostic approximations:

```text
EFT_CORE_SCALAR_PROXY
EFT_CORE_DIPOLE_WILSON_MO_APPROX
```

The scalar proxy uses the PRL bandwidth Wilson coefficient but is not a dipole
vdW response. The dipole approximation uses actual `l=1` transition dipoles but
does not yet derive them from the PRL Wilson coefficient.

### What Must Be Answered

1. Does the finite-imaginary-frequency response use the same coefficient as the
   bandwidth `z_core` derivative?
2. What is the correct `l=1` multipole-resolved Wilson coefficient?
3. Which static pieces are already absorbed into the PSP?
4. Is `alpha_cv` already included in PSP-RPA, or does it need an explicit EFT
   cross term?
5. What is the controlled small parameter for C6, where relevant `xi` may not be
   as small as the valence Fermi scale?

### Deliverable

```text
docs/eft_core_alpha_derivation.md
```

Update it from a working note to a derivation with a clearly stated first
implementation formula.

### Acceptance Criteria

- Explicit formula for `Delta alpha_core^EFT(i xi)`.
- Explicit double-counting subtraction rule.
- Clear distinction between scalar bandwidth Wilson coefficient and dipole
  vdW Wilson coefficient.
- Clear statement of what part, if any, approximates `alpha_cv`.

## Priority 2: Implement A Real EFT Dipole Wilson Coefficient

### Problem

The current `compute_dipole_wilson.py` computes all-electron MO dipole
transition channels:

```text
source = EFT_CORE_DIPOLE_WILSON_MO_APPROX
```

This is useful but not yet the EFT Wilson coefficient derived from PRL
`f_K^c`.

### Next Implementation Step

Implement a dedicated dipole Wilson builder:

```text
compute_dipole_wilson_from_core_form_factor.py
```

Target:

```text
atom, shell, Delta_E_Ha, l, m, dipole_weight, source
```

Start with:

```text
Mg 2s,2p
Ca 3s,3p
```

### Acceptance Criteria

- Does not use all-electron TDHF oscillator strengths as the correction.
- Uses atomic core/semi-core orbitals and an EFT-derived multipole form factor.
- Reduces to a backend-compatible oscillator channel.
- Gives the right sign for `C6_EFT - C6_PSP`.

## Priority 3: Preserve And Strengthen The Clean Mg q2 Benchmark

### Current Strongest Result

```text
Mg q2:
C6_PSP = 638.6202
C6_PSP + l=1 dipole = 647.6079
C6_all-e_PBE = 647.5881
double counting = clean
```

### Why It Matters

This is the first clean large-core benchmark where:

```text
explicit valence = 3s
EFT correction shells = 2s,2p
```

There is no shell overlap with the PSP valence space.

### Next Steps

1. Add `results/mg_q2_clean_benchmark.md`.
2. Record all inputs:
   - PSP
   - basis
   - active shells
   - correction shells
   - all-electron reference
3. Add a small reproduction command block.
4. Add a caveat that the current dipole correction is still an MO approximation.

### Acceptance Criteria

- Anyone can reproduce the Mg q2 row.
- Double-counting audit is shown.
- It is clear why this is stronger than the Ca adapted case.

## Priority 4: Find A Matched Ca Or Sr q2 Large-Core Route

### Current Ca Status

The promising Ca row is:

```text
Ca_q2_PBE_adapted:
C6_PSP = 1496.3122
C6_PSP + EFT = 1658.0317
C6_all-e_PBE = 2206.7588
closure = 22.76%
double counting = clean
```

But it uses:

```text
GTH-PBE-q2 pseudo + TZV2P-MOLOPT-PBE-GTH-q10 basis
```

This is a pseudo-basis mismatch and must remain diagnostic.

### Next Steps

1. Search for matched Ca q2 basis data.
2. If unavailable, search Sr q2:
   - explicit valence `5s`
   - EFT correction shells `4s,4p`
3. If still unavailable, document that the project needs external PSP
   generation or DFTK/CP2K support.

### Acceptance Criteria

- A matched large-core Ca or Sr PSP/basis pair runs SCF + TDDFT.
- The corrected semicore shells are absent from explicit PSP valence.
- Double-counting status is `clean`.

## Priority 5: Compare Against Prior-Art Core Polarization Benchmarks

### Motivation

Prior art already has PSP + core-polarization corrections. The novelty here is
not the existence of a core correction; it is the first-principles Wilson
coefficient and cross-observable connection to bandwidth `z_core`.

### Targets

Use prior-art numbers where possible:

- alkali `C6` and core fractions from Derevianko/Johnson/Safronova/Babb
- `alpha_c + alpha_cv + alpha_v` decomposition from Dutt et al.
- frozen-core polarizability failures from Fowler-Sadlej

### Deliverable

```text
results/prior_art_validation_targets.csv
```

### Acceptance Criteria

- The project states which prior-art benchmark each new EFT-core result targets.
- The Mg/Ca clean benchmark is not overclaimed as a first-of-its-kind core
  polarization idea.

## Priority 6: Only Then Start Screened `W_v`

### Do Not Start Yet

Do not implement:

```text
W_v = (v^-1 - chi_v)^-1
```

until the PSP+EFT-core correction is meaningful and double-counting-clean.

### Future Deliverables

```text
screened_pairwise_vdw.py
finite_system_log_vdw.py
```

### Acceptance Criteria

- Screened model is tested first on a known clean benchmark.
- Bare-limit result reproduces current C6-based tail.

## Immediate Next Step

Do these in order:

1. Write `results/mg_q2_clean_benchmark.md`.
2. Upgrade `docs/eft_core_alpha_derivation.md` from working note to explicit
   first-implementation ansatz.
3. Search for matched Ca/Sr q2 basis route.
4. Only after that, implement the true dipole Wilson builder.

The current most important message:

```text
Mg q2 is the strongest clean benchmark.
Ca q2 adapted is promising but diagnostic due to pseudo-basis mismatch.
The final EFT-core derivation remains open.
```
