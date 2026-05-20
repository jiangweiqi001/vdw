# All-Electron RKS-LDA TDDFT Reference: Current Status

## What Was Implemented

Milestone A has been started with a PySCF all-electron Kohn-Sham response leg:

```text
run_all_e_rpa_atom.py
run_all_e_rpa_convergence_ar.py
```

The current implementation runs all-electron `RKS` with `LDA,VWN` or `PBE`, then
PySCF `TDDFT` or `TDA`, exports oscillator-strength channels, and reuses the
existing `alpha(i xi) -> C6` backend.

## Ar Convergence Result

The first convergence scan was run for Ar with:

```text
xc = lda
method = TDDFT
basis = cc-pVTZ, aug-cc-pVTZ, aug-cc-pVQZ
nstates = 50, 100, 150, 200
```

Key rows:

```text
basis        nstates  alpha0      C6        alpha0_err   C6_err
cc-pVTZ      200      6.6500      36.0791   -40.09%      -43.89%
aug-cc-pVTZ  100      11.4197     67.1998   +2.88%       +4.51%
aug-cc-pVTZ  200      11.4610     68.0409   +3.25%       +5.82%
aug-cc-pVQZ  100      11.7709     67.8828   +6.04%       +5.57%
aug-cc-pVQZ  200      11.8181     68.7910   +6.47%       +6.98%
```

## Interpretation

The all-electron KS-LDA TDDFT route is technically working.

For Ar, diffuse basis functions are essential:

```text
cc-pVTZ severely underestimates alpha0 and C6.
aug-cc-pVTZ / aug-cc-pVQZ give C6 within about +5-7% of reference.
```

The `nstates=100` results are already close to the `nstates=200` values for the
diffuse bases, although the 200-state values are retained for consistency with
the TDHF benchmark.

## Ne/Kr/Mg/Ca LDA-TDDFT Checks

The same all-electron RKS-LDA TDDFT path was run for Ne, Kr, Mg, and Ca:

```text
atom  basis        nstates  alpha0      C6         alpha0_err   C6_err
Ne    aug-cc-pVQZ  200      2.8797      7.0935     +7.85%       +11.18%
Kr    aug-cc-pVQZ  200      17.6900     134.0930   +5.30%       +3.47%
Mg    aug-cc-pVQZ  200      70.8902     610.0957   -0.57%       -2.70%
Ca    cc-pVQZ      200      148.4496    1982.7700  -5.51%       -10.73%
```

Ca uses `cc-pVQZ` because PySCF's built-in basis library does not provide
`aug-cc-pVQZ` for Ca in this environment.

These results establish an all-electron KS-LDA TDDFT reference leg for the
current test set. Ne and Kr are on the high side, Mg is close, and Ca is lower
than the all-atom reference.

## PBE-TDDFT Checks

The same single-point checks were repeated with PBE:

```text
atom  basis        nstates  alpha0      C6         alpha0_err   C6_err
Ne    aug-cc-pVQZ  200      2.8774      7.0720     +7.77%       +10.85%
Ar    aug-cc-pVQZ  200      11.7357     68.2301    +5.73%       +6.11%
Kr    aug-cc-pVQZ  200      17.7114     134.6496   +5.42%       +3.90%
Mg    aug-cc-pVQZ  200      73.6424     647.5881   +3.29%       +3.28%
Ca    cc-pVQZ      200      159.6043    2206.7588  +1.59%       -0.64%
```

PBE is broadly similar to LDA for Ne/Ar/Kr, improves Ca substantially, and shifts
Mg from a slight LDA underestimation to a small overestimation. This gives a
basic XC sensitivity check for Milestone A.

## KS-TDDFT vs HF-TDHF Comparison

`results/all_e_rpa_vs_hf_tdhf_summary.csv` compares the KS-LDA/PBE TDDFT route
against the existing all-electron HF/TDHF route on the same atom set.

Key C6 rows:

```text
atom  route          C6        C6_err
Ne    LDA-TDDFT      7.0935    +11.18%
Ne    PBE-TDDFT      7.0720    +10.85%
Ne    HF-TDHF        5.4503    -14.57%

Ar    LDA-TDDFT      68.7910   +6.98%
Ar    PBE-TDDFT      68.2301   +6.11%
Ar    HF-TDHF        60.7303   -5.55%

Kr    LDA-TDDFT      134.0930  +3.47%
Kr    PBE-TDDFT      134.6496  +3.90%
Kr    HF-TDHF        121.8583  -5.97%

Mg    LDA-TDDFT      610.0957  -2.70%
Mg    PBE-TDDFT      647.5881  +3.28%
Mg    HF-TDHF        757.8355  +20.87%

Ca    LDA-TDDFT      1982.7700 -10.73%
Ca    PBE-TDDFT      2206.7588 -0.64%
Ca    HF-TDHF        2759.5097 +24.25%
```

The KS-TDDFT and HF-TDHF routes bracket several references from opposite sides.
For Mg and Ca, KS-PBE TDDFT is much closer to the all-atom reference than
HF-TDHF. This supports using KS-LDA/PBE TDDFT as the all-electron reference leg
for the PSP-RPA comparison.

## Is Milestone A Complete?

Yes for the current prototype scope.

Completed:

- All-electron RKS-LDA TDDFT path exists.
- PBE TDDFT checks were run for the current test set.
- Ar basis/nstates convergence was run.
- Ne, Kr, Mg, and Ca single-point checks were run.
- KS-LDA/PBE TDDFT was compared against the existing HF/TDHF route.
- Output format matches the roadmap.
- Existing oscillator-channel backend is reused.

Remaining optional extensions:

- Add larger-basis checks where PySCF provides the basis.
- Add more atoms when moving beyond the current test set.

Current status:

```text
Milestone A is complete for the current prototype test set.
The next milestone is a genuine PSP-valence RPA baseline.
```
