# Ar Benchmark Summary

## Current Recommendation

Use the TDHF oscillator route with `aug-cc-pVQZ` and `nstates=200` as the
current Ar prediction setting.

```text
alpha0 = 10.64941140 a.u.   error = -4.06%
C6     = 60.73027908 a.u.   error = -5.55%
```

The TDHF nstates convergence in `ar_tdhf_nstates_convergence.csv` shows that
20 and 50 states are insufficient. The 100, 150, and 200 state results are close
for vdW purposes, with 200 retained as the safer benchmark setting.

## Method Comparison

```text
method                 C6_ArAr    error      role
reference              64.30      0.00%      reference
calibrated             64.30      0.00%      fitted/control
EFT-MO aug-cc-pVQZ     76.14      +18.41%    independent-particle baseline
EFT-TDHF aug-cc-pVQZ   60.73      -5.55%     current prediction baseline
```

The independent-particle MO route is useful as a baseline, but it overestimates
Ar C6. The TDHF/RPA oscillator route corrects much of the low-frequency response
shape and is the current best prediction route.

## Radial Route Status

The radial shell-average route remains diagnostic only. The `aug-cc-pVQZ`
radial shell-average result violates the discrete oscillator-strength sum rule:

```text
sum_osc_discrete / N_core = 1.7118
```

The per-shell audit localizes the overcounting primarily to the Ar `3p` shell.
Do not use the radial shell-average basis-convergence table as a prediction
benchmark until the shell TRK issue is fixed.

## Ar2 Long-Range Tail

`ar2_tail_comparison.csv` compares

```text
E(R) = -C6 / R^6
```

for reference, calibrated, MO, and TDHF C6 values at
`R = 8, 10, 12, 15, 20, 30, 40` Bohr. Since this is a pure long-range tail, the
TDHF curve is uniformly 5.55% weaker than reference and the MO curve is uniformly
18.41% stronger.

## Next Step

Generalize the TDHF route to closed-shell noble gases, starting with Ne and then
Kr, before returning to frozen-core EFT decomposition.
