# Noble Gas TDHF Benchmark

## Current TDHF Route

The current prediction route is PySCF TDHF oscillator strengths passed through
the oscillator-channel `alpha(i xi)` and C6 backend.

For the noble gases tested so far:

```text
atom   setting                         C6_tdhf      C6_error
Ne     aug-cc-pV5Z, nstates=300        5.47617368   -14.17%
Ar     aug-cc-pVQZ, nstates=200        60.73027908  -5.55%
Kr     aug-cc-pVQZ, nstates=200        121.85826194 -5.97%
```

## Ne Convergence Check

Ne was checked with `aug-cc-pVTZ`, `aug-cc-pVQZ`, and `aug-cc-pV5Z` for
`nstates = 100, 200, 300`. The `200 -> 300` change is small:

```text
basis         nstates  C6_tdhf     C6_error
aug-cc-pVQZ   200      5.45025099  -14.57%
aug-cc-pVQZ   300      5.46498488  -14.34%
aug-cc-pV5Z   200      5.43116779  -14.87%
aug-cc-pV5Z   300      5.47617368  -14.17%
```

This indicates that the Ne C6 underestimation is not caused by insufficient
`nstates`; it is a basis/method-level systematic error for the current TDHF
route.

## Interpretation

Ar and Kr establish the current TDHF route as a usable long-range C6 baseline at
roughly 5-6% error. Ne is less accurate, with a converged C6 error near -14%.
Future work should investigate whether Ne needs a different response treatment,
larger/more specialized basis sets, or reference-data consistency checks.
