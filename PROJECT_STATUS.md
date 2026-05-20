# EFT-vdW Prototype: Project Status

## 1. What This Project Is

This repository implements an EFT-inspired prototype for computing long-range
van der Waals response from oscillator-channel data.

The central pipeline is:

```text
Delta_lambda, f_lambda^osc
    -> alpha(i xi)
    -> C6
    -> E(R) = -C6 / R^6
```

The long-term goal is to build a screened EFT-vdW framework in which frozen-core
or semicore electronic dynamics are treated as explicit response channels,
rather than being absorbed only into a static pseudopotential.

At the current stage, this is not a complete screened EFT-vdW functional. It is
a working response-to-C6 prototype with TDHF/RPA oscillator inputs and
core/valence/semicore decomposition diagnostics.

## 2. Physical Idea

Standard frozen-core pseudopotential calculations keep core electrons only
through a static effective potential. Their dynamic response is usually not
explicitly present in the valence response function.

This project tests the idea that frozen core or semicore electrons can be
represented as dynamic response channels:

```text
alpha_core/semicore(i xi)
    = sum_lambda f_lambda^osc / (Delta_lambda^2 + xi^2)
```

These response channels can then be combined with valence response to compute
long-range dispersion coefficients.

The key diagnostic is whether a valence-only response differs significantly from
a valence-plus-semicore response:

```text
Delta C6^semicore
    = C6[alpha_valence + alpha_semicore] - C6[alpha_valence]
```

## 3. What Has Been Implemented

### 3.1 Core Backend

Implemented:

- `alpha(i xi)` from oscillator strengths
- C6 numerical integration
- single-oscillator calibrated controls
- pairwise long-range tail `E(R) = -C6/R^6`
- CSV channel input/output
- reference comparison scripts
- unit tests for the response pipeline

### 3.2 Oscillator Input Routes

| Route | Status | Role |
|---|---|---|
| Calibrated single oscillator | Working | Fitted/control baseline |
| 3D MO oscillator strengths | Working | Independent-particle baseline |
| Radial shell-average | Diagnostic only | Has shell/TRK overcounting issues |
| TDHF/RPA oscillator strengths | Working | Current prediction route |

The current reliable prediction route is the PySCF TDHF/RPA oscillator-channel
route.

### 3.3 Noble Gas Benchmarks

TDHF/RPA, mostly `aug-cc-pVQZ`, with converged or near-converged `nstates`:

| Atom | Setting | alpha0 error | C6 error |
|---|---|---:|---:|
| Ne | `aug-cc-pV5Z`, `nstates=300` | about -11.6% | about -14.2% |
| Ar | `aug-cc-pVQZ`, `nstates=200` | about -4.1% | about -5.6% |
| Kr | `aug-cc-pVQZ`, `nstates=200` | about -3.3% | about -6.0% |

Ne convergence checks indicate that the remaining Ne C6 error is not mainly due
to insufficient TDHF states. Ar and Kr establish the current TDHF route as a
usable long-range C6 baseline at roughly 5-6% error.

### 3.4 Core, Valence, And Semicore Decomposition

The project supports projection of the TDHF transition dipole into valence,
core/semicore, and mixing contributions. This is done at the transition-dipole
amplitude level, not by assigning each TDHF excitation to a single orbital.

Current results:

| System | Partition | Relative C6 correction |
|---|---|---:|
| Ar | valence `3s,3p`; core `1s,2s,2p` | about -0.95% |
| Kr | valence `4s,4p`; core through `3d` | about -2.97% |
| Mg | valence `3s`; semicore `2s,2p` | about -1.95% |
| Ca | valence `4s`; semicore `3s,3p` | about -7.83% |

The most important current result is Ca.

For Ca at the TDHF/`cc-pVQZ` level:

```text
C6[4s]           = 2977.54
C6[4s + 3s3p]    = 2761.48
C6[all-electron] = 2759.51
```

Thus the `3s,3p` semicore dynamic response/mixing changes C6 by:

```text
Delta C6^(3s3p) = -216.06
relative effect = -7.83%
```

This indicates that a large-core-like `4s`-only response can overestimate the
Ca long-range C6, and that semicore dynamic mixing can bring the result close to
the all-electron TDHF response.

## 4. What These Results Mean

The Ca result is the clearest current evidence that semicore dynamics can matter
for long-range dispersion.

Importantly, the correction is not a large positive semicore-only C6. The direct
semicore-only C6 is small. The dominant effect is valence-semicore TDHF mixing
or screening, which reduces the valence-only response.

This is consistent with the EFT motivation: frozen core or semicore electrons
should not be treated only as a static potential; their dynamic response can
renormalize the long-range dispersion response.

## 5. What Is Not Implemented Yet

This repository is not yet a complete screened EFT-vdW functional.

Missing pieces include:

1. Valence-screened interaction:

```text
W_v(i xi) = [v^-1 - chi_v^irr(i xi)]^-1
```

The current implementation uses bare long-range C6/R6 tails.

2. Screened pairwise vdW:

```text
E_AB^(2) = -(1 / 2pi) int dxi Tr[
    alpha_A(i xi) T_AB^scr(i xi)
    alpha_B(i xi) T_BA^scr(i xi)
]
```

3. Full log/RPA or MBD-like energy:

```text
E_vdW ~ (1 / 2pi) int dxi Tr ln[1 - W_v chi_c]
```

4. Periodic implementation and forces.

The current code is atomic/dimer long-range focused. It does not yet provide
periodic q-space screening, total-energy corrections for solids, or forces.

## 6. Current PSP + EFT-Core Benchmark Status

The PSP + EFT-core work has moved beyond the original double-counting
diagnostic for Mg.

Current clean Mg q2 chain:

```text
PSP baseline:       Mg GTH-PBE-q2 / TZV2P-MOLOPT-SR-GTH-q2
explicit valence:   3s
EFT shells:         2s,2p
C6_PSP:             638.6202
C6_PSP+dipole_EFT:  647.6079
C6_all-e_PBE_TDDFT: 647.5881
double counting:    clean
```

This makes Mg q2 the first clean benchmark candidate. The correction is still
an unscreened `l=1` MO dipole Wilson approximation, so it should be described as
a benchmark loop rather than the final screened EFT-vdW method.

Ca also has a clean-by-shell-overlap q2 diagnostic:

```text
C6_PSP:             1972.5599
C6_PSP+dipole_EFT:  2149.5021
C6_all-e_LDA_TDDFT: 1982.7700
```

However, the current Ca q2 row uses `GTH-LDA-q2 + cc-pVQZ`, so it should remain
a diagnostic until a matched large-core Ca PSP basis is established.

Recent clean non-q2 and Be consistency checks are useful negative controls:

```text
case        explicit PSP shells  EFT shells       closure
Be q2 PBE   2s                   1s               2.21%
Be q2 LDA   2s                   1s               2.26%
Kr q8 PBE   4s,4p                core through 3d  4.28%
Ca q10 PBE  3s,3p,4s             1s,2s,2p         0.081%
```

These rows pass the shell-overlap audit, but they do not provide a second strong
closure benchmark. The Be LDA consistency check shows that Be fails as a strong
example because the frozen `1s` contribution is intrinsically small, not because
of the PBE/LDA mismatch in the Be q2 diagnostic. Kr q8 and Ca q10 similarly
show that deep-core-only additive corrections are too small for the current
benchmark goal.

The next production-quality benchmark therefore still needs a large-core
alkaline-earth case where the missing shell is semicore, not deep core:

```text
Ca q2: explicit PSP valence = 4s; EFT shells = 3s,3p
Sr q2: explicit PSP valence = 5s; EFT shells = 4s,4p
```

The local CP2K scan found Ca/Sr/Ba q2 pseudos but no matched q2 basis for these
targets. The active task is now to import or construct a clearly labeled
large-core q2 basis route for Ca or Sr.

## 7. Near-Term Roadmap

Next steps:

1. Import or construct a matched large-core q2 basis route for Ca or Sr:
   - target Ca q2: `GTH-*-q2` with explicit `4s`
   - target Sr q2: `GTH-*-q2` with explicit `5s`
   - require a PySCF build/RKS/TDDFT smoke test before benchmark use
2. Re-run the clean benchmark chain for the first successful Ca/Sr route:
   - PSP valence-only TDDFT
   - PSP + `l=1` dipole EFT semicore correction
   - all-electron TDDFT reference with the same XC backend where possible
3. Clean up documentation and make all current benchmarks reproducible.
4. Add plots for:
   - noble gas TDHF benchmark
   - semicore/core correction summary
   - Ca basis sensitivity
5. Add Ca2 long-range tail comparison:
   - all-electron TDHF
   - valence-only
   - valence+semicore
6. Begin screened-pairwise prototype with a model `W_v`:
   - bare
   - dielectric
   - Yukawa/Thomas-Fermi
7. Later: implement finite-system log/MBD energy.
8. Long term: connect to ab initio valence response `W_v` and periodic systems.

## 8. Current Status In One Sentence

This repository currently provides a working atomic response-to-C6 prototype
with a clean Mg q2 PSP -> PSP+dipole-EFT -> all-electron benchmark loop, while
the final screened EFT-vdW functional and production-quality Ca large-core PSP
benchmark remain future work.
