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

## 6. Near-Term Roadmap

Next steps:

1. Clean up documentation and make all current benchmarks reproducible.
2. Add plots for:
   - noble gas TDHF benchmark
   - semicore/core correction summary
   - Ca basis sensitivity
3. Add Ca2 long-range tail comparison:
   - all-electron TDHF
   - valence-only
   - valence+semicore
4. Optionally test K or Sr as shallow-semicore systems.
5. Begin screened-pairwise prototype with a model `W_v`:
   - bare
   - dielectric
   - Yukawa/Thomas-Fermi
6. Later: implement finite-system log/MBD energy.
7. Long term: connect to ab initio valence response `W_v` and periodic systems.

## 7. Current Status In One Sentence

This repository currently provides a working TDHF/RPA oscillator-response
prototype for EFT-inspired long-range C6 calculations, with evidence that
semicore dynamic mixing can change Ca C6 by about 8% relative to a valence-only
partition.
