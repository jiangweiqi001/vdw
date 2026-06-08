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
C6_PSP+dipole_EFT:  647.6996
C6_all-e_PBE_TDDFT: 647.5881
double counting:    clean
```

This makes Mg q2 the first clean benchmark candidate. The correction is still
an unscreened `l=1` MO dipole Wilson approximation, so it should be described as
a benchmark loop rather than the final screened EFT-vdW method.

The best current official-matched secondary candidate is Sn q4:

```text
PSP baseline:       Sn GTH-PBE-q4 / TZV2P-MOLOPT-PBE-GTH-q4
explicit valence:   5s,5p
EFT shell:          4d
C6_PSP:             474.0812
C6_PSP+dipole_EFT:  566.1352
C6_all-e_ANO_TDDFT: 576.9841  (nstates=150)
closure:            89.46%
double counting:    clean
```

Sn q4 is not yet as final as Mg q2 because its all-electron reference is more
delicate. Full TDDFT/ANO is stable from 80 to 150 states at the few-percent
level, but TDA gives a much larger reference C6. Therefore Sn q4 should be
reported as a strong official-matched candidate with a reference-method caveat.

Ge q4 is the best weaker official-matched control:

```text
PSP baseline:       Ge GTH-PBE-q4 / TZV2P-MOLOPT-PBE-GTH-q4
explicit valence:   4s,4p
EFT shell:          3d
C6_PSP:             306.6594
C6_PSP+dipole_EFT:  317.4539
C6_all-e_augQZ:     375.5721
closure:            15.66%
```

The current implementation also includes a finite-system model-screened
pairwise/logdet prototype:

```text
screened_pairwise_vdw.py
screened_eft_vdw.py
```

This supports bare, dielectric, and Yukawa model `W_v` kernels. It verifies that
the second-order logdet bare limit reproduces `-C6/R^6`. This is a model-screened
interface, not an ab initio `W_v = (v^-1 - chi_v)^-1` implementation.

Ca also has a clean-by-shell-overlap q2 diagnostic:

```text
C6_PSP:             1972.5599
C6_PSP+dipole_EFT:  2149.5021
C6_all-e_LDA_TDDFT: 1982.7700
```

However, the current Ca q2 rows remain diagnostic. One route uses
`GTH-LDA-q2 + cc-pVQZ`; the later adapted route uses `GTH-PBE-q2` with the
`TZV2P-MOLOPT-PBE-GTH-q10` UZH basis. This is a pseudo-basis mismatch:
the pseudopotential is q2 but the basis was optimized for q10. It should be
reported only as `diagnostic`, not as a final matched Ca q2 benchmark, until a
proper matched large-core Ca q2 basis is established.

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

The next production-quality alkaline-earth benchmark still needs a large-core
q2 route where the missing shell is semicore, not deep core:

```text
Ca q2: explicit PSP valence = 4s; EFT shells = 3s,3p
Sr q2: explicit PSP valence = 5s; EFT shells = 4s,4p
```

The local and online CP2K scans found Ca/Sr/Ba q2 pseudopotentials but no
reliable matched q2 basis blocks for these targets. Generated-basis work is now
guarded by `docs/generated_basis_protocol.md`, which requires freezing and
hashing a generated basis before any vdW validation.

### 6.1 Group XII Large-Core Semicore Status

The current Zn/Cd work targets the student-task direction: large-core q2
pseudopotentials with explicit outer `ns` valence and missing `d10` semicore
response.

Current candidates:

```text
Zn: GTH-LDA-q2 + TZV2P/TZVP/DZVP-MOLOPT-PBE-GTH-q12-as-q2-adapted
    explicit PSP shell = 4s; correction shell = 3d

Cd: GTH-LDA-q2 + QZVPP/TZVPP/SVP-MOLOPT-PBE-GTH-q20-as-q2-adapted
    explicit PSP shell = 5s; correction shell = 4d
```

No library-native matched q2 basis was found for Zn/Cd in the local CP2K files
or in the public CP2K `BASIS_MOLOPT_UZH`, `BASIS_MOLOPT`,
`BASIS_MOLOPT_UCL`, and `GTH_BASIS_SETS` files checked online. The current
Zn/Cd rows therefore remain adapted-basis diagnostics, not final matched-q2
benchmarks.

The reference policy is documented in `docs/semicore_reference_policy.md`.
For this proof-of-concept, the primary Group XII reference scale is:

```text
Zn2: 257.5 a.u. measurement, with 225 and 276 as alternate checks
Cd2: 493 a.u. present calculation and 466 a.u. measurement
```

The `359/686` QDO/CCSD(T)-table values are retained as high-reference stress
tests. Any go/no-go statement must report both raw finite-basis and
TRK-normalized diagnostics. Raw finite-basis data alone cannot justify a `go`.

Basis sensitivity:

```text
Zn raw finite-basis:
  TZV2P/TZVP: C6_PSP = 244.46, C6_PSP+core = 285.75
  DZVP:       C6_PSP = 240.21, C6_PSP+core = 281.12
  status: basis-stable, reference-sensitive

Cd raw finite-basis:
  QZVPP: C6_PSP = 271.09, C6_PSP+core = 508.39
  TZVPP: C6_PSP = 171.09, C6_PSP+core = 381.32
  SVP:   C6_PSP = 238.35, C6_PSP+core = 459.59
  status: promising against Group XII references in QZVPP/SVP, but
          adapted-basis-sensitive
```

The finite-basis Sternheimer diagnostic shows that the current d10 core channels
over-saturate the TRK sum in the useful bases:

```text
Zn 3d10: sum_osc / N_core about 1.9
Cd 4d10: sum_osc / N_core about 2.2 for def2-TZVPP/ano-rcc
```

TRK-normalized diagnostics reduce the corrections substantially:

```text
Zn: C6_PSP+core 285.75 -> 264.69
Cd: C6_PSP+core 508.39 -> 364.87
```

The current classification is `review` for both Zn2 and Cd2. Cd raw finite-basis
is close to the Group XII references, but TRK-normalized Cd is not:

```text
Cd raw finite-basis:       +3.1% vs 493, +9.1% vs 466
Cd TRK-normalized:        -26.0% vs 493, -21.7% vs 466
```

Cd all-electron closure is handled through literature proxy targets rather than
the local `ano-rcc` TDHF result. The local `ano-rcc` TDHF calculation gives
`alpha0 = 86.30` and `C6 = 811.86`, which is too high compared with the standard
Cd polarizability scale and the high-reference `686` C6. It is retained only as
a diagnostic control. The Cd literature closure proxy table is
`results/semicore_zn_cd_validation/cd_literature_closure_proxy.csv`.

A first radial-grid Sternheimer prototype now exists in
`compute_core_sternheimer_radial_grid.py`. It solves the imaginary-frequency
response directly in a radial box, so it is no longer a Gaussian virtual-MO sum.
The current `ano-rcc` smoke diagnostics are:

```text
Zn 3d10: cphf_alpha0_core = 3.6082, cphf_C6_self_core = 7.1652
         signed_sum_osc / N_core = 0.9921
         raw_positive_sum_osc / N_core = 1.5040
         cphf_positive_sum_osc / N_core = 1.0000
         operator_trk_pass = true, psp_go_no_go_ready = true
Cd 4d10: cphf_alpha0_core = 4.6156, cphf_C6_self_core = 18.1363
         signed_sum_osc / N_core = 0.9796
         raw_positive_sum_osc / N_core = 1.5486
         cphf_positive_sum_osc / N_core = 1.0000
         operator_trk_pass = true, psp_go_no_go_ready = true
```

This is now folded into the PSP C6 workflow as a Pauli-positive radial continuum
channel set with a shell TRK constrained CPHF/local-field scale. It passes the
oscillator gate, but the resulting C6 values are still `review`, not `go`:

```text
Zn2: C6_PSP = 244.46 -> C6_PSP+radial_CPHF = 320.62
     error vs 257.5 measurement = +24.5%

Cd2: C6_PSP = 271.09 -> C6_PSP+radial_CPHF = 397.21
     error vs 493 present calculation = -19.4%
```

The continuum response therefore fixes the oscillator-strength bookkeeping
problem but not the final Zn/Cd physics target. The next method step is to
improve the physical CPHF kernel/central field and validate it against dynamic
polarizability or better all-electron response data.

A static-alpha constrained local-field kernel is now implemented in
`run_radial_cphf_alpha_constrained_validation.py`. It uses the radial continuum
CPHF channel spectrum, then constrains the observable local-field response to
published static polarizabilities from
`reference_static_polarizability_group12.csv` rather than fitting C6 directly:

```text
Zn alpha0 reference: 37.95 +/- 0.77 a.u.
Cd alpha0 reference: 45.68 +/- 1.21 a.u.
```

The resulting kernel scales and C6 values are:

```text
Zn:
  alpha0_PSP = 37.7369
  alpha0_core_unconstrained = 3.6082
  local_field_scale = 0.0591
  C6_PSP = 244.46 -> C6_PSP+alpha_constrained_CPHF = 248.56
  error vs 257.5 measurement = -3.47%

Cd:
  alpha0_PSP = 37.4355
  alpha0_core_unconstrained = 4.6156
  local_field_scale = 1.7862
  C6_PSP = 271.09 -> C6_PSP+alpha_constrained_CPHF = 521.84
  error vs 493 present calculation = +5.85%
```

This is the first Zn/Cd radial-continuum result that enters the primary 10%
window. It should still be reported with reference sensitivity: Zn is `+10.47%`
versus the `225` Group XII present-calculation value, Cd is `+11.98%` versus the
`466` measurement, and both remain far below the high `359/686` QDO/CCSD(T)
stress references.

Sr is now a warning/control case for finite-basis core-response basis
sensitivity. The original def2 core response gave an unphysical correction
(`2733 -> 9268`, about `+192%` vs `3170`), while `ano-rcc` core channels give a
reasonable sanity check (`2733 -> 3120`, about `-1.6%`; TRK-normalized
`2733 -> 2977`, about `-6.1%`). This indicates that the Sr shell partition is
clean, but finite-basis core-response basis choice can dominate the correction.

## 7. Near-Term Roadmap

Next steps:

1. Decide how to present the second benchmark:
   - headline remains Mg q2
   - Sn q4 is the strongest official-matched secondary candidate, with a
     reference-method caveat
   - Ge q4 is a weaker official-matched control
2. Import or construct a matched large-core q2 basis route for Ca or Sr:
   - target Ca q2: `GTH-*-q2` with explicit `4s`
   - target Sr q2: `GTH-*-q2` with explicit `5s`
   - require a PySCF build/RKS/TDDFT smoke test before benchmark use
3. Re-run the clean benchmark chain for the first successful Ca/Sr route:
   - PSP valence-only TDDFT
   - PSP + `l=1` dipole EFT semicore correction
   - all-electron TDDFT reference with the same XC backend where possible
4. Clean up documentation and make all current benchmarks reproducible.
5. Add plots for:
   - noble gas TDHF benchmark
   - semicore/core correction summary
   - Ca basis sensitivity
6. Add Ca2 long-range tail comparison:
   - all-electron TDHF
   - valence-only
   - valence+semicore
7. Begin screened-pairwise prototype with a model `W_v`:
   - bare
   - dielectric
   - Yukawa/Thomas-Fermi
8. Later: implement finite-system log/MBD energy.
9. Long term: connect to ab initio valence response `W_v` and periodic systems.

## 8. Current Status In One Sentence

This repository currently provides a working atomic response-to-C6 prototype
with a clean Mg q2 PSP -> PSP+dipole-EFT -> all-electron benchmark loop, plus a
strong official-matched Sn q4 secondary candidate and a model-screened
pairwise/logdet interface. The final ab initio screened EFT-vdW functional,
periodic implementation, forces, and production-quality alkaline-earth
large-core PSP benchmark remain future work.
