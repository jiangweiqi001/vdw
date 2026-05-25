# Mg q2 Clean PSP + EFT-Core Benchmark

This is the current strongest clean benchmark for the project.

It demonstrates the central PSP+EFT-core idea at the atomic C6 level:

```text
PSP valence-only response
  + frozen semicore l=1 dipole response
  -> nearly all-electron C6
```

It is not yet a final screened periodic EFT-vdW functional.

## 1. Benchmark Definition

The benchmark compares three responses:

```text
all-electron PBE-TDDFT reference
q2 PSP-RPA valence-only baseline
q2 PSP-RPA + l=1 EFT-core dipole correction
```

The PSP leg is:

```text
atom: Mg
pseudopotential: GTH-PBE-q2
basis: TZV2P-MOLOPT-SR-GTH-q2
explicit PSP valence shell: 3s
```

The EFT correction targets the frozen semicore shells:

```text
EFT correction shells: 2s,2p
```

## 2. Double-Counting Audit

The shell overlap is empty:

```text
explicit_valence_shells ∩ eft_core_shells = empty
3s ∩ (2s,2p) = empty
```

The benchmark row is therefore clean:

```text
double_counting_status = clean
```

This is the main reason Mg q2 is stronger than the Ca q2 adapted diagnostic.

## 3. C6 Closure

Main result:

```text
C6_PSP                  = 638.62015545
C6_PSP + EFT_dipole     = 647.69960692
C6_all-electron PBE     = 647.58810041

Delta_C6_missing        = 8.96794495
Delta_C6_EFT            = 9.07945147
residual_C6             = -0.11150651
closure                 = 101.24%
```

The additive l=1 dipole correction essentially closes the PSP-to-all-electron
C6 gap for Mg q2.

## 4. Core-Response Variants

Two related implementations bracket the core correction:

```text
neutral-atom MO dipole approximation:
  C6_PSP + EFT = 647.69960692
  Delta_C6     = +9.07945147

core-ion TDHF transition-density proxy:
  C6_PSP + EFT = 643.73374504
  Delta_C6     = +5.11358959
```

The first is the current headline benchmark because it closes the all-electron
PBE C6 gap. The second is closer to a core-sector transition-density response
and gives a smaller correction. This difference is an estimate of the present
uncertainty in the l=1 dipole Wilson approximation.

The explicit multipole transition-density implementation is:

```text
compute_multipole_core_wilson.py
source = EFT_CORE_MULTIPOLE_TDENSITY_TDHF
```

The small-q form factor check is stored in:

```text
results/multipole_form_factors_mg_core.csv
```

It explicitly evaluates, along the z axis,

```text
tau_lambda(q) = i q d_lambda,z + O(q^3)
F_lambda(q) = 4 pi tau_lambda(q) / q^2
```

and verifies that:

```text
dipole_from_tau_z = tau_lambda(q) / (i q)
```

is independent of the small q value. This is the explicit transition-density /
small-q Wilson check corresponding to the derivation in
`docs/eft_core_alpha_derivation.md`.

## 5. Mg2 Long-Range Tail

The long-range tail comparison is stored in:

```text
results/mg/mg2_tail_comparison.csv
```

It compares:

```text
E_all-e(R)       = -C6_all-e / R^6
E_PSP(R)         = -C6_PSP / R^6
E_PSP+EFT(R)     = -C6_PSP+EFT / R^6
```

Since this is a pure C6 tail, the percentage error is independent of R:

```text
PSP error      = -1.3848%
PSP+EFT error  = +0.0172%
```

Thus the PSP-only tail is slightly too weak, while PSP+EFT essentially recovers
the all-electron PBE tail.

## 6. Model Screening Sensitivity

Model `W_v` sensitivity is stored in:

```text
results/mg/mg_q2_model_screening_sensitivity.csv
```

It includes:

```text
dielectric: epsilon = 1, 2, 4
yukawa:     kappa = 0.05, 0.1, 0.2
```

These are model screenings only:

```text
model W_v sensitivity; not ab initio screening
```

They should not be interpreted as a completed ab initio screened EFT-vdW
functional.

## 7. Screened / Logdet Interface Check

The same Mg q2 PSP+EFT channels were passed through the finite-system
second-order logdet interface:

```text
results/mg_q2_eft_logdet_bare_R20.csv
results/mg_q2_eft_logdet_dielectric2_R20.csv
```

For the bare model, the second-order logdet result recovers `-C6/R^6`. For
`epsilon=2`, the energy is reduced by `1/epsilon^2`, as expected for a model
screened dipole tensor.

## 8. Limitations

This benchmark is clean and numerically strong, but it is still not the final
screened EFT-vdW functional:

- The l=1 dipole correction is currently an MO/core-ion transition-density
  approximation.
- The best closure uses the neutral-atom MO dipole approximation, not yet a
  closed-form analytic core-only Wilson coefficient.
- The correction is additive and unscreened at the atomic C6 level.
- Ab initio `W_v = (v^-1 - chi_v)^-1` is not yet implemented.
- Periodic systems, forces, and production log-determinant calculations are not
  implemented.

## 9. Reproduction

The relevant pipeline is:

```bash
python3 run_psp_rpa_atom.py \
  --atom Mg \
  --psp GTH-PBE-q2 \
  --basis placeholder \
  --xc pbe \
  --nstates 100 \
  --method TDDFT \
  --basis-file external_data/cp2k/BASIS_MOLOPT_UCL \
  --basis-name TZV2P-MOLOPT-SR-GTH-q2 \
  --pseudo-file external_data/cp2k/GTH_POTENTIALS \
  --pseudo-name GTH-PBE-q2

python3 compute_dipole_wilson.py \
  --atom Mg \
  --basis aug-cc-pVQZ \
  --shell 2s \
  --shell 2p \
  --output results/eft_core_dipole_wilson_channels.csv

python3 run_eft_core_dipole_validation.py
python3 run_mg_q2_model_screening_sensitivity.py
```

The clean benchmark rows are recorded in:

```text
results/mg_q2/summary.csv
results/eft_core_dipole_validation_summary.csv
results/mg/mg2_tail_comparison.csv
results/mg/mg_q2_model_screening_sensitivity.csv
```
