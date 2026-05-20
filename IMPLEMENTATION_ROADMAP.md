# EFT-vdW Implementation Roadmap

This document reframes the project around the implementation route in issue #1:

```text
all-electron RPA reference
PSP-valence RPA baseline
EFT core/semicore correction
validation against reference C6 and prior-art decompositions
```

The existing TDHF/MO/radial work remains useful, but it should be treated as
diagnostic infrastructure and reference-data generation, not yet as the final
PSP+EFT-core method.

## 0. Current Baseline And Reframing

Current code already provides:

- CSV oscillator-channel backend: `delta_Ha, osc -> alpha(i xi) -> C6`
- calibrated controls
- all-electron PySCF TDHF oscillator export
- all-electron TDHF transition-dipole projection diagnostics
- noble-gas and Ca/Mg semicore diagnostic tables
- test coverage for the response pipeline

Current code does not yet provide:

- a pseudopotential valence-only RPA/TDDFT response
- an EFT-derived `alpha_core(i xi)` from PRL Wilson coefficients `f_K^c`
- a screened `W_v`
- a finite-system or periodic log/RPA vdW energy

Important terminology:

- **All-electron TDHF/RPA**: reference or diagnostic.
- **Projected all-electron TDHF partition**: diagnostic only.
- **PSP-RPA**: the intended cheap valence-only leg.
- **EFT-core correction**: the future Wilson-coefficient correction from
  `f_K^c` and `Delta E_c`; not yet implemented.

## 1. Milestone A: All-Electron RPA Reference

Goal: establish a consistent all-electron response reference for atoms.

### A1. Add KS-LDA / KS-PBE all-electron response

Current implementation uses spherical all-electron HF + TDHF. The EFT matching
language is closer to Kohn-Sham LDA/PBE, so add a KS-based reference leg.

Deliverable:

```text
run_all_e_rpa_atom.py
results/all_e_rpa_summary.csv
```

Suggested CLI:

```bash
python3 run_all_e_rpa_atom.py --atom Ar --xc lda --basis aug-cc-pVQZ --nstates 200
python3 run_all_e_rpa_atom.py --atom Ca --xc lda --basis cc-pVQZ --nstates 200
```

Required output columns:

```text
atom, xc, basis, nstates, alpha0, C6, alpha0_ref, alpha0_err, C6_ref, C6_err
```

Acceptance criteria:

- Reproduces existing HF/TDHF order of magnitude.
- Provides stable `alpha0` and `C6` for Ne/Ar/Kr/Mg/Ca.
- Records whether the response is TDHF, TDA, TDDFT, or RPA.

### A2. Keep existing HF/TDHF as secondary check

Existing files:

```text
pyscf_export_ar_tdhf_oscillators.py
run_tdhf_atom.py
results/noble_gas_tdhf_summary.csv
```

Keep them, but label them as all-electron TDHF references/diagnostics.

## 2. Milestone B: PSP-Valence RPA Baseline

Goal: compute the response that a frozen-core pseudopotential calculation would
actually see.

This is the biggest missing piece.

### B1. Pseudo-atom setup

Use a pseudopotential/ECP calculation for the atom with only valence electrons
active. Candidate routes:

- PySCF GTH/HGH pseudo atoms if available.
- DFTK/Julia for GTH pseudopotentials.
- A small custom pseudo-atom wrapper if PySCF basis support is limiting.

Deliverable:

```text
run_psp_rpa_atom.py
results/psp_rpa_summary.csv
```

Suggested CLI:

```bash
python3 run_psp_rpa_atom.py --atom Ca --psp gth-lda --basis gth-dzvp --zval large-core
python3 run_psp_rpa_atom.py --atom Ca --psp gth-lda --basis gth-dzvp --zval small-core
```

Required output columns:

```text
atom, psp, basis, zval, active_shells, alpha0_psp, C6_psp
```

Acceptance criteria:

- Produces a genuine PSP valence-only response, not an all-electron projected
  response.
- Can run at least for Ar/Ne/Kr if suitable PSPs exist, and for Mg/Ca with
  large-core/small-core partition if possible.
- Stores enough metadata to know which shells are explicit valence.

### B2. Compare all-electron and PSP response

Deliverable:

```text
results/all_e_vs_psp_rpa_summary.csv
```

Required output columns:

```text
atom, partition, C6_all_e, C6_psp, delta_C6_missing, relative_missing_pct
```

Acceptance criteria:

- Quantifies what frozen-core PSP actually misses.
- Separates "large-core" from "small-core" cases for Mg/Ca and later K/Sr.

## 3. Milestone C: EFT Core/Semicore Correction Derivation

Goal: derive `Delta alpha_core^EFT(i xi)` from the same Wilson coefficients used
in the PRL bandwidth correction.

This is primarily a theory milestone before coding.

### C1. Starting point

The PRL dynamic core term is:

```text
delta V_pp^dyn(K,K'; i omega)
    = sum_c f_K^c f_K'^c / (i omega + Delta E_c)
```

For bandwidth narrowing, the relevant quantity is the derivative near the Fermi
energy, yielding `z_core`.

For vdW, the target is a correction along the imaginary axis:

```text
Delta alpha_core^EFT(i xi)
```

### C2. Derivation questions to settle

1. Which static piece is already absorbed into the pseudopotential?
2. How does the one-particle dynamic pole map to a density response kernel?
3. What is the exact long-wavelength dipole limit of `f_K^c`?
4. Does the finite-`xi` matching use the same coefficient as the `xi -> 0`
   bandwidth derivative?
5. Is the core-valence cross term already included in PSP-RPA, or does it need a
   separate Wilson term?

Deliverable:

```text
docs/eft_core_alpha_derivation.md
```

Acceptance criteria:

- States the double-counting subtraction.
- Gives a formula for `Delta alpha_core^EFT(i xi)`.
- States the expected small parameter and validity domain.
- Connects the static/low-frequency limit to `z_core`.

## 4. Milestone D: Wilson Coefficient Evaluation

Goal: compute the PRL Wilson coefficients and excitation energies needed for
the EFT correction.

### D1. Atomic core quantities

Inputs:

```text
u_c(r), V_H,c(r), J_c, Delta E_c
```

Core form factor:

```text
f_K^c = sqrt(4 pi) / K * integral u_c(r) [V_H,c(r) - J_c] sin(K r) dr
```

Deliverables:

```text
compute_core_wilson.py
core_wilson_coefficients.csv
```

Suggested columns:

```text
atom, shell, Delta_E_Ha, Jc_Ha, K_grid, fK_values, source
```

Acceptance criteria:

- Reproduces the qualitative PRL trends for Na/K/Ca/Mg.
- Can compute at least the dominant semicore shell for Mg and Ca.
- Includes tests for normalization and limiting behavior.

### D2. Convert Wilson coefficients to EFT alpha correction

Deliverable:

```text
build_eft_core_alpha_channels.py
eft_core_channels.csv
```

Suggested columns:

```text
atom, shell, channel, delta_Ha, osc_or_weight, source
```

Acceptance criteria:

- Provides a backend-compatible oscillator/channel representation, or a clearly
  documented direct `alpha(i xi)` correction.
- Does not reuse all-electron TDHF transition strengths as the EFT correction.

## 5. Milestone E: PSP + EFT-Core Validation

Goal: test whether the EFT correction closes the gap between PSP-RPA and
all-electron/reference response.

Compute:

```text
alpha_EFT(i xi) = alpha_PSP-RPA(i xi) + Delta alpha_core^EFT(i xi)
C6_EFT = Casimir-Polder[alpha_EFT]
```

Deliverable:

```text
results/eft_core_validation_summary.csv
```

Required columns:

```text
atom, partition, C6_ref, C6_all_e, C6_psp, C6_eft,
delta_psp, delta_eft, recovered_fraction, note
```

Acceptance criteria:

1. `C6_EFT - C6_PSP` has the right sign and order of magnitude.
2. `C6_EFT` approaches `C6_all_e` for noble gases.
3. For Mg/Ca large-core partitions, the EFT correction captures a meaningful
   fraction of the missing semicore response.

## 6. Milestone F: Long-Range Tail And Method Comparison

Goal: translate C6 differences into long-range vdW tails.

Deliverables:

```text
results/<atom>/<atom>2_tail_comparison.csv
results/method_comparison_summary.csv
```

Compare:

- reference
- all-electron RPA/TDHF
- PSP-RPA
- PSP + EFT-core
- calibrated control
- optional D4/MBD if available

Acceptance criteria:

- Only compare long-range `-C6/R^6` tails first.
- Do not claim full binding curves until damping, exchange repulsion, and
  short-range matching are handled.

## 7. Milestone G: Screened Pairwise Prototype

Goal: replace bare Coulomb dipole tensor with a model screened interaction.

Candidate models:

- bare Coulomb
- dielectric constant screening
- Yukawa / Thomas-Fermi screening

Deliverable:

```text
screened_pairwise_vdw.py
```

Acceptance criteria:

- Demonstrates how `W_v` modifies the long-range pairwise interaction.
- Keeps model screening clearly separate from ab initio `W_v`.

## 8. Milestone H: Full Screened EFT-vdW

Long-term target:

```text
W_v(i xi) = [v^-1 - chi_v^irr(i xi)]^-1

E_vdW ~ (1 / 2pi) int dxi Tr ln[1 - W_v chi_c]
```

This requires:

- valence response in finite or periodic systems
- screened interaction in real or reciprocal space
- double-counting projector for local/static pieces
- periodic implementation and eventually forces

This is not part of the immediate implementation.

## 9. Immediate Next Tasks

Recommended order:

1. Add `run_all_e_rpa_atom.py` for KS-LDA/PBE all-electron response.
2. Add `run_psp_rpa_atom.py` for genuine pseudopotential valence response.
3. Write `docs/eft_core_alpha_derivation.md`.
4. Implement `compute_core_wilson.py` for PRL `f_K^c`.
5. Validate `PSP + EFT-core` against all-electron/reference C6.

Do not use the current all-electron TDHF partition numbers as the EFT-core
correction itself. They are diagnostics and motivation for the PSP+EFT-core
implementation.
