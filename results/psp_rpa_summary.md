# PSP-Valence RPA Baseline: Current Status

## What Was Implemented

Milestone B has been started with a PySCF GTH pseudopotential response leg:

```text
run_psp_rpa_atom.py
run_all_e_vs_psp_rpa_summary.py
```

The implementation runs a pseudo-atom `RKS` calculation with a GTH
pseudopotential and GTH basis, then PySCF `TDDFT`, exports oscillator-strength
channels, and reuses the existing `alpha(i xi) -> C6` backend.

## PySCF PSP Basis Availability

`probe_pyscf_psp_basis.py` scans the available PySCF GTH pseudopotential/basis
combinations. The resulting table is:

```text
results/psp_basis_availability.csv
```

Main findings:

```text
Ne, Ar, Mg: standard gth-qzv2p is available.
Kr, Ca: standard gth-* basis files are mostly unavailable, but
        gth-dzvp-molopt-sr is available.
gth-hf pseudopotentials are unavailable for the current test atoms.
Additional CP2K data were imported from the official CP2K repository:

```text
external_data/cp2k/BASIS_MOLOPT_UZH
external_data/cp2k/GTH_POTENTIALS
```

These provide larger PBE GTH basis/pseudopotential combinations for Kr and Ca:

```text
Kr: TZV2P-MOLOPT-PBE-GTH-q8 + GTH-PBE-q8
Ca: TZV2P-MOLOPT-PBE-GTH-q10 + GTH-PBE-q10
```

Julia/DFTK is not available in this environment, so DFTK has not been tested.
```

## First PSP-RPA Smoke Results

The first smoke run used:

```text
psp = gth-lda
basis = gth-dzvp
xc = lda
method = TDDFT
nstates = 50
```

Available results:

```text
atom  active_electrons  alpha0_psp  C6_psp
Ne    8                 0.5954      0.7804
Ar    8                 4.1363      16.8213
Mg    10                55.8837     441.2003
```

Unavailable in the current PySCF GTH-DZVP library:

```text
Kr, Ca
```

## Larger Available PSP-RPA Cases

Using the largest available PySCF GTH basis for each atom in this environment:

```text
atom  psp/basis                         active_electrons  alpha0_psp  C6_psp
Ne    gth-lda / gth-qzv2p               8                 1.5158      2.6915
Ar    gth-lda / gth-qzv2p               8                 8.0405      41.3880
Mg    gth-lda / gth-qzv2p               10                69.5109     603.2473
Kr    GTH-PBE-q8 / TZV2P-MOLOPT-PBE     8                 9.9190      66.3555
Ca    GTH-PBE-q10 / TZV2P-MOLOPT-PBE    10                119.5029    1424.6731
```

## Clean Large-Core Mg q2 Case

A clean Mg large-core PSP diagnostic is now present:

```text
atom  psp/basis                              active_shells  alpha0_psp  C6_psp
Mg    GTH-PBE-q2 / TZV2P-MOLOPT-SR-GTH-q2    3s             72.1168     638.6202
```

This is the first clean benchmark candidate because the PSP explicit valence
space is `3s`, while the EFT dipole correction shells are `2s,2p`. The current
validation chain is:

```text
C6_PSP              = 638.6202
C6_PSP+dipole_EFT   = 647.6079
C6_all-e_PBE_TDDFT  = 647.5881
double counting     = clean
```

The Mg q2 result was generated from CP2K GTH q2 pseudo/basis data already stored
under `external_data/cp2k`. It has been freshly reproduced with WSL Ubuntu
PySCF 2.12.1, and is now wrapped by `run_mg_q2_benchmark.py`.

## All-Electron vs PSP-RPA Gap

Using the matching all-electron TDDFT reference from
`results/all_e_rpa_summary.csv`:

```text
atom  best current PSP basis         C6_all_e    C6_psp    missing fraction
Ne    gth-pbe/gth-qzv2p              7.0720      2.7319    61.4%
Ar    gth-pbe/gth-qzv2p              68.2301     41.8306   38.7%
Mg    GTH-PBE-q2/TZV2P-MOLOPT-SR     647.5881    638.6202  1.4%
Kr    TZV2P-MOLOPT-PBE-GTH-q8        134.6496    66.3555   50.7%
Ca    TZV2P-MOLOPT-PBE-GTH-q10       2206.7588   1424.6731 35.4%
```

These PSP-RPA values are intentionally crude first baselines. They show that a
genuine frozen-core pseudopotential response can miss a large amount of
all-electron oscillator strength, especially for the noble-gas-like GTH valence
spaces. The Mg q2 row is now the clean large-core benchmark candidate rather
than a small-core diagnostic.

## Interpretation

This is now a real PSP-valence response calculation, unlike the earlier
all-electron TDHF projection diagnostics.

This milestone is complete as a prototype baseline, but not as a final
quantitative PSP-RPA benchmark:

- The original GTH-DZVP basis is too small for quantitative C6.
- Larger GTH-QZV2P improves Ne/Ar/Mg substantially and is the current best
  PySCF-internal choice for these atoms.
- Kr and Ca require external CP2K basis files to access the larger UZH
  `TZV2P-MOLOPT-PBE` GTH basis; these are now wired into PySCF and working.
- The original Mg PBE GTH-QZV2P row is a 10-electron small-core setup, but a
  clean 2-electron Mg q2 row is now available.
- Ca has both a 10-electron GTH-PBE q10 setup and a clean-by-shell-overlap q2
  diagnostic. The q2 Ca row uses `GTH-LDA-q2 + cc-pVQZ`, so it should not yet be
  treated as a matched production-quality GTH basis benchmark.
- Julia/DFTK is not installed in this environment, so the DFTK route has not
  been tested yet.
- Larger GTH bases, alternate PSP libraries, installing DFTK, or importing a
  true large-core Mg/Ca PSP may still be needed for final quantitative PSP-RPA.

## Milestone B Closeout

What is done:

1. A genuine PSP-valence TDDFT response route exists.
2. Ne, Ar, Mg, Kr, and Ca have working PSP-RPA rows.
3. External CP2K basis/pseudopotential data are wired in for Kr and Ca.
4. `results/all_e_vs_psp_rpa_summary.csv` quantifies the PSP missing-C6 gap.

What is not done:

1. Production-quality large-core Ca q2 with a matched PSP basis is not yet
   available.
2. PSP basis convergence is not complete.
3. DFTK/periodic PSP-RPA is not installed or tested.

Current decision:

```text
Milestone B is closed for prototype purposes.
The Mg q2 row is sufficient to start a clean large-core PSP + EFT-core benchmark
loop, but the broader PSP-RPA baseline is still not sufficient for final
production-level quantitative claims.
```

Next milestone:

```text
Milestone C/D: derive and implement the EFT-core oscillator correction from
the PRL Wilson coefficients, then test whether PSP + EFT-core closes the
missing-C6 gap.
```
