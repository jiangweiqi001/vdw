# EFT-vdW Prototype

This repository prototypes an EFT-style atomic channel pipeline for computing
imaginary-frequency polarizabilities, C6 coefficients, and long-range pairwise
dispersion tails.

## Setup

Install the Python dependencies and run the test suite:

```bash
python3 -m pip install -r requirements.txt
python3 -m unittest discover
```

## Minimal Ar Reproduction

Reproduce the current Ar TDHF benchmark and Ar2 long-range tail:

```bash
python3 pyscf_export_ar_tdhf_oscillators.py --basis aug-cc-pVQZ --nstates 200 --output ar_tdhf_channels.csv
python3 run_alpha_table.py --input ar_tdhf_channels.csv --output ar_tdhf_alpha_c6_table.csv
python3 run_c6_table.py --input ar_tdhf_channels.csv --output ar_tdhf_c6_table.csv
python3 compare_alpha_c6.py --eft ar_tdhf_alpha_c6_table.csv
python3 run_ar2_tail.py
```

## Current Ar Prediction Baseline

The PySCF 3D MO oscillator route is the current prediction baseline for Ar.
Its summary is stored in `results/ar/ar_mo_basis_convergence.csv`.

The radial shell-average route is diagnostic only until the per-shell TRK
overcounting is fixed. In particular, diffuse radial shell-average runs can
violate the oscillator-strength sum rule even when the corresponding 3D MO
dipole-integral route does not.

## Ar Basis Convergence Caveat

`results/ar/ar_basis_convergence.csv` includes three standard PySCF basis sets:

- `cc-pVTZ`
- `aug-cc-pVTZ`
- `aug-cc-pVQZ`

It also includes `d-aug-cc-pVTZ-local`. This is not an official Basis Set
Exchange `d-aug-cc-pVTZ` basis. It is a locally generated, even-tempered second
diffuse augmentation built from PySCF's `aug-cc-pVTZ` basis by adding one extra
diffuse primitive per angular momentum channel. Treat it as a reproducible local
sensitivity test, not as a published basis-set benchmark.

## Status And Limitations

Implemented:

- Oscillator-channel `alpha(i xi)` and C6 backend
- Calibrated, MO, TDHF, and radial diagnostic routes
- Ar TDHF benchmark and Ar2 long-range tail
- Ar basis, TDHF nstates, and oscillator-strength audit tables

Not yet implemented:

- Valence/core decomposition
- Valence-screened `W_v`
- RPA/log-determinant many-body energy
- Periodic systems and forces
