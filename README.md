# EFT-vdW Prototype

This repository prototypes an EFT-style atomic channel pipeline for computing
imaginary-frequency polarizabilities, C6 coefficients, and long-range pairwise
dispersion tails.

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
