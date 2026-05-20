# Import Large-Core q2 Basis Task

Created: 2026-05-21

## Goal

Find and integrate one clean large-core q2 PSP+basis route for a second
semicore benchmark beyond Mg q2.

Preferred targets:

```text
Ca q2: explicit PSP valence = 4s; EFT shells = 3s,3p
Sr q2: explicit PSP valence = 5s; EFT shells = 4s,4p
```

## Motivation

Mg q2 is currently the only strong clean closure case:

```text
C6_PSP -> C6_PSP+EFT -> C6_all-e
638.6202 -> 647.6079 -> 647.5881
```

Clean deep-core controls do not close the gap:

```text
Be q2 LDA  closure = 2.26%
Kr q8 PBE  closure = 4.28%
Ca q10 PBE closure = 0.081%
```

This means the second benchmark should target missing semicore response, not
only missing deep-core response.

## Current Local Data Status

The current local CP2K scan found q2 pseudos for Ca/Sr/Ba/Zn/Cd/Hg, but no
matched q2 basis for these target atoms.

Known positive control:

```text
Mg GTH-PBE-q2 / TZV2P-MOLOPT-SR-GTH-q2 -> TDDFT smoke OK
```

Current blocker:

```text
Ca/Sr q2 pseudo exists, but no local matched Ca/Sr q2 basis has been found.
```

An additional official CP2K master-data check found Ca/Sr `GTH-PBE-q2` blocks
in `POTENTIAL_UZH`, but still did not find Ca/Sr q2 basis blocks in the checked
CP2K basis files.

The extracted Ca/Sr PBE q2 pseudo blocks are stored in:

```text
external_data/cp2k/POTENTIAL_UZH_CASR_Q2
```

Source:

```text
https://raw.githubusercontent.com/cp2k/cp2k/master/data/POTENTIAL_UZH
```

## Adapted Basis Smoke Test

Because no library-native Ca/Sr q2 basis was found, the scanner now supports
explicitly listed imported/adapted candidate bases with provenance metadata.

Candidate list:

```text
external_data/cp2k/large_core_q2_basis_candidates.csv
```

Current adapted candidates use UZH q10 MOLOPT basis blocks as q2 candidates.
This is not a library-native matched q2 basis and must remain labeled
`adapted_from_q10`.

Smoke-test outputs:

```text
results/import_large_core_q2_basis_scan.csv
results/import_large_core_q2_basis_pbe_scan.csv
```

Passing PBE-consistent smoke tests:

```text
Ca GTH-PBE-q2 / TZV2P-MOLOPT-PBE-GTH-q10-as-q2-adapted -> TDDFT smoke OK
Sr GTH-PBE-q2 / TZV2P-MOLOPT-PBE-GTH-q10-as-q2-adapted -> TDDFT smoke OK
Sr GTH-PBE-q2 / QZVPP-MOLOPT-PBE-GTH-q10-as-q2-adapted -> TDDFT smoke OK
```

All passing rows are diagnostic candidates only because the basis provenance is
`adapted_from_q10`.

## First Adapted Benchmark Diagnostic

The first post-smoke benchmark was run for Ca:

```text
case:        Ca_q2_PBE_adapted
pseudo:      GTH-PBE-q2 from POTENTIAL_UZH_CASR_Q2
basis:       TZV2P-MOLOPT-PBE-GTH-q10-as-q2-adapted
explicit:    4s
EFT shells:  3s,3p
backend:     PBE-TDDFT
audit:       pass
```

Result:

```text
C6_PSP        = 1496.31224087
C6_PSP+EFT    = 1658.03166276
C6_all-e_PBE  = 2206.75882541
closure       = 22.7631%
residual_C6   = 548.72716265
```

Interpretation:

This is a useful large-core semicore diagnostic and is much more relevant than
the deep-core controls, but it is not yet the production second benchmark
because the basis is adapted from q10 rather than a native or optimized q2
basis.

## Acceptance Criteria

A candidate Ca/Sr q2 route is acceptable only if:

1. The basis provenance is recorded:
   - library-native
   - imported external
   - constructed/adapted
2. PySCF can build the atom with the q2 pseudo and candidate basis.
3. RKS converges.
4. TDDFT smoke test runs.
5. The PSP summary reports the intended active shell:
   - Ca q2: `4s`
   - Sr q2: `5s`
6. The EFT correction shells do not overlap the explicit PSP shells.

## Implementation Plan

1. Extend the q2 scanner to accept imported candidate basis files.
2. Add candidate basis files under `external_data/` with provenance notes.
3. Run smoke tests:

```text
probe_large_core_q2_candidates.py
```

4. For the first passing Ca/Sr route, run:

```text
PSP-RPA valence-only TDDFT
all-electron TDDFT reference
PSP + l=1 dipole EFT semicore correction
benchmark audit
```

5. Write the result to a dedicated summary file rather than mixing it with
deep-core controls.

Current implementation status:

```text
scanner imported/adapted basis support: done
Ca/Sr adapted q2 smoke tests: done
Ca q2 PBE adapted diagnostic benchmark: done
native/optimized Ca/Sr q2 basis: still open
```

## Non-Goals

- Do not treat Be q2, Kr q8, or Ca q10 as the second strong benchmark.
- Do not call constructed/adapted bases library-native.
- Do not add shells already explicit in the PSP valence space.
- Do not claim a production Ca/Sr benchmark until the q2 basis provenance and
  smoke tests are documented.
