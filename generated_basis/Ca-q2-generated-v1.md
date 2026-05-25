# Ca q2 Generated Basis v1

## Status

```text
benchmark_label = generated_protocol_frozen
vdW targets used = false
official matched basis = false
```

This is a project-local generated basis for testing the Ca `GTH-PBE-q2`
pseudopotential path. It is not an official CP2K/MOLOPT matched q2 basis.

## Generation Procedure

The v1 procedure is fixed and non-vdW:

```text
source basis file = external_data/cp2k/BASIS_MOLOPT_UZH
source basis name = TZV2P-MOLOPT-PBE-GTH-q10
generated name    = Ca-q2-generated-v1
rule              = copy numeric CP2K seed block unchanged and replace header
```

This procedure was selected before running vdW validation. It does not use
`alpha0`, `C6`, PSP+EFT closure, or vdW tail energies as optimization targets.

The purpose of v1 is to establish a reproducible generated-basis path and freeze
record. It should be treated as project-local generated/protocol-frozen, not as
an official matched q2 MOLOPT basis.

## Freeze Record

```text
basis file = generated_basis/Ca-q2-generated-v1.bas
freeze record = generated_basis/Ca-q2-generated-v1.freeze.json
provenance record = generated_basis/Ca-q2-generated-v1.provenance.json
sha256 = 0bd7511717e0714e095d0938838c7602cde9ac5c0201875cfdac186397a8dacc
```

Allowed non-vdW targets:

- PSP SCF convergence
- occupied orbital eigenvalue stability
- basis parse/build reproducibility

Forbidden targets were not used:

- `alpha0`
- `C6`
- PSP+EFT closure
- vdW tail energy

## Non-vdW SCF Audit

Audit file:

```text
generated_basis/Ca-q2-generated-v1.scf_audit.csv
```

Summary:

| case | charge | spin | nelectron | nao | SCF converged | gap Ha | status |
|---|---:|---:|---:|---:|---|---:|---|
| Ca | 0 | 0 | 2 | 23 | true | 0.0895766039 | pass |
| Ca+ | 1 | 1 | 1 | 23 | true | 0.0910916949 | pass |

No TDDFT, `alpha0`, `C6`, or PSP+EFT validation was used in this audit.

## Next Validation Step

Only after this freeze point should the project run:

```text
PSP-RPA -> EFT-core correction -> all-electron comparison -> C6 closure
```

If the basis is edited after looking at those vdW quantities, this v1 validation
must be treated as exploratory and a new generated basis version must be frozen.
