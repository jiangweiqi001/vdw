# Ca q2 Generated Basis v1 Validation

## Status

```text
case = Ca-q2-generated-v1
basis label = generated_protocol_frozen
official matched q2 basis = false
vdW targets used in basis generation = false
double counting = clean
```

This validation uses the frozen project-local generated basis:

```text
generated_basis/Ca-q2-generated-v1.bas
sha256 = 0bd7511717e0714e095d0938838c7602cde9ac5c0201875cfdac186397a8dacc
freeze record = generated_basis/Ca-q2-generated-v1.freeze.json
```

The numeric basis block was copied unchanged from the official Ca
`TZV2P-MOLOPT-PBE-GTH-q10` seed and renamed to `Ca-q2-generated-v1` before any
vdW validation. This is not an official CP2K/MOLOPT matched q2 basis.

## Validation Setup

```text
PSP                 = GTH-PBE-q2
PSP basis           = Ca-q2-generated-v1
PSP explicit shell  = 4s
EFT correction      = 3s, 3p l=1 MO dipole Wilson approximation
all-electron ref    = PBE TDDFT, cc-pVQZ, nstates=200
```

The shell partition is clean:

```text
explicit PSP valence shells = 4s
EFT correction shells       = 3s, 3p
overlap                     = none
```

## C6 Closure

| quantity | C6 au |
|---|---:|
| PSP-RPA | 1496.31224087 |
| PSP-RPA + EFT-core | 1658.03166276 |
| all-electron PBE TDDFT | 2206.75882541 |

```text
missing C6 gap = 710.44658454
EFT correction = 161.71942189
residual C6    = 548.72716265
closure        = 22.76306557%
```

## Output Files

```text
results/generated_basis_validation/Ca-q2-generated-v1/summary.csv
results/generated_basis_validation/Ca-q2-generated-v1/psp_alpha_c6_table.csv
results/generated_basis_validation/Ca-q2-generated-v1/psp_plus_eft_alpha_c6_table.csv
results/generated_basis_validation/Ca-q2-generated-v1/ca_3s3p_dipole_wilson_channels.csv
results/generated_basis_validation/Ca-q2-generated-v1/ca_psp_plus_eft_channels.csv
```

## Interpretation

This is a protocol-frozen generated-basis validation, not an official matched
q2 benchmark. The result is numerically identical to the earlier Ca q2 adapted
diagnostic because v1 deliberately uses the same numeric Ca q10 seed block.

The important change is procedural: the basis is now frozen, hashed, and
audited before vdW validation. Within that protocol, the additive EFT
semicore correction closes about 22.8% of the PSP/all-electron C6 gap.

This is a useful Ca semicore diagnostic, but Mg q2 remains the strongest clean
official-matched benchmark.
