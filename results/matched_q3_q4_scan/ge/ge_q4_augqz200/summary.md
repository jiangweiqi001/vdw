# Ge q4 Improved Reference Closure

## Setup

```text
case = Ge_q4_augqz200_reference
PSP = GTH-PBE-q4
PSP basis = TZV2P-MOLOPT-PBE-GTH-q4
PSP nstates = 200
explicit PSP shells = 4s, 4p
EFT correction shell = 3d
all-electron reference = PBE TDDFT, aug-cc-pVQZ, nstates=200
```

This uses an official matched Ge q4 PSP/basis pair. The all-electron reference
is the best currently available Ge reference from the convergence scan, but it
is still not fully converged with respect to `nstates`.

## Result

| quantity | C6 au |
|---|---:|
| PSP-RPA | 306.65944141 |
| PSP-RPA + EFT 3d | 317.45391425 |
| all-electron PBE TDDFT | 375.57206997 |

```text
missing C6 gap = 68.91262856
EFT 3d adds    = 10.79447284
residual C6    = 58.11815572
closure        = 15.66399812%
```

## Interpretation

Using a more reasonable diffuse all-electron reference changes Ge q4 from a
negative-closure smoke result into a positive but modest correction.

The frozen `3d` EFT channel closes about 15.7% of the PSP/all-electron C6 gap.
That is useful as an official-matched control, but it is much weaker than Mg q2
and should not replace Mg q2 as the headline benchmark.

Current label:

```text
Ge q4 = official-matched control / secondary candidate, not strong benchmark
```
