# Sn q4 Reference Repair

## Goal

Sn q4 is physically attractive because the official matched PSP keeps `5s,5p`
explicit while freezing the relatively shallow `4d` shell:

```text
PSP explicit shells = 5s, 5p
EFT shell           = 4d
```

The original Sn q4 smoke result was unusable because the all-electron
`def2-TZVP/TZVPP` references produced nearly zero C6.

## Diagnosis

The def2 all-electron references are not stable:

| basis | SCF converged | C6 |
|---|---|---:|
| def2-SVP | false | 190.77137368 |
| def2-TZVP | true | 0.19510512 |
| def2-TZVPP | true | 0.19469673 |
| def2-QZVPP | true | 162.20705558 |

This is a reference problem, not evidence that Sn q4 physics is weak.

## ANO Reference Probe

PySCF's built-in `ano` basis gives a much more plausible reference:

```text
Sn all-electron PBE TDDFT
basis   = ano
nstates = 100
alpha0  = 56.62325050
C6      = 569.96394236
```

Partial nstates scan:

| nstates | alpha0 | C6 |
|---:|---:|---:|
| 50 | 52.60337059 | 470.39898966 |
| 60 | 53.33479312 | 487.19834872 |
| 80 | 56.58522693 | 568.86558419 |
| 100 | 56.62325050 | 569.96394236 |
| 120 | 56.62405341 | 569.98822345 |
| 150 | 56.83959940 | 576.98408849 |

The 80, 100, and 120 state results are essentially stable. The 150-state result
is modestly higher, by about 1.2% relative to 120 states, so the full-TDDFT ANO
reference is close to converged for the low-frequency C6 response.

## Near-100 TDA Check

TDA was used as a cheaper stability probe around 100 states:

| method | nstates | alpha0 | C6 |
|---|---:|---:|---:|
| TDA | 50 | 59.45612695 | 625.02150383 |
| TDA | 80 | 66.03916899 | 819.66912503 |
| TDA | 100 | 66.27449672 | 827.91519868 |
| TDA | 120 | 66.36259933 | 831.03236415 |

The TDA values stabilize near 100 states, but they are far above the full-TDDFT
100-state result. This method dependence means the 96% full-TDDFT closure cannot
yet be treated as a stable benchmark result.

## Sn q4 Closure With ANO Reference

Using the official matched PSP side:

```text
PSP        = GTH-PBE-q4
PSP basis  = TZV2P-MOLOPT-PBE-GTH-q4
C6_PSP     = 474.08122525
```

The PSP-RPA side is stable with respect to `nstates`:

| nstates | C6 PSP |
|---:|---:|
| 50 | 474.07715190 |
| 100 | 474.08123434 |
| 150 | 474.08122537 |
| 200 | 474.08122525 |

Two EFT `4d` channel constructions were tested:

| EFT basis | C6 PSP+EFT | EFT delta | closure |
|---|---:|---:|---:|
| def2-TZVP | 760.53316244 | 286.45601135 | 298.74% |
| ano | 566.13108046 | 92.05392937 | 96.00% |

The `def2-TZVP` 4d channel overcorrects badly. The `ano` 4d channel is much
better matched to the ANO all-electron reference and nearly closes the gap at
`nstates=100`.

Using the converged PSP value and the ANO `4d` EFT delta, the closure depends
strongly on the reference:

| reference | C6 all-e | closure |
|---|---:|---:|
| TDDFT ANO nstates=100 | 569.96394236 | 96.01% |
| TDDFT ANO nstates=150 | 576.98408849 | 89.46% |
| TDA ANO nstates=100 | 827.91519868 | 26.02% |
| TDA ANO nstates=120 | 831.03236415 | 25.79% |

Thus the PSP side is not the bottleneck. The benchmark is limited by the
all-electron reference method and state convergence.

## Interpretation

Sn q4 is now the most physically promising no-generated-basis candidate found so
far.

However, it should still be reported with a reference-method caveat. The
full-TDDFT ANO reference gives strong closure from 100 to 150 states, while TDA
gives a much larger all-electron C6 and therefore a smaller closure fraction.

Current label:

```text
Sn q4 = strongest official-matched secondary candidate, reference-method caveat
```

## Recommended Next Step

Optional next confirmation, if compute time allows:

```text
nstates = 200
```

Then recompute closure with the `ano` 4d EFT channels. The existing 80-150
state full-TDDFT trend is already sufficient to label Sn q4 as the strongest
official-matched secondary candidate found so far.

For now, the honest range is:

```text
Sn q4 closure is about 89-96% under full TDDFT ANO between 100 and 150 states,
but about 26% under TDA. The full-TDDFT trend supports Sn q4 as a strong
candidate, while the TDA discrepancy remains a reference-method caveat.
```
