# Official Matched Candidate Scan

## Scope

This scan only uses existing basis and pseudopotential data. No generated basis
is used.

Candidates:

```text
Al q3, Si q4, P q5, S q6, Cl q7
Ge q4, Sn q4, Pb q4
```

The PSP side uses CP2K/GTH official matched PBE basis/pseudopotential pairs when
available. The all-electron side uses PySCF built-in all-electron basis sets.

## Main Group Ne-Core Candidates

These candidates keep `3s/3p` explicit and treat `2s/2p` as the frozen core
response. Because Al, P, and Cl are open-shell odd-electron atoms, their
correction was evaluated with the closed-shell core-ion multipole TDHF route.

| atom | q | C6 PSP | C6 PSP+EFT | C6 all-e | EFT closure |
|---|---:|---:|---:|---:|---:|
| Al | 3 | 243.05351530 | 245.41012255 | 579.04752745 | 0.70% |
| Si | 4 | 149.19506801 | 150.63999198 | 337.70144156 | 0.77% |
| P | 5 | 125.04243104 | 126.24552611 | 201.16269394 | 1.58% |
| S | 6 | 70.38294960 | 71.02337170 | 139.09741354 | 0.93% |
| Cl | 7 | 80.32314934 | 80.84330247 | 68.88274213 | -4.55% |

These are clean official-matched controls, but the frozen Ne-core response is
too deep and too compact to close much of the C6 gap.

## Ge q4

Ge q4 remains the best official-matched control in this scan.

Using the better all-electron reference from the Ge convergence study:

```text
PSP basis       = TZV2P-MOLOPT-PBE-GTH-q4
all-e basis     = aug-cc-pVQZ
all-e nstates   = 200
EFT shell       = 3d
```

Result:

```text
C6 PSP      = 306.65944141
C6 PSP+EFT  = 317.45391425
C6 all-e    = 375.57206997
closure     = 15.66%
```

This is a positive official-matched result, but it is still much weaker than Mg
q2.

## Sn And Pb q4

Sn and Pb are physically tempting because the frozen `4d`/`5d` shells are
shallower than Ge `3d`.

Sn q4 has now been repaired with PySCF's built-in `ano` all-electron basis and
is the strongest official-matched secondary candidate found so far:

```text
Sn GTH-PBE-q4 / TZV2P-MOLOPT-PBE-GTH-q4
explicit PSP shells = 5s,5p
EFT shell           = 4d
C6_PSP              = 474.0812
C6_PSP+EFT          = 566.1352
C6_all-e_TDDFT150   = 576.9841
closure             = 89.46%
```

The caveat is reference-method dependence: TDA/ANO near 100 states gives a much
larger all-electron C6 and a closure closer to 26%. Full TDDFT/ANO is stable
from 80 to 150 states at the few-percent level.

Sn all-electron PBE TDDFT probe:

| basis | SCF converged | C6 |
|---|---|---:|
| def2-SVP | false | 190.77137368 |
| def2-TZVP | true | 0.19510512 |
| def2-TZVPP | true | 0.19469673 |
| def2-QZVPP | true | 162.20705558 |

Pb all-electron PBE TDDFT probe:

| basis | SCF converged | C6 |
|---|---|---:|
| def2-SVP | error | n/a |
| def2-TZVP | false | 302.56342218 |
| def2-TZVPP | false | 403.59550566 |
| def2-QZVPP | true | 15.97704256 |

These are not usable benchmark references without a separate heavy-element
reference strategy.

## Ranking

| rank | candidate | label | reason |
|---:|---|---|---|
| 1 | Sn q4 | strongest official-matched secondary candidate | clean, strong full-TDDFT closure, reference-method caveat |
| 2 | Ge q4 | best official-matched control | clean, positive closure, reference partially audited |
| 3 | P q5 | clean weak control | largest Ne-core closure among Al-Cl, but only 1.6% |
| 4 | Si/S q4/q6 | clean weak controls | stable but about 1% closure |
| 5 | Pb q4 | physically tempting, not ready | relativistic/heavy-reference issues |

## Recommendation

Best use:

```text
Sn q4 = strongest official-matched secondary candidate, reference-method caveat
Ge q4 = official-matched secondary control
Mg q2 = primary strong benchmark
Ca/Sr/Ba q2 = physically best route, but requires generated/imported q2 basis
```

If the goal is a strong benchmark rather than a control, the search should
return to alkaline-earth large-core q2 with a protocol-frozen generated or
imported basis.
