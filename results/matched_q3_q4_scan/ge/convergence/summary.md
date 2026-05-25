# Ge q4 All-Electron PBE TDDFT Convergence

## Purpose

This checks whether Ge q4 can serve as an official-matched second benchmark.
The PSP side is clean and matched:

```text
pseudopotential = GTH-PBE-q4
PSP basis       = TZV2P-MOLOPT-PBE-GTH-q4
explicit shell  = 4s, 4p
candidate EFT shell = 3d
```

The issue is the all-electron PBE TDDFT reference.

## Basis And Diffuse Scan

All rows use open-shell `spin=2`, PBE TDDFT, `nstates=100`.

| basis | nao | alpha0 | C6 |
|---|---:|---:|---:|
| def2-SVP | 32 | 22.16874373 | 148.59528424 |
| def2-TZVP | 48 | 29.28790225 | 241.07363762 |
| def2-TZVPP | 48 | 29.28606189 | 241.06648617 |
| cc-pVTZ | 43 | 29.48770005 | 238.40342410 |
| aug-cc-pVTZ | 59 | 44.28623330 | 368.83606287 |
| cc-pVQZ | 68 | 36.24707916 | 306.47193557 |
| aug-cc-pVQZ | 93 | 44.38304460 | 361.12518226 |

Diffuse functions strongly change the reference. Triple-zeta non-augmented
bases give `C6 ~ 238-241`, while augmented QZ gives `C6 ~ 361` at the same
state count.

## nstates Scan

### cc-pVQZ

| nstates | alpha0 | C6 |
|---:|---:|---:|
| 50 | 35.80528989 | 295.44889956 |
| 100 | 36.24708042 | 306.47193588 |
| 150 | 36.36786338 | 309.67864769 |
| 200 | 36.43406222 | 311.50890305 |

### aug-cc-pVQZ

| nstates | alpha0 | C6 |
|---:|---:|---:|
| 50 | 41.75569786 | 307.74900715 |
| 100 | 44.38304444 | 361.12516558 |
| 150 | 44.60893704 | 367.09938494 |
| 200 | 44.91520257 | 375.57206997 |

The augmented QZ reference is not fully converged at 200 states. The change from
150 to 200 states is about 2.3% in C6.

## Implication For Ge q4 Benchmark

The original smoke test used a non-augmented all-electron reference and found:

```text
C6_PSP      = 306.59756732
C6_PSP+EFT  = 317.38856671
C6_all-e    = 243.20214317
closure     = -17.0%
```

With the better augmented QZ all-electron reference, the reference C6 is closer
to:

```text
C6_all-e ~ 375.57 at aug-cc-pVQZ, nstates=200
```

Against that reference, the same PSP+EFT result would close only about:

```text
(317.38856671 - 306.59756732) / (375.57206997 - 306.59756732)
  = 15.65%
```

This is a positive correction, but still much weaker than Mg q2.

## Recommendation

Do not promote Ge q4 to the second strong benchmark yet.

Ge q4 is the best official-matched q3/q4 follow-up candidate found so far, but
it needs more reference work:

1. Extend `aug-cc-pVQZ` to higher `nstates` if affordable.
2. Check whether `aug-cc-pV5Z` exists and is tractable.
3. Audit the `3d` EFT correction route, because the correction is modest.
4. Keep Mg q2 as the headline benchmark.

Current label:

```text
Ge q4 = official-matched candidate, promising control, not strong benchmark
```
