# Matched q3/q4 Candidate Scan

## Web And Local Availability

Online CP2K documentation and data indicate that `qXX` suffixes denote the
matching GTH pseudopotential valence. The relevant official sources are:

- CP2K `BASIS_MOLOPT`
- CP2K `GTH_POTENTIALS`
- CP2K basis-set documentation

The local CP2K data confirms official matched q3/q4 availability for:

```text
Al q3
Si q4
Ge q4
Sn q4
Pb q4
```

The local availability table previously listed Ga/In/Tl q3 candidates, but the
current local files only contain their q3 pseudopotentials, not usable matched
q3 basis blocks. They should not be treated as matched candidates unless a real
basis source is imported.

## PSP-RPA Smoke Results

Matched PSP-TDDFT smoke tests ran for:

| atom | q | status | note |
|---|---:|---|---|
| Al | 3 | OK | open-shell q3 |
| Si | 4 | OK | open-shell q4 |
| Ge | 4 | OK | d-shell frozen candidate |
| Sn | 4 | OK | d-shell frozen candidate |
| Pb | 4 | OK | d-shell frozen candidate |

Al and Si do not expose a shallow frozen d semicore shell, so they are lower
priority for the EFT-core benchmark story.

## Semicore Smoke Results

The d-shell q4 candidates were tested with additive l=1 MO dipole Wilson
corrections:

| atom | frozen shell | C6 PSP | C6 PSP+EFT | C6 all-e | closure |
|---|---|---:|---:|---:|---:|
| Ge | 3d | 306.59756732 | 317.38856671 | 243.20214317 | -17.02% |
| Sn | 4d | 474.07715109 | 760.53316244 | 0.57567417 | invalid |
| Pb | 5d | 510.25085626 | 510.25085626 | 524.42956310 | 0.00% |

Interpretation:

- `Ge q4` is clean/matched but not a strong positive closure case in the first
  smoke test. The PSP C6 is already larger than the non-aug all-electron C6, so
  the d-shell correction moves in the wrong direction.
- `Sn q4` all-electron `def2-TZVP` TDDFT produced an unphysical C6 and should be
  marked invalid until the all-electron reference is fixed.
- `Pb q4` produced no useful 5d correction through the current MO-shell route.

## Ge Reference Sensitivity

Additional Ge all-electron checks show significant basis sensitivity:

| basis | C6 all-e |
|---|---:|
| def2-TZVP | 241.06648751 |
| def2-TZVPP | 241.06649404 |
| cc-pVTZ | 238.40341881 |
| aug-cc-pVTZ | 368.83617137 |

Thus Ge q4 cannot be promoted without a focused all-electron reference
convergence study.

## Current Recommendation

No q3/q4 candidate currently beats Mg q2 as a second strong benchmark.

Best next options:

1. Treat `Ge q4` as the most plausible official-matched follow-up, but first run
   a dedicated all-electron basis/diffuse/nstates convergence study.
2. Keep `Sn q4` and `Pb q4` as lower priority until all-electron reference and
   shell-mapping issues are fixed.
3. Do not use Ga/In/Tl q3 as matched candidates unless matched q3 basis blocks
   are imported from an official source.
