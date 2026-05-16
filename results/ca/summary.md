# Ca Semicore Partition Benchmark

## Basis Note

PySCF's built-in basis library does not provide `aug-cc-pVTZ` or
`aug-cc-pVQZ` for Ca in this environment. The Ca smoke test used `cc-pVTZ`, and
the partition benchmark used all-electron `cc-pVQZ`.

## TDHF All-Electron Check

```text
Ca cc-pVTZ, nstates=200
alpha0 = 184.64905323
C6     = 2757.89012057
```

## cc-pVQZ Partition Result

The partition definition is:

```text
valence  = 4s
semicore = 3s,3p
deep_core = 1s,2s,2p
```

Key C6 values from `ca_core_valence_decomposition.csv`:

```text
C6_all                  = 2759.50968183
C6_valence_4s           = 2977.54172206
C6_semicore_3s3p        = 9.69475239
C6_deep_core            ~ 0
C6_valence_plus_semicore = 2761.47855918
```

The semicore effect relative to a 4s-only valence partition is:

```text
Delta_C6_semicore = C6[valence+semicore] - C6[valence]
                  = -216.06316288

relative_semicore_contribution = Delta_C6_semicore / C6_all
                               = -7.83%
```

Interpretation: in the TDHF projected response, Ca semicore mixing is materially
larger than the Ar/Kr noble-gas core correction. The effect is not a large
standalone semicore polarizability; it is mainly valence-semicore mixing that
reduces the valence-only C6.
