# Alkali q1 Candidate Scan

This scan checks whether alkali q1 pseudopotentials can provide a second strong
PSP+EFT-core benchmark after Mg q2.

## Result

K and Rb q1 pseudopotentials are available in the current CP2K/PySCF data:

```text
K  GTH-LDA-q1
Rb GTH-LDA-q1
```

However, no matched q1 basis was found in the available CP2K/PySCF data. The
available basis sets are q9:

```text
TZV2P-MOLOPT-PBE-GTH-q9
TZVP-MOLOPT-PBE-GTH-q9
DZVP-MOLOPT-PBE-GTH-q9
TZV2P-MOLOPT-SR-GTH-q9
```

Using q1 pseudopotentials with q9 basis sets gives an adapted diagnostic route,
similar in spirit to the current Ca q2 adapted diagnostic.

## Smoke Test

K and Rb q1 adapted setups run with open-shell UKS + TDDFT/TDA:

```text
K  q1 pseudo + q9 basis: nelec = 1, UKS converges, TDDFT/TDA runs
Rb q1 pseudo + q9 basis: nelec = 1, UKS converges, TDDFT/TDA runs
```

Cs q1 with DZVP q9 basis runs, but TZV/TZV2P q9 smoke tests produced near-zero
oscillator strength, so Cs is not a good immediate candidate.

## Interpretation

K/Rb q1 are physically attractive because:

```text
K  explicit valence = 4s; EFT shells = 3s,3p
Rb explicit valence = 5s; EFT shells = 4s,4p
```

This is closer to the PRL alkali frozen-core dynamics than the noble-gas tests.

But there are two caveats:

1. The basis is q9, not matched q1.
2. The atoms are open-shell, so the response path requires UKS/ROKS and careful
   TDDFT/TDA spin handling.

## Recommendation

K/Rb q1 should be treated as the next adapted diagnostic direction, not yet as a
clean production benchmark.

Recommended next steps:

1. Implement an open-shell PSP response runner for UKS + TDA/TDDFT.
2. Run K q1 and Rb q1 adapted diagnostics.
3. Add double-counting metadata:

```text
K  explicit valence = 4s; EFT shells = 3s,3p
Rb explicit valence = 5s; EFT shells = 4s,4p
```

4. Search for or generate matched q1 basis if the adapted diagnostics are
   promising.
