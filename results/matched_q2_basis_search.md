# Matched q2 Basis Search

## Summary

No matched Ca/Sr q2 basis was found in the currently available CP2K/PySCF data.

The current conclusion is:

```text
matched Ca/Sr q2 basis not found in available CP2K/PySCF data
```

The project should therefore treat Mg q2 as the main clean benchmark for now.

## What Was Found

### Mg q2

Mg has a matched q2 route:

```text
pseudopotential: GTH-PBE-q2
basis: TZV2P-MOLOPT-SR-GTH-q2
explicit valence: 3s
EFT correction shells: 2s,2p
double counting: clean
```

This is the current strongest benchmark:

```text
C6_PSP = 638.6202
C6_PSP+EFT = 647.6996
C6_all-e_PBE = 647.5881
```

### Ca / Sr / Ba q2

Available pseudopotentials include:

```text
Ca GTH-LDA-q2 / GTH-PADE-q2
Sr GTH-LDA-q2 / GTH-PADE-q2
Ba GTH-LDA-q2 / GTH-PADE-q2 / GTH-BLYP-q2
```

But matched q2 MOLOPT basis files were not found in the available CP2K/PySCF
data. Available basis sets are generally q10:

```text
Ca q10 basis
Sr q10 basis
Ba q10 basis
```

Thus the current Ca q2 adapted row remains diagnostic:

```text
GTH-PBE-q2 pseudo + TZV2P-MOLOPT-PBE-GTH-q10 basis
```

This is a pseudo-basis mismatch and should not be reported as a final matched
Ca q2 benchmark.

## Current Main Benchmark

Use Mg q2 as the main clean result:

```text
results/mg_q2_clean_benchmark.md
results/mg_q2/summary.csv
results/mg_q2/benchmark.png
```

## Candidate Next Directions

### Option A: K / Rb q1

Potentially strong physics:

```text
K q1: explicit valence = 4s; frozen shell = 3s,3p
Rb q1: explicit valence = 5s; frozen shell = 4s,4p
```

Pros:

- Closest to the alkali/semicore physics emphasized by the PRL.
- Frozen shallow core response may be large.
- Could provide a strong second benchmark.

Cons:

- Open-shell / single-valence response is more complicated.
- Current scripts mostly assume closed-shell RKS/TDDFT.
- Requires careful spin treatment and possibly ROHF/UKS/TDDFT variants.

### Option B: Zn / Cd q2

Potentially clean q2 candidates:

```text
Zn q2
Cd q2
```

Pros:

- q2 pseudopotentials exist in available data.
- Could be closed-shell-like in the explicit valence sector.

Cons:

- d-shell physics is more complex.
- The frozen d-shell or shallow semicore response may require more careful
  partitioning.
- Need to verify matched basis availability and response stability.

### Option C: Generate Ca / Sr q2 Basis

Most directly aligned with the desired Ca/Sr story:

```text
Ca q2: explicit valence = 4s; EFT shells = 3s,3p
Sr q2: explicit valence = 5s; EFT shells = 4s,4p
```

Pros:

- Best physical match to the semicore EFT-core benchmark story.
- Would provide a second strong result if successful.

Cons:

- Requires basis generation or external basis import.
- CP2K `optimize_basis` / ATOM route may be needed.
- More engineering effort than using existing data.

## Recommended Next Step

Do not spend more time blindly searching for Ca/Sr q2 basis files.

Recommended order:

1. Keep Mg q2 as the headline benchmark.
2. Add prior-art comparison and figures around Mg q2.
3. For a second example, first scan whether Zn/Cd q2 can run cleanly with matched
   available data.
4. If a stronger physical example is needed, prototype K/Rb q1 with explicit
   open-shell treatment.
5. Generate a matched Ca/Sr q2 basis only if the project specifically needs an
   alkaline-earth semicore production benchmark.
