# EFT-vdW Implementation Roadmap

This roadmap reflects the current project state after the Mg q2 clean benchmark.

The project has moved from general response prototyping to a concrete
PSP+EFT-core benchmark:

```text
PSP valence response
  + l=1 frozen-core dipole channels
  -> alpha(i xi)
  -> C6
  -> long-range tail / model-screened pairwise energy
```

The current strongest result is:

```text
Mg q2:
C6_PSP = 638.62015545
C6_PSP + EFT = 647.69960692
C6_all-e PBE = 647.58810040
double counting = clean
```

The roadmap below focuses on making this result independently reproducible,
then finding a second strong large-core benchmark.

## 1. Solidify Mg q2 As The Main Result

Status: partly done.

Already available:

```text
results/mg_q2_clean_benchmark.md
results/mg/mg2_tail_comparison.csv
results/mg/mg_q2_model_screening_sensitivity.csv
```

Next deliverables:

```text
run_mg_q2_clean_benchmark.py
plot_mg_q2_results.py
results/mg_q2/
```

`run_mg_q2_clean_benchmark.py` should regenerate from scratch:

```text
results/mg_q2/summary.csv
results/mg_q2/tail.csv
results/mg_q2/screening_sensitivity.csv
```

`plot_mg_q2_results.py` should produce:

```text
results/mg_q2/benchmark.png
```

Acceptance criteria:

- One command regenerates PSP q2 response, EFT dipole channels, PSP+EFT C6,
  all-electron PBE reference, Mg2 tail, and model-screening sensitivity.
- `results/mg_q2/summary.csv` contains the headline closure value.
- `README.md` links to the Mg q2 result page.

## 2. Find A Matched Ca Or Sr q2 Large-Core Benchmark

Status: open.

Current Ca adapted result is promising but not final:

```text
Ca_q2_PBE_adapted:
GTH-PBE-q2 pseudo + TZV2P-MOLOPT-PBE-GTH-q10 basis
closure = 22.76%
double counting = clean
```

The issue is a pseudo-basis mismatch: q2 pseudopotential with q10 basis.

Target clean benchmark:

```text
Ca q2: explicit valence = 4s; EFT shells = 3s,3p
Sr q2: explicit valence = 5s; EFT shells = 4s,4p
```

Next options:

1. Find or import a matched Ca q2 basis.
2. Search for a matched Sr q2 pseudopotential/basis pair.
3. If neither exists, document that an external PSP generation route is needed.

Acceptance criteria:

- SCF + TDDFT runs for a matched q2 PSP/basis.
- Double-counting status is `clean`.
- The benchmark closes a meaningful part of the PSP-to-all-electron C6 gap.

## 3. Add Prior-Art Comparison

Status: not done.

Create:

```text
results/prior_art_comparison.md
```

Include:

- Stuttgart CPP
- Marinescu/Dalgarno/Babb model-potential + CPP
- Derevianko alkali C6 core contribution
- Dutt RCCSD `alpha_c + alpha_cv + alpha_v`
- Current Mg q2 result

Purpose:

```text
The novelty is not "core polarizability affects C6".
The novelty is PSP + first-principles EFT-core correction + double-counting control.
```

Acceptance criteria:

- The document clearly states what prior art already did.
- The Mg q2 result is positioned as an EFT-style implementation benchmark, not
  as the first observation of core polarizability.

## 4. Model W_v Expansion

Status: interface exists.

Already implemented:

```text
screened_pairwise_vdw.py
screened_eft_vdw.py
```

Current model kernels:

```text
bare
dielectric
yukawa
```

Next optional deliverables:

```text
plot_model_screening_sensitivity.py
results/mg_q2/screening_sensitivity.png
```

Acceptance criteria:

- Plots clearly state "model W_v sensitivity; not ab initio screening".
- Bare limit reproduces `-C6/R^6`.

## 5. Ab Initio W_v

Status: future work.

Do not start until Mg q2 is fully packaged and at least one additional matched
large-core benchmark is attempted.

Target:

```text
W_v(i xi) = [v^-1 - chi_v^irr(i xi)]^-1
```

This requires:

- a valence response representation
- real-space or reciprocal-space screened interaction
- local-field effects
- double-counting projection
- eventually periodic systems and forces

Acceptance criteria:

- Starts from the same PSP valence response used in the benchmark.
- Reduces to the current bare/model-screened result in the appropriate limit.

## 6. Current Priority

The next concrete task is:

```text
Implement run_mg_q2_clean_benchmark.py.
```

This will turn the strongest result into a one-command reproducible experiment.
