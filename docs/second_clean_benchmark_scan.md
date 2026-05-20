# Second Clean Benchmark Scan

Scan date: 2026-05-21

Goal: find a second clean large-core q2 PSP benchmark so Mg q2 is not the only
clean case.

## Scanner

The scan is implemented in:

```text
probe_large_core_q2_candidates.py
```

It searches CP2K potential and basis files for exact `q2` aliases, excluding
`q20`, `q21`, etc. For each candidate it records whether PySCF can build the
molecule, converge RKS, and run a small TDDFT smoke test when a matched q2 basis
exists.

Primary output:

```text
results/large_core_q2_candidate_scan.csv
```

Positive-control output:

```text
results/large_core_q2_candidate_scan_mg_control.csv
```

## Positive Control

Mg q2 is found in the current CP2K data and passes TDDFT smoke for several q2
bases:

```text
Mg GTH-PBE-q2 / TZV2P-MOLOPT-SR-GTH-q2 -> tddft_smoke_ok
Mg GTH-LDA-q2 / TZV2P-MOLOPT-SR-GTH-q2 -> tddft_smoke_ok
Mg GTH-PADE-q2 / TZV2P-MOLOPT-SR-GTH-q2 -> tddft_smoke_ok
Mg GTH-BLYP-q2 / TZV2P-MOLOPT-SR-GTH-q2 -> tddft_smoke_ok
```

This confirms the scanner can find a real matched q2 pseudo/basis route when
one exists.

## Target Candidate Results

Atoms scanned:

```text
Ca, Sr, Ba, Zn, Cd, Hg
```

Current result:

```text
atom  q2 pseudo found  matched q2 basis found  status
Ca    yes              no                      no_matched_q2_basis
Sr    yes              no                      no_matched_q2_basis
Ba    yes              no                      no_matched_q2_basis
Zn    yes              no                      no_matched_q2_basis
Cd    yes              no                      no_matched_q2_basis
Hg    yes              no                      no_matched_q2_basis
```

Important bug fixed during the scan:

```text
q20 basis aliases were initially being matched as q2.
```

The scanner now uses exact q-token matching, so `q20` is not accepted as `q2`.

## Conclusion

The local CP2K data currently do not provide a second matched q2 large-core
PSP+basis benchmark for Ca/Sr/Ba/Zn/Cd/Hg.

This means the next clean benchmark cannot be obtained by simply reusing the
currently imported CP2K files. The best next options are:

1. Import an external matched large-core q2 basis for Ca, Sr, or Ba.
2. Generate or adapt a q2-compatible basis and keep it clearly labeled as
   constructed, not library-native.
3. Use the existing Ca q2 clean-by-shell-overlap result only as a diagnostic,
   not as the second production clean benchmark.

Recommended next target remains Ca or Sr because their semicore physics is
closest to the Mg/Ca motivation. Cd/Hg q2 pseudos exist, but no matched q2 basis
was found locally, and their physical interpretation is less direct.
