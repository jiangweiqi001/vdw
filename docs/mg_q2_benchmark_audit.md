# Mg q2 PSP + EFT-Core Benchmark Audit

Audit date: 2026-05-21

## Benchmark Chain

```text
C6_PSP             = 638.62015545
C6_PSP+dipole_EFT  = 647.60794451
C6_all-e_PBE_TDDFT = 647.58810040
closure_fraction   = 1.002212782316
residual_C6         = -0.01984411
```

## Audit Status

```text
audit_status                  = pass
psp_path                      = results/psp_rpa/mg/GTH-PBE-q2_TZV2P-MOLOPT-SR-GTH-q2_pbe_tddft/mg_psp_channels.csv
placeholder_path_used          = false
active_shells                 = 3s
eft_shells                    = 2p;2s
double_counting_status         = clean
```

The shell and virtual-cutoff sensitivity table is written to:

```text
results/mg_q2_sensitivity_summary.csv
```

The all-electron/PSP numerical stability audit is written to:

```text
results/mg_q2_stability_audit.csv
docs/mg_q2_sensitivity_audit.md
```

## Interpretation

Mg q2 remains a clean benchmark candidate: the PSP response has 2 active
electrons in the `3s` pseudo-valence space, while the added EFT dipole channels
come from `2s,2p`. The result should still be described as an unscreened `l=1`
MO dipole approximation, not as the final screened EFT-vdW functional.
