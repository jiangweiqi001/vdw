# Semicore C6 Reproducibility Manifest

This manifest separates durable inputs/results from regenerable run
intermediates for the current semicore C6 work.

## Durable Scripts

Keep these scripts in the repository:

- `compute_core_sternheimer.py`
- `semicore_c6_targets.py`
- `probe_large_core_q2_candidates.py`
- `run_psp_rpa_atom.py`
- `run_semicore_c6_workflow.py`
- `run_semicore_c6_validation.py`
- `run_semicore_sr_validation.py`
- `run_semicore_reference_errors.py`
- `run_semicore_basis_sensitivity.py`
- `run_core_sternheimer_basis_diagnostic.py`
- `compute_core_sternheimer_radial_grid.py`
- `run_radial_cphf_alpha_constrained_validation.py`
- `run_alpha_constrained_uncertainty.py`
- `run_semicore_diagnostic_summary.py`
- `scripts/setup_wsl_pyscf_env.sh`
- `scripts/run_semicore_sr_smoke.sh`
- `scripts/run_semicore_sr_validation.sh`

## Durable Input Tables

Keep these compact input/reference tables:

- `external_data/cp2k/large_core_q2_basis_candidates_zn_cd.csv`
- `reference_pair_c6.csv`
- `reference_pair_c6_alternates.csv`
- `reference_static_polarizability_group12.csv`

## Durable Documentation

Keep these documents:

- `PROJECT_STATUS.md`
- `docs/student_task_semicore_c6.md`
- `docs/student_task_semicore_c6_checklist.md`
- `docs/semicore_reference_policy.md`
- `docs/finite_basis_sternheimer_note.md`
- `docs/semicore_reproducibility_manifest.md`

## Durable Result Summaries

Keep these summary CSVs because they are cited by the status/policy docs:

- `results/semicore_c6_q2_candidate_scan_wsl.csv`
- `results/semicore_c6_workflow_targets_wsl.csv`
- `results/semicore_c6_q2_candidate_scan_zn_cd_adapted.csv`
- `results/semicore_c6_workflow_targets_zn_cd_adapted.csv`
- `results/semicore_zn_cd_basis_sensitivity/summary.csv`
- `results/semicore_zn_cd_basis_sensitivity/reference_errors.csv`
- `results/semicore_zn_cd_validation/diagnostic_summary.csv`
- `results/semicore_zn_cd_validation/multi_reference_errors.csv`
- `results/semicore_zn_cd_validation/trk_normalized_summary.csv`
- `results/semicore_zn_cd_validation/cd_literature_closure_proxy.csv`
- `results/semicore_zn_cd_validation/all_e_controls/attempt_summary.csv`
- `results/semicore_zn_cd_validation/validation/zn/summary.csv`
- `results/semicore_zn_cd_validation/validation/cd/summary.csv`
- `results/semicore_zn_cd_validation/validation/zn_trk_normalized/summary.csv`
- `results/semicore_zn_cd_validation/validation/cd_trk_normalized/summary.csv`
- `results/core_sternheimer_basis_diagnostic/summary.csv`
- `results/core_sternheimer_sr_diagnostic/summary.csv`
- `results/semicore_sr_validation/diagnostic_summary.csv`
- `results/semicore_sr_validation/trk_normalized_summary.csv`
- `results/semicore_sr_validation/sr_sanity_summary.csv`
- `results/semicore_sr_validation/validation/summary.csv`
- `results/semicore_sr_validation/validation_trk_normalized/summary.csv`
- `results/semicore_sr_validation/validation_ano_rcc_core/summary.csv`
- `results/semicore_sr_validation/validation_ano_rcc_core_trk_normalized/summary.csv`
- `results/radial_grid_sternheimer/zn_3d_summary.csv`
- `results/radial_grid_sternheimer/cd_4d_summary.csv`
- `results/radial_grid_sternheimer/zn_3d_cphf_alpha_c6_table.csv`
- `results/radial_grid_sternheimer/cd_4d_cphf_alpha_c6_table.csv`
- `results/radial_grid_sternheimer/psp_rpa_summary.csv`
- `results/radial_grid_sternheimer/validation/zn/summary.csv`
- `results/radial_grid_sternheimer/validation/cd/summary.csv`
- `results/radial_grid_sternheimer/multi_reference_errors.csv`
- `results/radial_grid_sternheimer/alpha_constrained/zn/kernel_summary.csv`
- `results/radial_grid_sternheimer/alpha_constrained/zn/validation/summary.csv`
- `results/radial_grid_sternheimer/alpha_constrained/cd/kernel_summary.csv`
- `results/radial_grid_sternheimer/alpha_constrained/cd/validation/summary.csv`
- `results/radial_grid_sternheimer/alpha_constrained/multi_reference_errors.csv`
- `results/radial_grid_sternheimer/alpha_constrained/uncertainty_summary.csv`

## Regenerable Intermediates

Do not keep these in git unless debugging a specific numerical issue:

- PSP channel directories under `results/**/psp_rpa/`
- core-channel directories under `results/**/core_sternheimer/`
- diagnostic channel directories under `results/**/channels/`
- all-electron channel run directories under
  `results/**/all_e_controls/*/`
- validation channel and alpha/C6 intermediate tables under
  `results/**/validation/**/`
- radial-grid CPHF channel files
  `results/radial_grid_sternheimer/*_cphf_channels.csv` unless a specific
  audit requires retaining the discretized continuum channels

These are ignored by `.gitignore` and can be recreated with the durable scripts.

## Current Cleanup Decision

The repository should track compact summaries and policy documents, not the
large channel-level CSVs. If a future result needs channel-level auditability,
regenerate it from the scripts and include only the specific file needed for the
audit.
