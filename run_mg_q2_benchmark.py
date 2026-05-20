import argparse
import csv
from pathlib import Path

from compute_dipole_wilson import compute_dipole_wilson_channels, write_channels as write_dipole_channels
from run_all_e_rpa_atom import run_atom as run_all_e_atom
from run_all_e_rpa_atom import write_summary as write_all_e_summary
from run_all_e_vs_psp_rpa_summary import build_summary as build_all_e_vs_psp_summary
from run_all_e_vs_psp_rpa_summary import write_summary as write_all_e_vs_psp_summary
from run_eft_core_dipole_validation import BEST_PSP, run_validation
from run_psp_rpa_atom import run_atom as run_psp_atom
from run_psp_rpa_atom import write_summary as write_psp_summary


MG_Q2_PSP_PATH = "results/psp_rpa/mg/GTH-PBE-q2_TZV2P-MOLOPT-SR-GTH-q2_pbe_tddft/mg_psp_channels.csv"


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_rows(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def upsert_row(path, row, keys, writer):
    existing = read_rows(path) if Path(path).exists() else []
    filtered = [old for old in existing if not all(str(old.get(key)) == str(row.get(key)) for key in keys)]
    filtered.append(row)
    writer(path, filtered)


def select_one(rows, **criteria):
    matches = [
        row for row in rows
        if all(str(row.get(key)) == str(value) for key, value in criteria.items())
    ]
    if not matches:
        raise ValueError(f"No row matching {criteria}")
    if len(matches) > 1:
        raise ValueError(f"Multiple rows matching {criteria}")
    return matches[0]


def select_mg_q2_psp_row(rows):
    return select_one(
        rows,
        atom="Mg",
        psp="GTH-PBE-q2",
        basis="TZV2P-MOLOPT-SR-GTH-q2",
        xc="pbe",
        method="TDDFT",
        nstates="100",
    )


def select_mg_pbe_all_e_row(rows):
    return select_one(
        rows,
        atom="Mg",
        xc="pbe",
        basis="aug-cc-pVQZ",
        method="TDDFT",
        nstates="200",
    )


def sorted_shells(value):
    return sorted(shell for shell in str(value).split(";") if shell)


def benchmark_row(psp_row, all_e_row, eft_row):
    c6_psp = float(psp_row["C6_psp"])
    c6_all_e = float(all_e_row["C6"])
    c6_eft = float(eft_row["C6_psp_plus_dipole"])
    delta_psp = c6_all_e - c6_psp
    delta_eft = c6_all_e - c6_eft
    closure = (c6_eft - c6_psp) / delta_psp if delta_psp else 0.0
    clean = (
        str(eft_row.get("double_counting_status")) == "clean"
        and set(sorted_shells(eft_row.get("correction_shells", ""))).isdisjoint(
            set(sorted_shells(eft_row.get("psp_explicit_valence_shells", "")))
        )
    )
    return {
        "atom": "Mg",
        "psp": psp_row["psp"],
        "psp_basis": psp_row["basis"],
        "xc": psp_row["xc"],
        "method": psp_row["method"],
        "active_electrons": int(psp_row["active_electrons"]),
        "active_shells": psp_row["active_shells"],
        "eft_shells": eft_row["correction_shells"],
        "C6_psp": c6_psp,
        "C6_psp_plus_eft": c6_eft,
        "C6_all_e": c6_all_e,
        "delta_C6_missing": delta_psp,
        "residual_C6": delta_eft,
        "closure_fraction": closure,
        "closure_pct": 100.0 * closure,
        "double_counting_status": eft_row["double_counting_status"],
        "benchmark_status": "clean_candidate" if clean else "needs_audit",
        "note": "unscreened l=1 MO dipole EFT approximation",
    }


def audit_row(psp_path, psp_row, all_e_row, eft_row):
    correction_shells = set(sorted_shells(eft_row.get("correction_shells", "")))
    explicit_shells = set(sorted_shells(eft_row.get("psp_explicit_valence_shells", "")))
    shell_overlap = bool(correction_shells & explicit_shells)
    checks = {
        "mg_q2_two_electron": str(psp_row.get("active_electrons")) == "2",
        "psp_active_shell_3s": str(psp_row.get("active_shells")) == "3s",
        "correction_shells_2s2p": correction_shells == {"2s", "2p"},
        "no_shell_overlap": not shell_overlap,
        "psp_pbe_tddft": psp_row.get("xc") == "pbe" and psp_row.get("method") == "TDDFT",
        "all_e_pbe_tddft": all_e_row.get("xc") == "pbe" and all_e_row.get("method") == "TDDFT",
        "placeholder_path_not_used": "placeholder" not in str(psp_path),
        "double_counting_clean": eft_row.get("double_counting_status") == "clean",
    }
    return {
        "atom": "Mg",
        "psp_path": str(psp_path),
        "placeholder_path_used": str(not checks["placeholder_path_not_used"]).lower(),
        "shell_overlap": str(shell_overlap).lower(),
        "audit_status": "pass" if all(checks.values()) else "fail",
        **{key: str(value).lower() for key, value in checks.items()},
    }


def sensitivity_rows(shell_rows, virtual_rows):
    rows = []
    for row in shell_rows:
        if row.get("atom") != "Mg":
            continue
        rows.append(
            {
                "atom": "Mg",
                "case_type": "shell",
                "case": row["shell"],
                "C6_psp": float(row["C6_psp"]),
                "Delta_C6": float(row["Delta_C6_shell"]),
            }
        )
    for row in virtual_rows:
        if row.get("atom") != "Mg":
            continue
        rows.append(
            {
                "atom": "Mg",
                "case_type": "virtual_cutoff",
                "case": row["max_delta_Ha"],
                "C6_psp": float(row["C6_psp"]),
                "Delta_C6": float(row["Delta_C6_cutoff"]),
            }
        )
    return rows


def stability_row(category, case, value, reference, tolerance_pct, quantity):
    value = float(value)
    reference = float(reference)
    delta = value - reference
    delta_pct = 100.0 * delta / reference if reference else 0.0
    return {
        "category": category,
        "case": case,
        "quantity": quantity,
        "value": value,
        "reference": reference,
        "delta": delta,
        "delta_pct": delta_pct,
        "tolerance_pct": float(tolerance_pct),
        "status": "pass" if abs(delta_pct) <= float(tolerance_pct) else "review",
    }


def write_stability_markdown(path, rows):
    review_rows = [row for row in rows if row["status"] != "pass"]
    lines = [
        "# Mg q2 Sensitivity Audit",
        "",
        "This audit checks whether the near-complete Mg q2 closure is stable under",
        "basic numerical choices for the all-electron reference and PSP baseline.",
        "",
        "## Summary",
        "",
        f"Total cases: {len(rows)}",
        f"Review cases: {len(review_rows)}",
        "",
        "## Cases",
        "",
        "```text",
        "category,case,quantity,value,reference,delta_pct,status",
    ]
    for row in rows:
        lines.append(
            f"{row['category']},{row['case']},{row['quantity']},"
            f"{row['value']:.8f},{row['reference']:.8f},{row['delta_pct']:.6f},{row['status']}"
        )
    lines.extend(["```", ""])
    if review_rows:
        lines.extend([
            "Cases marked `review` are not automatic failures; they identify numerical",
            "choices that shift the benchmark beyond the configured tolerance and should",
            "not be used as the headline value without explanation.",
            "",
        ])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def run_full_sensitivity(bench):
    rows = []
    all_e_reference = float(bench["C6_all_e"])
    psp_reference = float(bench["C6_psp"])

    for nstates in [100, 150, 200]:
        row = run_all_e_atom("Mg", "pbe", "aug-cc-pVQZ", nstates, "TDDFT")
        upsert_row(
            "results/all_e_rpa_summary.csv",
            row,
            keys=["atom", "xc", "basis", "nstates", "method"],
            writer=write_all_e_summary,
        )
        rows.append(stability_row("all_e_nstates", f"aug-cc-pVQZ_nstates_{nstates}", row["C6"], all_e_reference, 0.2, "C6_all_e"))

    for basis in ["aug-cc-pVTZ", "aug-cc-pVQZ"]:
        row = run_all_e_atom("Mg", "pbe", basis, 200, "TDDFT")
        upsert_row(
            "results/all_e_rpa_summary.csv",
            row,
            keys=["atom", "xc", "basis", "nstates", "method"],
            writer=write_all_e_summary,
        )
        rows.append(stability_row("all_e_basis", f"{basis}_nstates_200", row["C6"], all_e_reference, 0.5, "C6_all_e"))

    for nstates in [20, 50, 100]:
        row = run_psp_atom(
            atom="Mg",
            psp="GTH-PBE-q2",
            basis="TZV2P-MOLOPT-SR-GTH-q2",
            xc="pbe",
            nstates=nstates,
            method="TDDFT",
            basis_file="external_data/cp2k/BASIS_MOLOPT_UCL",
            basis_name="TZV2P-MOLOPT-SR-GTH-q2",
            pseudo_file="external_data/cp2k/GTH_POTENTIALS",
            pseudo_name="GTH-PBE-q2",
        )
        upsert_row(
            "results/psp_rpa_summary.csv",
            row,
            keys=["atom", "psp", "basis", "xc", "nstates", "method"],
            writer=write_psp_summary,
        )
        rows.append(stability_row("psp_nstates", f"nstates_{nstates}", row["C6_psp"], psp_reference, 0.2, "C6_psp"))

    write_rows(
        "results/mg_q2_stability_audit.csv",
        rows,
        ["category", "case", "quantity", "value", "reference", "delta", "delta_pct", "tolerance_pct", "status"],
    )
    write_stability_markdown("docs/mg_q2_sensitivity_audit.md", rows)
    return rows


def write_benchmark_markdown(path, bench, audit):
    text = f"""# Mg q2 PSP + EFT-Core Benchmark Audit

Audit date: 2026-05-21

## Benchmark Chain

```text
C6_PSP             = {bench['C6_psp']:.8f}
C6_PSP+dipole_EFT  = {bench['C6_psp_plus_eft']:.8f}
C6_all-e_PBE_TDDFT = {bench['C6_all_e']:.8f}
closure_fraction   = {bench['closure_fraction']:.12f}
residual_C6         = {bench['residual_C6']:.8f}
```

## Audit Status

```text
audit_status                  = {audit['audit_status']}
psp_path                      = {audit['psp_path']}
placeholder_path_used          = {audit['placeholder_path_used']}
active_shells                 = {bench['active_shells']}
eft_shells                    = {bench['eft_shells']}
double_counting_status         = {bench['double_counting_status']}
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
"""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8")


def load_current_rows(
    psp_summary="results/psp_rpa_summary.csv",
    all_e_summary="results/all_e_rpa_summary.csv",
    eft_summary="results/eft_core_dipole_validation_summary.csv",
):
    psp_row = select_mg_q2_psp_row(read_rows(psp_summary))
    all_e_row = select_mg_pbe_all_e_row(read_rows(all_e_summary))
    eft_row = select_one(read_rows(eft_summary), atom="Mg")
    return psp_row, all_e_row, eft_row


def rerun_inputs():
    psp_row = run_psp_atom(
        atom="Mg",
        psp="GTH-PBE-q2",
        basis="TZV2P-MOLOPT-SR-GTH-q2",
        xc="pbe",
        nstates=100,
        method="TDDFT",
        basis_file="external_data/cp2k/BASIS_MOLOPT_UCL",
        basis_name="TZV2P-MOLOPT-SR-GTH-q2",
        pseudo_file="external_data/cp2k/GTH_POTENTIALS",
        pseudo_name="GTH-PBE-q2",
    )
    upsert_row(
        "results/psp_rpa_summary.csv",
        psp_row,
        keys=["atom", "psp", "basis", "xc", "nstates", "method"],
        writer=write_psp_summary,
    )

    all_e_row = run_all_e_atom("Mg", "pbe", "aug-cc-pVQZ", 200, "TDDFT")
    upsert_row(
        "results/all_e_rpa_summary.csv",
        all_e_row,
        keys=["atom", "xc", "basis", "nstates", "method"],
        writer=write_all_e_summary,
    )

    dipole_rows = []
    for atom, basis, shells in [
        ("Mg", "cc-pVQZ", {"2s", "2p"}),
        ("Ca", "cc-pVQZ", {"3s", "3p"}),
    ]:
        dipole_rows.extend(compute_dipole_wilson_channels(atom, basis, shells))
    write_dipole_channels("results/eft_core_dipole_wilson_channels.csv", dipole_rows)

    comparison_rows = build_all_e_vs_psp_summary()
    write_all_e_vs_psp_summary("results/all_e_vs_psp_rpa_summary.csv", comparison_rows)
    run_validation()


def run_benchmark(skip_rerun=False):
    if not skip_rerun:
        rerun_inputs()
    psp_row, all_e_row, eft_row = load_current_rows()
    bench = benchmark_row(psp_row, all_e_row, eft_row)
    audit = audit_row(MG_Q2_PSP_PATH, psp_row, all_e_row, eft_row)
    write_rows(
        "results/mg_q2_benchmark_summary.csv",
        [bench],
        [
            "atom",
            "psp",
            "psp_basis",
            "xc",
            "method",
            "active_electrons",
            "active_shells",
            "eft_shells",
            "C6_psp",
            "C6_psp_plus_eft",
            "C6_all_e",
            "delta_C6_missing",
            "residual_C6",
            "closure_fraction",
            "closure_pct",
            "double_counting_status",
            "benchmark_status",
            "note",
        ],
    )
    write_rows(
        "results/mg_q2_benchmark_audit.csv",
        [audit],
        [
            "atom",
            "psp_path",
            "placeholder_path_used",
            "shell_overlap",
            "audit_status",
            "mg_q2_two_electron",
            "psp_active_shell_3s",
            "correction_shells_2s2p",
            "no_shell_overlap",
            "psp_pbe_tddft",
            "all_e_pbe_tddft",
            "placeholder_path_not_used",
            "double_counting_clean",
        ],
    )
    write_benchmark_markdown("docs/mg_q2_benchmark_audit.md", bench, audit)
    sens = sensitivity_rows(
        read_rows("results/eft_core_dipole_shell_summary.csv"),
        read_rows("results/eft_core_dipole_virtual_convergence.csv"),
    )
    write_rows(
        "results/mg_q2_sensitivity_summary.csv",
        sens,
        ["atom", "case_type", "case", "C6_psp", "Delta_C6"],
    )
    return bench, audit


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run and audit the clean Mg q2 PSP + EFT-core benchmark.")
    parser.add_argument("--skip-rerun", action="store_true", help="Use existing CSV artifacts instead of rerunning PySCF jobs.")
    parser.add_argument("--full-sensitivity", action="store_true", help="Run Mg q2 all-e/PSP numerical stability checks.")
    args = parser.parse_args(argv)
    bench, audit = run_benchmark(skip_rerun=args.skip_rerun)
    if args.full_sensitivity:
        run_full_sensitivity(bench)
    print("atom,C6_psp,C6_psp_plus_eft,C6_all_e,closure_fraction,residual_C6,audit_status,benchmark_status")
    print(
        f"{bench['atom']},{bench['C6_psp']:.8f},{bench['C6_psp_plus_eft']:.8f},"
        f"{bench['C6_all_e']:.8f},{bench['closure_fraction']:.12f},"
        f"{bench['residual_C6']:.8f},{audit['audit_status']},{bench['benchmark_status']}"
    )


if __name__ == "__main__":
    main()
