import argparse
import csv
from pathlib import Path

from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows


BEST_PSP = {
    "Mg": "results/psp_rpa/mg/GTH-PBE-q2_TZV2P-MOLOPT-SR-GTH-q2_pbe_tddft/mg_psp_channels.csv",
    "Ca": "results/psp_rpa/ca/gth-lda_cc-pVQZ_lda_tddft/ca_psp_channels.csv",
}

PSP_EXPLICIT_VALENCE_SHELLS = {
    # q2 large-core Mg keeps only 3s explicit.
    "Mg": {"3s"},
    # q2 large-core Ca keeps only 4s explicit.
    "Ca": {"4s"},
}

DELTA_CUTOFFS_HA = [5.0, 10.0, 20.0, 50.0]


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_combined(path, rows):
    fieldnames = ["atom", "channel", "delta_Ha", "osc", "is_core", "source"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})


def correction_shells(rows):
    return {row.get("from_shell", "") for row in rows if row.get("from_shell", "")}


def double_counting_status(correction_shell_set, explicit_valence_shells):
    if explicit_valence_shells is None:
        return "unknown"
    return (
        "diagnostic_double_counting"
        if set(correction_shell_set) & set(explicit_valence_shells)
        else "clean"
    )


def enforce_double_counting(atom, correction_rows, allow_diagnostic):
    status = double_counting_status(correction_shells(correction_rows), PSP_EXPLICIT_VALENCE_SHELLS.get(atom))
    if status == "diagnostic_double_counting" and not allow_diagnostic:
        overlap = sorted(correction_shells(correction_rows) & PSP_EXPLICIT_VALENCE_SHELLS.get(atom, set()))
        raise ValueError(
            f"{atom} correction shells overlap PSP explicit valence shells: {overlap}. "
            "Use --diagnostic-double-counting to generate diagnostic output."
        )
    return status


def c6_for(path, prefix):
    alpha_rows = build_alpha_rows(path)
    write_alpha_rows(f"{prefix}_alpha_c6_table.csv", alpha_rows)
    write_c6_rows(f"{prefix}_c6_table.csv", build_c6_rows(path))
    return {row["atom"]: float(row["C6_self_au"]) for row in alpha_rows}


def _validate_atom(atom, psp_path, correction_rows, output_root):
    combined_path = output_root / atom.lower() / f"{atom.lower()}_psp_plus_dipole_channels.csv"
    combined = read_csv(psp_path) + correction_rows
    write_combined(combined_path, combined)
    c6_psp = c6_for(psp_path, output_root / atom.lower() / "psp")[atom]
    c6_total = c6_for(combined_path, output_root / atom.lower() / "psp_plus_dipole")[atom]
    return c6_psp, c6_total


def shell_summary(atom, psp_path, correction_rows, output_root):
    rows = []
    c6_psp = c6_for(psp_path, output_root / atom.lower() / "psp_shell_base")[atom]
    for shell in sorted(correction_shells(correction_rows)):
        shell_rows = [row for row in correction_rows if row.get("from_shell") == shell]
        _psp, c6_total = _validate_atom(atom, psp_path, shell_rows, output_root / "shells" / shell)
        rows.append(
            {
                "atom": atom,
                "shell": shell,
                "n_channels": len(shell_rows),
                "sum_osc": sum(float(row["osc"]) for row in shell_rows),
                "C6_psp": c6_psp,
                "C6_psp_plus_shell": c6_total,
                "Delta_C6_shell": c6_total - c6_psp,
            }
        )
    return rows


def virtual_convergence(atom, psp_path, correction_rows, output_root):
    rows = []
    c6_psp = c6_for(psp_path, output_root / atom.lower() / "psp_virtual_base")[atom]
    for cutoff in DELTA_CUTOFFS_HA:
        filtered = [
            row for row in correction_rows
            if float(row["delta_Ha"]) <= cutoff
        ]
        _psp, c6_total = _validate_atom(atom, psp_path, filtered, output_root / "virtual_cutoffs" / f"{cutoff:g}")
        rows.append(
            {
                "atom": atom,
                "max_delta_Ha": cutoff,
                "n_channels": len(filtered),
                "sum_osc": sum(float(row["osc"]) for row in filtered),
                "C6_psp": c6_psp,
                "C6_psp_plus_cutoff": c6_total,
                "Delta_C6_cutoff": c6_total - c6_psp,
            }
        )
    return rows


def write_table(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_validation(
    input_path="results/eft_core_dipole_wilson_channels.csv",
    output="results/eft_core_dipole_validation_summary.csv",
    allow_diagnostic_double_counting=False,
):
    dipole_rows = read_csv(input_path)
    by_atom = {}
    for row in dipole_rows:
        by_atom.setdefault(row["atom"], []).append(row)
    summary = []
    shell_rows = []
    virtual_rows = []
    output_root = Path("results/eft_core_dipole")
    for atom, correction_rows in sorted(by_atom.items()):
        psp_path = BEST_PSP.get(atom)
        if not psp_path or not Path(psp_path).exists():
            continue
        status = enforce_double_counting(atom, correction_rows, allow_diagnostic_double_counting)
        c6_psp, c6_total = _validate_atom(atom, psp_path, correction_rows, output_root)
        summary.append(
            {
                "atom": atom,
                "C6_psp": c6_psp,
                "C6_psp_plus_dipole": c6_total,
                "Delta_C6_dipole": c6_total - c6_psp,
                "n_dipole_channels": len(correction_rows),
                "correction_shells": ";".join(sorted(correction_shells(correction_rows))),
                "psp_explicit_valence_shells": ";".join(sorted(PSP_EXPLICIT_VALENCE_SHELLS.get(atom, []))),
                "double_counting_status": status,
                "note": "l=1 MO dipole Wilson approximation; unscreened additive test",
            }
        )
        shell_rows.extend(shell_summary(atom, psp_path, correction_rows, output_root))
        virtual_rows.extend(virtual_convergence(atom, psp_path, correction_rows, output_root))
    write_table(
        output,
        summary,
        [
            "atom",
            "C6_psp",
            "C6_psp_plus_dipole",
            "Delta_C6_dipole",
            "n_dipole_channels",
            "correction_shells",
            "psp_explicit_valence_shells",
            "double_counting_status",
            "note",
        ],
    )
    write_table(
        "results/eft_core_dipole_shell_summary.csv",
        shell_rows,
        ["atom", "shell", "n_channels", "sum_osc", "C6_psp", "C6_psp_plus_shell", "Delta_C6_shell"],
    )
    write_table(
        "results/eft_core_dipole_virtual_convergence.csv",
        virtual_rows,
        ["atom", "max_delta_Ha", "n_channels", "sum_osc", "C6_psp", "C6_psp_plus_cutoff", "Delta_C6_cutoff"],
    )
    return summary


def run_core_tdhf_validation(
    input_path="results/eft_core_tdhf_wilson_channels.csv",
    output="results/eft_core_tdhf_validation_summary.csv",
):
    rows = read_csv(input_path)
    by_atom = {}
    for row in rows:
        by_atom.setdefault(row["atom"], []).append(row)
    summary = []
    output_root = Path("results/eft_core_tdhf")
    for atom, correction_rows in sorted(by_atom.items()):
        psp_path = BEST_PSP.get(atom)
        if not psp_path or not Path(psp_path).exists():
            continue
        c6_psp, c6_total = _validate_atom(atom, psp_path, correction_rows, output_root)
        status = double_counting_status(set(), PSP_EXPLICIT_VALENCE_SHELLS.get(atom))
        summary.append(
            {
                "atom": atom,
                "C6_psp": c6_psp,
                "C6_psp_plus_core_tdhf": c6_total,
                "Delta_C6_core_tdhf": c6_total - c6_psp,
                "n_channels": len(correction_rows),
                "double_counting_status": status,
                "note": "core-ion TDHF transition-density proxy for EFT dipole Wilson",
            }
        )
    write_table(
        output,
        summary,
        ["atom", "C6_psp", "C6_psp_plus_core_tdhf", "Delta_C6_core_tdhf", "n_channels", "double_counting_status", "note"],
    )
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate l=1 dipole Wilson channels by adding them to PSP-RPA channels.")
    parser.add_argument("--input", default="results/eft_core_dipole_wilson_channels.csv")
    parser.add_argument("--output", default="results/eft_core_dipole_validation_summary.csv")
    parser.add_argument("--diagnostic-double-counting", action="store_true")
    args = parser.parse_args(argv)

    rows = run_validation(args.input, args.output, args.diagnostic_double_counting)
    print("atom,C6_psp,C6_psp_plus_dipole,Delta_C6_dipole,n_dipole_channels,double_counting_status,note")
    for row in rows:
        print(
            f"{row['atom']},{row['C6_psp']},{row['C6_psp_plus_dipole']},"
            f"{row['Delta_C6_dipole']},{row['n_dipole_channels']},"
            f"{row['double_counting_status']},{row['note']}"
        )


if __name__ == "__main__":
    main()
