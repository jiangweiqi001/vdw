import argparse
import csv
from pathlib import Path

from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from semicore_c6_targets import audit_semicore_target


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_rows(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def combine_channel_files(psp_channels_path, core_channels_path, output_path):
    psp_rows = read_rows(psp_channels_path)
    core_rows = read_rows(core_channels_path)
    fieldnames = sorted(set().union(*(row.keys() for row in psp_rows + core_rows)))
    write_rows(output_path, psp_rows + core_rows, fieldnames)
    return output_path


def validation_summary_row(atom, psp_row, corrected_row, active_electrons, active_shells, all_e_row=None, reference_c6=None):
    c6_psp = float(psp_row["C6_self_au"])
    c6_corrected = float(corrected_row["C6_self_au"])
    delta_core = c6_corrected - c6_psp
    audit = audit_semicore_target(atom, active_electrons, active_shells)

    row = {
        "atom": atom,
        "dimer": audit["dimer"],
        "active_electrons": int(active_electrons),
        "active_shells": audit["active_shells"],
        "correction_shells": audit["correction_shells"],
        "C6_PSP": c6_psp,
        "C6_PSP_plus_sternheimer": c6_corrected,
        "Delta_C6_core": delta_core,
        "double_counting_status": "clean" if audit["audit_status"] == "pass" else "fail",
        "target_audit_status": audit["audit_status"],
        "C6_all_e": "",
        "closure_pct": "",
        "C6_reference": "",
        "reference_error_pct": "",
        "go_no_go": "",
    }

    if all_e_row:
        c6_all_e = float(all_e_row["C6_self_au"])
        missing = c6_all_e - c6_psp
        row["C6_all_e"] = c6_all_e
        row["closure_pct"] = 100.0 * delta_core / missing if missing else 0.0

    if reference_c6 is not None and str(reference_c6).strip():
        ref = float(reference_c6)
        row["C6_reference"] = ref
        row["reference_error_pct"] = 100.0 * (c6_corrected - ref) / ref
        row["go_no_go"] = "go" if atom in {"Zn", "Cd"} and abs(row["reference_error_pct"]) <= 10.0 else "review"

    return row


def run_validation(atom, psp_channels, core_channels, output_root, active_electrons, active_shells, all_e_channels=None, reference_c6=None):
    output_root = Path(output_root)
    combined_path = output_root / f"{atom.lower()}_psp_plus_sternheimer_channels.csv"
    combine_channel_files(psp_channels, core_channels, combined_path)

    psp_alpha_path = output_root / "psp_alpha_c6_table.csv"
    corrected_alpha_path = output_root / "psp_plus_sternheimer_alpha_c6_table.csv"
    corrected_c6_path = output_root / "psp_plus_sternheimer_c6_table.csv"
    write_alpha_rows(psp_alpha_path, build_alpha_rows(psp_channels))
    corrected_alpha_rows = build_alpha_rows(combined_path)
    write_alpha_rows(corrected_alpha_path, corrected_alpha_rows)
    write_c6_rows(corrected_c6_path, build_c6_rows(combined_path))

    all_e_row = None
    if all_e_channels:
        all_e_rows = build_alpha_rows(all_e_channels)
        all_e_row = next(row for row in all_e_rows if row["atom"] == atom)

    psp_row = next(row for row in build_alpha_rows(psp_channels) if row["atom"] == atom)
    corrected_row = next(row for row in corrected_alpha_rows if row["atom"] == atom)
    summary = validation_summary_row(
        atom=atom,
        psp_row=psp_row,
        corrected_row=corrected_row,
        active_electrons=active_electrons,
        active_shells=active_shells,
        all_e_row=all_e_row,
        reference_c6=reference_c6,
    )
    write_rows(output_root / "summary.csv", [summary], list(summary.keys()))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate large-core PSP + non-empirical semicore Sternheimer C6.")
    parser.add_argument("--atom", required=True, choices=["Sr", "Zn", "Cd"])
    parser.add_argument("--psp-channels", required=True)
    parser.add_argument("--core-channels", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--active-electrons", required=True, type=int)
    parser.add_argument("--active-shells", required=True)
    parser.add_argument("--all-e-channels")
    parser.add_argument("--reference-c6")
    args = parser.parse_args(argv)

    summary = run_validation(
        atom=args.atom,
        psp_channels=args.psp_channels,
        core_channels=args.core_channels,
        output_root=args.output_root,
        active_electrons=args.active_electrons,
        active_shells=args.active_shells,
        all_e_channels=args.all_e_channels,
        reference_c6=args.reference_c6,
    )
    print(",".join(summary.keys()))
    print(",".join(str(value) for value in summary.values()))
    return summary


if __name__ == "__main__":
    main()
