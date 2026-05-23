import argparse
import csv
from pathlib import Path

from compute_dipole_wilson import compute_dipole_wilson_channels, write_channels as write_dipole_channels
from run_all_e_rpa_atom import run_atom as run_all_e_atom
from run_all_e_rpa_atom import write_summary as write_all_e_summary
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from run_eft_core_dipole_validation import run_validation
from run_mg_q2_model_screening_sensitivity import build_sensitivity_rows, write_rows as write_screening_rows
from run_psp_rpa_atom import run_atom as run_psp_atom
from run_psp_rpa_atom import write_summary as write_psp_summary


MG_Q2_PSP_PATH = "results/psp_rpa/mg/GTH-PBE-q2_TZV2P-MOLOPT-SR-GTH-q2_pbe_tddft/mg_psp_channels.csv"


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def upsert_row(path, row, keys, writer):
    existing = read_rows(path) if Path(path).exists() else []
    filtered = [old for old in existing if not all(str(old.get(key)) == str(row.get(key)) for key in keys)]
    filtered.append(row)
    writer(path, filtered)


def summary_row(c6_psp, c6_eft, c6_all_e, double_counting_status):
    c6_psp = float(c6_psp)
    c6_eft = float(c6_eft)
    c6_all_e = float(c6_all_e)
    missing = c6_all_e - c6_psp
    correction = c6_eft - c6_psp
    residual = c6_all_e - c6_eft
    return {
        "case": "Mg_q2",
        "C6_PSP": c6_psp,
        "C6_PSP_EFT": c6_eft,
        "C6_all_e": c6_all_e,
        "Delta_C6_missing": missing,
        "Delta_C6_EFT": correction,
        "residual_C6": residual,
        "closure_fraction": correction / missing if missing else 0.0,
        "closure_pct": 100.0 * correction / missing if missing else 0.0,
        "double_counting_status": double_counting_status,
    }


def write_summary(path, rows):
    fieldnames = [
        "case",
        "C6_PSP",
        "C6_PSP_EFT",
        "C6_all_e",
        "Delta_C6_missing",
        "Delta_C6_EFT",
        "residual_C6",
        "closure_fraction",
        "closure_pct",
        "double_counting_status",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_tail(path, c6_all, c6_psp, c6_eft):
    HARTREE_TO_MEV = 27211.386245988
    rows = []
    for r in [8, 10, 12, 15, 20, 30, 40]:
        row = {"R_bohr": r}
        for name, c6 in [("all_e_pbe", c6_all), ("psp_q2", c6_psp), ("psp_plus_eft", c6_eft)]:
            e = -float(c6) / r**6
            row[f"E_{name}_Ha"] = e
            row[f"E_{name}_meV"] = e * HARTREE_TO_MEV
        row["err_psp_pct"] = 100.0 * (float(c6_psp) - float(c6_all)) / float(c6_all)
        row["err_eft_pct"] = 100.0 * (float(c6_eft) - float(c6_all)) / float(c6_all)
        rows.append(row)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_benchmark(output_dir="results/mg_q2"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_e_row = run_all_e_atom("Mg", "pbe", "aug-cc-pVQZ", 200, "TDDFT")
    upsert_row(
        "results/all_e_rpa_summary.csv",
        all_e_row,
        keys=["atom", "xc", "basis", "nstates", "method"],
        writer=write_all_e_summary,
    )

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

    correction_rows = compute_dipole_wilson_channels("Mg", "aug-cc-pVQZ", {"2s", "2p"})
    write_dipole_channels("results/eft_core_dipole_wilson_channels.csv", correction_rows)
    eft_row = next(row for row in run_validation() if row["atom"] == "Mg")

    summary = summary_row(
        c6_psp=eft_row["C6_psp"],
        c6_eft=eft_row["C6_psp_plus_dipole"],
        c6_all_e=all_e_row["C6"],
        double_counting_status=eft_row["double_counting_status"],
    )
    write_summary(output_dir / "summary.csv", [summary])
    write_tail(output_dir / "tail.csv", summary["C6_all_e"], summary["C6_PSP"], summary["C6_PSP_EFT"])
    write_screening_rows(output_dir / "screening_sensitivity.csv", build_sensitivity_rows(summary["C6_PSP_EFT"]))

    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Reproduce the Mg q2 clean PSP+EFT-core benchmark.")
    parser.add_argument("--output-dir", default="results/mg_q2")
    args = parser.parse_args(argv)

    row = run_benchmark(args.output_dir)
    print("case,C6_PSP,C6_PSP_EFT,C6_all_e,closure_pct,double_counting_status")
    print(
        f"{row['case']},"
        f"{row['C6_PSP']:.8f},"
        f"{row['C6_PSP_EFT']:.8f},"
        f"{row['C6_all_e']:.8f},"
        f"{row['closure_pct']:.6f},"
        f"{row['double_counting_status']}"
    )


if __name__ == "__main__":
    main()
