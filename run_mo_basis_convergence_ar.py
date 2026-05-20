import argparse
import csv
from pathlib import Path

from compare_alpha_c6 import compare_tables
from pyscf_export_ar_mo_oscillators import BASIS_LIST, export_ar_mo_oscillators, write_mo_channels
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows


def _comparison_for_ar(alpha_table):
    rows = compare_tables(alpha_table, "reference_alpha_c6.csv")
    return next(row for row in rows if row["atom"] == "Ar")


def run_one_basis(basis, output_root):
    result_dir = Path(output_root) / basis.lower()
    result_dir.mkdir(parents=True, exist_ok=True)
    channels_path = result_dir / "ar_mo_channels.csv"
    alpha_path = result_dir / "alpha_c6_table.csv"
    c6_path = result_dir / "c6_table.csv"

    channel_rows, trk = export_ar_mo_oscillators(basis)
    write_mo_channels(channels_path, channel_rows)
    alpha_rows = build_alpha_rows(channels_path)
    write_alpha_rows(alpha_path, alpha_rows)
    write_c6_rows(c6_path, build_c6_rows(channels_path))
    comparison = _comparison_for_ar(alpha_path)

    return {
        "basis": basis,
        "sum_f_mo_over_N": trk["sum_f_mo_over_N_electrons"],
        "sum_f_mo": trk["sum_f_mo"],
        "N_electrons": trk["N_electrons"],
        "alpha0": comparison["alpha0_eft"],
        "C6": comparison["C6_eft"],
        "alpha0_error_pct": comparison["err_alpha_pct"],
        "C6_error_pct": comparison["err_C6_pct"],
        "n_occ": trk["n_occ"],
        "n_virt": trk["n_virt"],
        "n_channels": len(channel_rows),
    }


def write_summary(path, rows):
    fieldnames = [
        "basis",
        "sum_f_mo_over_N",
        "sum_f_mo",
        "N_electrons",
        "alpha0",
        "C6",
        "alpha0_error_pct",
        "C6_error_pct",
        "n_occ",
        "n_virt",
        "n_channels",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_mo_basis_convergence(bases=None, output_root="results/ar/mo"):
    rows = [run_one_basis(basis, output_root) for basis in (bases or BASIS_LIST)]
    write_summary(Path(output_root).parent / "ar_mo_basis_convergence.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ar 3D MO oscillator-strength basis convergence.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST)
    parser.add_argument("--output-root", default="results/ar/mo")
    args = parser.parse_args(argv)

    rows = run_mo_basis_convergence(args.basis or BASIS_LIST, args.output_root)
    print("basis,sum_f_mo/N,alpha0,C6,alpha0_error_pct,C6_error_pct,n_occ,n_virt")
    for row in rows:
        print(
            f"{row['basis']},"
            f"{row['sum_f_mo_over_N']:.8f},"
            f"{row['alpha0']:.8f},"
            f"{row['C6']:.8f},"
            f"{row['alpha0_error_pct']:.6f},"
            f"{row['C6_error_pct']:.6f},"
            f"{row['n_occ']},"
            f"{row['n_virt']}"
        )


if __name__ == "__main__":
    main()
