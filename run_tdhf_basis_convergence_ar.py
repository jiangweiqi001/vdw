import argparse
import csv
from pathlib import Path

from compare_alpha_c6 import compare_tables
from pyscf_export_ar_tdhf_oscillators import BASIS_LIST, export_ar_tdhf_oscillators, write_tdhf_channels
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows


def _comparison_for_ar(alpha_table):
    rows = compare_tables(alpha_table, "reference_alpha_c6.csv")
    return next(row for row in rows if row["atom"] == "Ar")


def run_one_basis(basis, output_root, nstates=100, min_osc=1e-10):
    result_dir = Path(output_root) / basis.lower()
    result_dir.mkdir(parents=True, exist_ok=True)
    channels_path = result_dir / "ar_tdhf_channels.csv"
    alpha_path = result_dir / "alpha_c6_table.csv"
    c6_path = result_dir / "c6_table.csv"

    channel_rows, tdhf_summary = export_ar_tdhf_oscillators(basis, nstates=nstates, min_osc=min_osc)
    write_tdhf_channels(channels_path, channel_rows)
    alpha_rows = build_alpha_rows(channels_path)
    write_alpha_rows(alpha_path, alpha_rows)
    write_c6_rows(c6_path, build_c6_rows(channels_path))
    comparison = _comparison_for_ar(alpha_path)

    return {
        "basis": basis,
        "nstates_requested": tdhf_summary["nstates_requested"],
        "nstates_returned": tdhf_summary["nstates_returned"],
        "n_channels": tdhf_summary["n_channels"],
        "sum_osc": tdhf_summary["sum_osc"],
        "alpha0": comparison["alpha0_eft"],
        "C6": comparison["C6_eft"],
        "alpha0_error_pct": comparison["err_alpha_pct"],
        "C6_error_pct": comparison["err_C6_pct"],
    }


def write_summary(path, rows):
    fieldnames = [
        "basis",
        "nstates_requested",
        "nstates_returned",
        "n_channels",
        "sum_osc",
        "alpha0",
        "C6",
        "alpha0_error_pct",
        "C6_error_pct",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_tdhf_basis_convergence(bases=None, output_root="results/ar/tdhf", nstates=100, min_osc=1e-10):
    rows = [run_one_basis(basis, output_root, nstates=nstates, min_osc=min_osc) for basis in (bases or BASIS_LIST)]
    write_summary(Path(output_root).parent / "ar_tdhf_basis_convergence.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ar TDHF oscillator-strength basis convergence.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST)
    parser.add_argument("--output-root", default="results/ar/tdhf")
    parser.add_argument("--nstates", type=int, default=100)
    parser.add_argument("--min-osc", type=float, default=1e-10)
    args = parser.parse_args(argv)

    rows = run_tdhf_basis_convergence(args.basis or BASIS_LIST, args.output_root, args.nstates, args.min_osc)
    print("basis,nstates_returned,n_channels,sum_osc,alpha0,C6,alpha0_error_pct,C6_error_pct")
    for row in rows:
        print(
            f"{row['basis']},"
            f"{row['nstates_returned']},"
            f"{row['n_channels']},"
            f"{row['sum_osc']:.8f},"
            f"{row['alpha0']:.8f},"
            f"{row['C6']:.8f},"
            f"{row['alpha0_error_pct']:.6f},"
            f"{row['C6_error_pct']:.6f}"
        )


if __name__ == "__main__":
    main()
