import argparse
import csv
from pathlib import Path

from compare_alpha_c6 import compare_tables
from pyscf_export_ar_tdhf_oscillators import export_ar_tdhf_oscillators, write_tdhf_channels
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows


BASIS_LIST = ["aug-cc-pVTZ", "aug-cc-pVQZ"]
NSTATES_LIST = [20, 50, 100, 150, 200]


def summarize_nstates_result(basis, nstates, tdhf_summary, comparison, n_electrons=18.0):
    return {
        "basis": basis,
        "nstates": int(nstates),
        "n_channels": int(tdhf_summary["n_channels"]),
        "sum_osc": float(tdhf_summary["sum_osc"]),
        "sum_osc_over_N": float(tdhf_summary["sum_osc"]) / n_electrons if n_electrons else 0.0,
        "alpha0": float(comparison["alpha0_eft"]),
        "C6": float(comparison["C6_eft"]),
        "alpha0_err": float(comparison["err_alpha_pct"]),
        "C6_err": float(comparison["err_C6_pct"]),
    }


def _comparison_for_ar(alpha_table):
    rows = compare_tables(alpha_table, "reference_alpha_c6.csv")
    return next(row for row in rows if row["atom"] == "Ar")


def run_one_case(basis, nstates, output_root, min_osc=1e-10):
    case_dir = Path(output_root) / basis.lower() / f"nstates_{nstates}"
    case_dir.mkdir(parents=True, exist_ok=True)
    channels_path = case_dir / "ar_tdhf_channels.csv"
    alpha_path = case_dir / "alpha_c6_table.csv"
    c6_path = case_dir / "c6_table.csv"

    channel_rows, tdhf_summary = export_ar_tdhf_oscillators(basis, nstates=nstates, min_osc=min_osc)
    write_tdhf_channels(channels_path, channel_rows)
    write_alpha_rows(alpha_path, build_alpha_rows(channels_path))
    write_c6_rows(c6_path, build_c6_rows(channels_path))
    return summarize_nstates_result(basis, nstates, tdhf_summary, _comparison_for_ar(alpha_path))


def write_summary(path, rows):
    fieldnames = [
        "basis",
        "nstates",
        "n_channels",
        "sum_osc",
        "sum_osc_over_N",
        "alpha0",
        "C6",
        "alpha0_err",
        "C6_err",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_nstates_convergence(bases=None, nstates_values=None, output_root="results/ar/tdhf_nstates", min_osc=1e-10):
    rows = []
    for basis in bases or BASIS_LIST:
        for nstates in nstates_values or NSTATES_LIST:
            rows.append(run_one_case(basis, int(nstates), output_root, min_osc=min_osc))
    write_summary(Path(output_root).parent / "ar_tdhf_nstates_convergence.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ar TDHF nstates convergence.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST)
    parser.add_argument("--nstates", action="append", type=int)
    parser.add_argument("--output-root", default="results/ar/tdhf_nstates")
    parser.add_argument("--min-osc", type=float, default=1e-10)
    args = parser.parse_args(argv)

    rows = run_nstates_convergence(args.basis or BASIS_LIST, args.nstates or NSTATES_LIST, args.output_root, args.min_osc)
    print("basis,nstates,alpha0,C6,alpha0_err,C6_err,sum_osc,sum_osc/N")
    for row in rows:
        print(
            f"{row['basis']},"
            f"{row['nstates']},"
            f"{row['alpha0']:.8f},"
            f"{row['C6']:.8f},"
            f"{row['alpha0_err']:.6f},"
            f"{row['C6_err']:.6f},"
            f"{row['sum_osc']:.8f},"
            f"{row['sum_osc_over_N']:.8f}"
        )


if __name__ == "__main__":
    main()
