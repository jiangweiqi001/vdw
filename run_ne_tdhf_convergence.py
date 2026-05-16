import argparse
import csv
from pathlib import Path

from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from run_noble_gas_tdhf import load_references, percent_error
from run_tdhf_atom import export_tdhf_atom, write_channels


BASIS_LIST = ["aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z"]
NSTATES_LIST = [100, 200, 300]


def summarize_case(
    basis,
    nstates,
    alpha0,
    c6,
    n_channels,
    sum_osc,
    alpha0_ref,
    c6_ref,
):
    return {
        "atom": "Ne",
        "basis": basis,
        "nstates": int(nstates),
        "alpha0": float(alpha0),
        "C6": float(c6),
        "alpha0_err": percent_error(alpha0, alpha0_ref),
        "C6_err": percent_error(c6, c6_ref),
        "sum_osc": float(sum_osc),
        "n_channels": int(n_channels),
        "status": "ok",
        "note": "",
    }


def unavailable_row(basis, nstates, reason):
    return {
        "atom": "Ne",
        "basis": basis,
        "nstates": int(nstates),
        "alpha0": "",
        "C6": "",
        "alpha0_err": "",
        "C6_err": "",
        "sum_osc": "",
        "n_channels": "",
        "status": "unavailable",
        "note": str(reason),
    }


def run_one_case(basis, nstates, output_root, refs):
    case_dir = Path(output_root) / basis.lower() / f"nstates_{nstates}"
    case_dir.mkdir(parents=True, exist_ok=True)
    channels_path = case_dir / "ne_tdhf_channels.csv"
    alpha_path = case_dir / "alpha_c6_table.csv"
    c6_path = case_dir / "c6_table.csv"
    try:
        channel_rows, tdhf_summary = export_tdhf_atom("Ne", basis, nstates=nstates)
    except Exception as exc:
        return unavailable_row(basis, nstates, f"{type(exc).__name__}: {exc}")

    write_channels(channels_path, channel_rows)
    alpha_rows = build_alpha_rows(channels_path)
    write_alpha_rows(alpha_path, alpha_rows)
    write_c6_rows(c6_path, build_c6_rows(channels_path))
    alpha_row = next(row for row in alpha_rows if row["atom"] == "Ne")
    return summarize_case(
        basis=basis,
        nstates=nstates,
        alpha0=float(alpha_row["alpha0_au"]),
        c6=float(alpha_row["C6_self_au"]),
        n_channels=int(alpha_row["n_channels"]),
        sum_osc=tdhf_summary["sum_osc"],
        alpha0_ref=refs["Ne"]["alpha0_ref"],
        c6_ref=refs["Ne"]["C6_ref"],
    )


def write_summary(path, rows):
    fieldnames = [
        "atom",
        "basis",
        "nstates",
        "alpha0",
        "C6",
        "alpha0_err",
        "C6_err",
        "sum_osc",
        "n_channels",
        "status",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_ne_convergence(bases=None, nstates_values=None, output_root="results/ne/tdhf_convergence"):
    refs = load_references()
    rows = []
    for basis in bases or BASIS_LIST:
        for nstates in nstates_values or NSTATES_LIST:
            rows.append(run_one_case(basis, int(nstates), output_root, refs))
    write_summary(Path(output_root).parent / "ne_tdhf_convergence.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ne TDHF basis/nstates convergence.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST)
    parser.add_argument("--nstates", action="append", type=int)
    parser.add_argument("--output-root", default="results/ne/tdhf_convergence")
    args = parser.parse_args(argv)

    rows = run_ne_convergence(args.basis or BASIS_LIST, args.nstates or NSTATES_LIST, args.output_root)
    print("atom,basis,nstates,alpha0,C6,alpha0_err,C6_err,sum_osc,n_channels,status,note")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['basis']},"
            f"{row['nstates']},"
            f"{row['alpha0']},"
            f"{row['C6']},"
            f"{row['alpha0_err']},"
            f"{row['C6_err']},"
            f"{row['sum_osc']},"
            f"{row['n_channels']},"
            f"{row['status']},"
            f"{row['note']}"
        )


if __name__ == "__main__":
    main()
