import argparse
import csv
from pathlib import Path

from run_all_e_rpa_atom import run_atom


BASIS_LIST = ["cc-pVTZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]
NSTATES_LIST = [50, 100, 150, 200]


def convergence_row(basis, nstates, summary, status="ok", note=""):
    return {
        "basis": basis,
        "nstates": int(nstates),
        "alpha0": summary.get("alpha0", ""),
        "C6": summary.get("C6", ""),
        "alpha0_err": summary.get("alpha0_err", ""),
        "C6_err": summary.get("C6_err", ""),
        "n_channels": summary.get("n_channels", ""),
        "status": status,
        "note": note,
    }


def unavailable_row(basis, nstates, reason):
    return convergence_row(basis, nstates, {}, status="unavailable", note=str(reason))


def run_one_case(basis, nstates, xc="lda", method="TDDFT", output_root="results/all_e_rpa_convergence/ar"):
    try:
        summary = run_atom(
            atom="Ar",
            xc=xc,
            basis=basis,
            nstates=nstates,
            method=method,
            output_root=Path(output_root) / basis.lower() / f"nstates_{nstates}",
        )
        return convergence_row(basis, nstates, summary)
    except Exception as exc:
        return unavailable_row(basis, nstates, f"{type(exc).__name__}: {exc}")


def write_summary(path, rows):
    fieldnames = [
        "basis",
        "nstates",
        "alpha0",
        "C6",
        "alpha0_err",
        "C6_err",
        "n_channels",
        "status",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_convergence(bases=None, nstates_values=None, xc="lda", method="TDDFT", output_root="results/all_e_rpa_convergence/ar"):
    rows = []
    for basis in bases or BASIS_LIST:
        for nstates in nstates_values or NSTATES_LIST:
            rows.append(run_one_case(basis, int(nstates), xc=xc, method=method, output_root=output_root))
    write_summary(Path(output_root).parent / "ar_all_e_rpa_convergence.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ar all-electron RKS TDDFT basis/nstates convergence.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST)
    parser.add_argument("--nstates", action="append", type=int)
    parser.add_argument("--xc", choices=["lda", "pbe"], default="lda")
    parser.add_argument("--method", choices=["TDDFT", "TDA"], default="TDDFT")
    parser.add_argument("--output-root", default="results/all_e_rpa_convergence/ar")
    args = parser.parse_args(argv)

    rows = run_convergence(args.basis or BASIS_LIST, args.nstates or NSTATES_LIST, args.xc, args.method, args.output_root)
    print("basis,nstates,alpha0,C6,alpha0_err,C6_err,n_channels,status,note")
    for row in rows:
        print(
            f"{row['basis']},"
            f"{row['nstates']},"
            f"{row['alpha0']},"
            f"{row['C6']},"
            f"{row['alpha0_err']},"
            f"{row['C6_err']},"
            f"{row['n_channels']},"
            f"{row['status']},"
            f"{row['note']}"
        )


if __name__ == "__main__":
    main()
