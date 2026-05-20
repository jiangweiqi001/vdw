import argparse
import csv
from pathlib import Path

from run_partition_decomposition import run_partition_decomposition


BASIS_LIST = ["cc-pVTZ", "cc-pVQZ", "cc-pVQZ-aug-local"]


def _read_decomposition(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return {row["component"]: row for row in csv.DictReader(fp)}


def run_basis_case(basis, nstates=200, output_root="results/ca/basis_sensitivity"):
    case_root = Path(output_root) / basis.lower()
    run_partition_decomposition(
        atom="Ca",
        basis=basis,
        nstates=nstates,
        partition_file="partition_definitions.csv",
        output_root=case_root,
    )
    rows = _read_decomposition(case_root / "ca" / "ca_core_valence_decomposition.csv")
    c6_all = float(rows["all"]["C6"])
    c6_val = float(rows["valence"]["C6"])
    c6_val_sc = float(rows["valence_plus_semicore"]["C6"])
    delta = float(rows["Delta_C6_semicore"]["C6"])
    return {
        "basis": basis,
        "nstates": nstates,
        "C6_all": c6_all,
        "C6_valence_only": c6_val,
        "C6_valence_plus_semicore": c6_val_sc,
        "Delta_C6_semicore": delta,
        "relative_semicore_contribution": delta / c6_all,
        "note": "local even-tempered diffuse augmentation; not official basis" if basis.endswith("-aug-local") else "",
    }


def write_summary(path, rows):
    fieldnames = [
        "basis",
        "nstates",
        "C6_all",
        "C6_valence_only",
        "C6_valence_plus_semicore",
        "Delta_C6_semicore",
        "relative_semicore_contribution",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_sensitivity(bases=None, nstates=200, output="results/ca/ca_partition_basis_sensitivity.csv"):
    rows = [run_basis_case(basis, nstates=nstates) for basis in (bases or BASIS_LIST)]
    write_summary(output, rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ca partition basis sensitivity.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST)
    parser.add_argument("--nstates", type=int, default=200)
    parser.add_argument("--output", default="results/ca/ca_partition_basis_sensitivity.csv")
    args = parser.parse_args(argv)

    rows = run_sensitivity(args.basis or BASIS_LIST, args.nstates, args.output)
    print("basis,nstates,C6_all,C6_valence_only,C6_valence_plus_semicore,Delta_C6_semicore,relative_semicore_contribution,note")
    for row in rows:
        print(
            f"{row['basis']},{row['nstates']},{row['C6_all']},{row['C6_valence_only']},"
            f"{row['C6_valence_plus_semicore']},{row['Delta_C6_semicore']},"
            f"{row['relative_semicore_contribution']},{row['note']}"
        )


if __name__ == "__main__":
    main()
