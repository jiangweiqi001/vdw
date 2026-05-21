import argparse
import csv
from pathlib import Path

from screened_pairwise_vdw import build_tail_rows


R_VALUES = [8, 10, 12, 15, 20, 30, 40]
DEFAULT_C6 = 647.60794451


def build_sensitivity_rows(c6=DEFAULT_C6):
    rows = []
    for epsilon in [1.0, 2.0, 4.0]:
        for row in build_tail_rows(c6, R_VALUES, model="dielectric", epsilon=epsilon):
            rows.append({**row, "screening_family": "model_dielectric", "note": "model W_v sensitivity; not ab initio screening"})
    for kappa in [0.05, 0.1, 0.2]:
        for row in build_tail_rows(c6, R_VALUES, model="yukawa", kappa=kappa):
            rows.append({**row, "screening_family": "model_yukawa", "note": "model W_v sensitivity; not ab initio screening"})
    return rows


def write_rows(path, rows):
    fieldnames = [
        "screening_family",
        "R_bohr",
        "model",
        "epsilon",
        "kappa",
        "C6",
        "E_Ha",
        "E_meV",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Mg q2 model W_v screening sensitivity.")
    parser.add_argument("--c6", type=float, default=DEFAULT_C6)
    parser.add_argument("--output", default="results/mg/mg_q2_model_screening_sensitivity.csv")
    args = parser.parse_args(argv)

    rows = build_sensitivity_rows(args.c6)
    write_rows(args.output, rows)
    print("screening_family,R_bohr,model,epsilon,kappa,C6,E_Ha,E_meV,note")
    for row in rows:
        print(
            f"{row['screening_family']},"
            f"{row['R_bohr']},"
            f"{row['model']},"
            f"{row['epsilon']},"
            f"{row['kappa']},"
            f"{row['C6']},"
            f"{row['E_Ha']:.12e},"
            f"{row['E_meV']:.12e},"
            f"{row['note']}"
        )


if __name__ == "__main__":
    main()
