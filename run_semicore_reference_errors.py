import argparse
import csv
from pathlib import Path


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def _pair_key(row):
    return tuple(sorted((row["A"], row["B"])))


def load_references(path):
    refs = {}
    for row in read_rows(path):
        if not row.get("C6_ref"):
            continue
        refs.setdefault(_pair_key(row), []).append(row)
    return refs


def model_rows_from_summary(summary_row):
    atom = summary_row["atom"]
    return [
        {
            "atom": atom,
            "dimer": summary_row["dimer"],
            "model": "psp",
            "C6_model": float(summary_row["C6_PSP"]),
        },
        {
            "atom": atom,
            "dimer": summary_row["dimer"],
            "model": "psp_plus_sternheimer",
            "C6_model": float(summary_row["C6_PSP_plus_sternheimer"]),
        },
    ]


def build_error_rows(summary_paths, reference_path):
    refs = load_references(reference_path)
    rows = []
    for summary_path in summary_paths:
        for summary in read_rows(summary_path):
            pair_refs = refs.get((summary["atom"], summary["atom"]), [])
            for model in model_rows_from_summary(summary):
                for ref in pair_refs:
                    c6_ref = float(ref["C6_ref"])
                    c6_model = float(model["C6_model"])
                    rows.append(
                        {
                            "atom": model["atom"],
                            "dimer": model["dimer"],
                            "model": model["model"],
                            "C6_model": c6_model,
                            "C6_ref": c6_ref,
                            "reference_label": ref["reference_label"],
                            "source": ref["source"],
                            "error_pct": 100.0 * (c6_model - c6_ref) / c6_ref,
                        }
                    )
    return rows


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "dimer",
        "model",
        "C6_model",
        "C6_ref",
        "reference_label",
        "source",
        "error_pct",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare semicore validation C6 against multiple references.")
    parser.add_argument("--summary", action="append", required=True)
    parser.add_argument("--references", default="reference_pair_c6_alternates.csv")
    parser.add_argument("--output", default="results/semicore_zn_cd_validation/multi_reference_errors.csv")
    args = parser.parse_args(argv)

    rows = build_error_rows(args.summary, args.references)
    write_rows(args.output, rows)
    print("atom,dimer,model,C6_model,C6_ref,reference_label,error_pct")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['dimer']},"
            f"{row['model']},"
            f"{row['C6_model']},"
            f"{row['C6_ref']},"
            f"{row['reference_label']},"
            f"{row['error_pct']}"
        )


if __name__ == "__main__":
    main()
