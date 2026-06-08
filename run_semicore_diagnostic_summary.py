import argparse
import csv
from pathlib import Path


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def _to_float(value):
    if value in {"", None}:
        return None
    return float(value)


def load_validation(path, variant):
    row = read_rows(path)[0]
    return {
        "atom": row["atom"],
        "dimer": row["dimer"],
        "variant": variant,
        "C6_PSP": _to_float(row["C6_PSP"]),
        "C6_corrected": _to_float(row["C6_PSP_plus_sternheimer"]),
        "Delta_C6_core": _to_float(row["Delta_C6_core"]),
        "C6_all_e": _to_float(row.get("C6_all_e")),
        "closure_pct": _to_float(row.get("closure_pct")),
        "reference_primary": _to_float(row.get("C6_reference")),
        "reference_primary_error_pct": _to_float(row.get("reference_error_pct")),
        "go_no_go": row.get("go_no_go", ""),
    }


def load_references(path):
    refs = {}
    for row in read_rows(path):
        if row.get("C6_ref"):
            refs.setdefault(row["A"], []).append(row)
    return refs


def load_trk_rows(path):
    rows = {}
    if not path:
        return rows
    for row in read_rows(path):
        rows[row["atom"]] = {
            "raw_sum_osc": _to_float(row.get("raw_sum_osc")),
            "trk_scale": _to_float(row.get("trk_scale")),
        }
    return rows


def reference_error(c6_model, c6_ref):
    if c6_model is None or c6_ref in {None, 0.0}:
        return None
    return 100.0 * (c6_model - c6_ref) / c6_ref


def build_rows(raw_summaries, trk_summaries, references_path, trk_summary_path=None):
    references = load_references(references_path)
    trk_meta = load_trk_rows(trk_summary_path)
    variants = []
    for path in raw_summaries:
        variants.append(load_validation(path, "raw_finite_basis"))
    for path in trk_summaries:
        variants.append(load_validation(path, "trk_normalized"))

    rows = []
    for variant in variants:
        atom = variant["atom"]
        refs = references.get(atom, [])
        meta = trk_meta.get(atom, {}) if variant["variant"] == "trk_normalized" else {}
        for ref in refs:
            c6_ref = float(ref["C6_ref"])
            rows.append(
                {
                    "atom": atom,
                    "dimer": variant["dimer"],
                    "variant": variant["variant"],
                    "C6_PSP": variant["C6_PSP"],
                    "C6_corrected": variant["C6_corrected"],
                    "Delta_C6_core": variant["Delta_C6_core"],
                    "relative_delta_pct": 100.0 * variant["Delta_C6_core"] / variant["C6_PSP"]
                    if variant["C6_PSP"]
                    else "",
                    "C6_all_e": variant["C6_all_e"],
                    "closure_pct": variant["closure_pct"],
                    "raw_sum_osc": meta.get("raw_sum_osc", ""),
                    "trk_scale": meta.get("trk_scale", ""),
                    "C6_ref": c6_ref,
                    "reference_label": ref["reference_label"],
                    "source": ref["source"],
                    "corrected_error_pct": reference_error(variant["C6_corrected"], c6_ref),
                    "psp_error_pct": reference_error(variant["C6_PSP"], c6_ref),
                    "go_no_go": variant["go_no_go"],
                }
            )
    return rows


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "dimer",
        "variant",
        "C6_PSP",
        "C6_corrected",
        "Delta_C6_core",
        "relative_delta_pct",
        "C6_all_e",
        "closure_pct",
        "raw_sum_osc",
        "trk_scale",
        "C6_ref",
        "reference_label",
        "source",
        "psp_error_pct",
        "corrected_error_pct",
        "go_no_go",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Summarize raw/TRK semicore diagnostics against references.")
    parser.add_argument("--raw-summary", action="append", required=True)
    parser.add_argument("--trk-summary", action="append", default=[])
    parser.add_argument("--trk-summary-table")
    parser.add_argument("--references", default="reference_pair_c6_alternates.csv")
    parser.add_argument("--output", default="results/semicore_zn_cd_validation/diagnostic_summary.csv")
    args = parser.parse_args(argv)

    rows = build_rows(
        raw_summaries=args.raw_summary,
        trk_summaries=args.trk_summary,
        references_path=args.references,
        trk_summary_path=args.trk_summary_table,
    )
    write_rows(args.output, rows)
    print("atom,variant,C6_PSP,C6_corrected,Delta_C6_core,closure_pct,C6_ref,reference_label,corrected_error_pct")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['variant']},"
            f"{row['C6_PSP']},"
            f"{row['C6_corrected']},"
            f"{row['Delta_C6_core']},"
            f"{row['closure_pct']},"
            f"{row['C6_ref']},"
            f"{row['reference_label']},"
            f"{row['corrected_error_pct']}"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
