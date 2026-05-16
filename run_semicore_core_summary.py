import argparse
import csv
import re
from pathlib import Path


PARTITION_LABELS = {
    "Ar": {
        "valence_partition": "3s;3p",
        "semicore_partition": "",
        "deep_core_partition": "1s;2s;2p",
        "basis": "aug-cc-pVQZ",
        "note": "noble-gas core/mixing diagnostic",
    },
    "Kr": {
        "valence_partition": "4s;4p",
        "semicore_partition": "",
        "deep_core_partition": "1s;2s;2p;3s;3p;3d",
        "basis": "aug-cc-pVQZ",
        "note": "noble-gas core/mixing diagnostic",
    },
    "Ca": {
        "valence_partition": "4s",
        "semicore_partition": "3s;3p",
        "deep_core_partition": "1s;2s;2p",
        "basis": "cc-pVQZ",
        "note": "Ca semicore partition benchmark; aug-cc-pVQZ unavailable in PySCF",
    },
    "Mg": {
        "valence_partition": "3s",
        "semicore_partition": "2s;2p",
        "deep_core_partition": "1s",
        "basis": "aug-cc-pVQZ",
        "note": "Mg semicore partition benchmark",
    },
}


def _read_component_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return {row["component"]: row for row in csv.DictReader(fp)}


def _relative_from_note(note):
    match = re.search(r"relative_(?:semicore|core)_contribution=([-+0-9.eE]+)", note or "")
    return float(match.group(1)) * 100.0 if match else ""


def summary_row(system, path, method="TDHF", nstates=200):
    rows = _read_component_rows(path)
    labels = PARTITION_LABELS[system]
    all_row = rows["all"]
    valence_row = rows["valence"]
    if "valence_plus_semicore" in rows:
        valence_plus_semicore = rows["valence_plus_semicore"]
        delta_row = rows["Delta_C6_semicore"]
        delta = float(delta_row["C6"])
        relative = _relative_from_note(delta_row["note"])
    else:
        valence_plus_semicore = rows["all"]
        delta_row = rows.get("Delta_C6_core")
        delta = float(delta_row["C6"]) if delta_row else float(all_row["C6"]) - float(valence_row["C6"])
        relative = _relative_from_note(delta_row["note"] if delta_row else "")
    return {
        "system": system,
        "method": method,
        "basis": labels["basis"],
        "nstates": nstates,
        "valence_partition": labels["valence_partition"],
        "semicore_partition": labels["semicore_partition"],
        "deep_core_partition": labels["deep_core_partition"],
        "C6_all": float(all_row["C6"]),
        "C6_valence_only": float(valence_row["C6"]),
        "C6_valence_plus_semicore": float(valence_plus_semicore["C6"]),
        "delta_C6_semicore": delta,
        "relative_delta_pct": relative,
        "note": labels["note"],
    }


def write_summary(path, rows):
    fieldnames = [
        "system",
        "method",
        "basis",
        "nstates",
        "valence_partition",
        "semicore_partition",
        "deep_core_partition",
        "C6_all",
        "C6_valence_only",
        "C6_valence_plus_semicore",
        "delta_C6_semicore",
        "relative_delta_pct",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_summary(output="results/semicore_core_correction_summary.csv"):
    inputs = {
        "Ar": "results/ar/ar_core_valence_decomposition.csv",
        "Kr": "results/kr/kr_core_valence_decomposition.csv",
        "Ca": "results/ca/ca_core_valence_decomposition.csv",
        "Mg": "results/mg/mg_core_valence_decomposition.csv",
    }
    rows = [summary_row(system, path) for system, path in inputs.items() if Path(path).exists()]
    write_summary(output, rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Summarize semicore/core C6 corrections.")
    parser.add_argument("--output", default="results/semicore_core_correction_summary.csv")
    args = parser.parse_args(argv)

    rows = build_summary(args.output)
    print("system,method,basis,nstates,C6_all,C6_valence_only,C6_valence_plus_semicore,delta_C6_semicore,relative_delta_pct,note")
    for row in rows:
        print(
            f"{row['system']},{row['method']},{row['basis']},{row['nstates']},"
            f"{row['C6_all']},{row['C6_valence_only']},{row['C6_valence_plus_semicore']},"
            f"{row['delta_C6_semicore']},{row['relative_delta_pct']},{row['note']}"
        )


if __name__ == "__main__":
    main()
