import argparse
import csv
import importlib.util
from pathlib import Path


REFERENCE_C6 = 64.3


def percent_error(value, reference=REFERENCE_C6):
    return 100.0 * (float(value) - float(reference)) / float(reference)


def method_row(method, c6, reference=REFERENCE_C6, note=""):
    if c6 is None:
        return {"method": method, "C6_ArAr": "", "error_pct": "", "note": note}
    return {
        "method": method,
        "C6_ArAr": float(c6),
        "error_pct": percent_error(c6, reference),
        "note": note,
    }


def optional_package_note(package_names):
    for package in package_names:
        if importlib.util.find_spec(package) is not None:
            return "package_available_manual_parameter_extraction_needed"
    return "not_available"


def build_method_rows():
    return [
        method_row("reference", REFERENCE_C6, note="reference"),
        method_row("calibrated", REFERENCE_C6, note="fitted/control"),
        method_row("EFT-MO", 76.13748508, note="independent-particle aug-cc-pVQZ"),
        method_row("EFT-TDHF", 60.73027908, note="TDHF aug-cc-pVQZ nstates=200"),
        method_row("D4", None, note=optional_package_note(["dftd4"])),
        method_row("MBD", None, note=optional_package_note(["pymbd", "mbd"])),
    ]


def write_rows(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["method", "C6_ArAr", "error_pct", "note"])
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare long-range Ar-Ar C6 methods.")
    parser.add_argument("--output", default="results/ar/ar_c6_method_comparison.csv")
    args = parser.parse_args(argv)

    rows = build_method_rows()
    write_rows(args.output, rows)
    print("method,C6_ArAr,error_pct,note")
    for row in rows:
        print(f"{row['method']},{row['C6_ArAr']},{row['error_pct']},{row['note']}")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
