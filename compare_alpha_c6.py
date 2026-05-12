import argparse
import csv


def _read_by_atom(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return {row["atom"].strip(): row for row in csv.DictReader(fp) if row.get("atom", "").strip()}


def _percent_error(value, reference):
    return 100.0 * (value - reference) / reference


def compare_tables(eft_path="alpha_c6_table.csv", reference_path="reference_alpha_c6.csv"):
    eft_rows = _read_by_atom(eft_path)
    ref_rows = _read_by_atom(reference_path)
    rows = []

    for atom in sorted(ref_rows):
        ref = ref_rows[atom]
        if atom not in eft_rows:
            continue
        if not ref.get("alpha0_ref", "").strip() or not ref.get("C6_ref", "").strip():
            continue

        eft = eft_rows[atom]
        alpha0_eft = float(eft["alpha0_au"])
        c6_eft = float(eft["C6_self_au"])
        alpha0_ref = float(ref["alpha0_ref"])
        c6_ref = float(ref["C6_ref"])
        rows.append(
            {
                "atom": atom,
                "alpha0_eft": alpha0_eft,
                "alpha0_ref": alpha0_ref,
                "err_alpha_pct": _percent_error(alpha0_eft, alpha0_ref),
                "C6_eft": c6_eft,
                "C6_ref": c6_ref,
                "err_C6_pct": _percent_error(c6_eft, c6_ref),
                "source": ref.get("source", ""),
            }
        )
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare EFT alpha0/C6 against reference atomic data.")
    parser.add_argument("--eft", default="alpha_c6_table.csv")
    parser.add_argument("--reference", default="reference_alpha_c6.csv")
    args = parser.parse_args(argv)

    rows = compare_tables(args.eft, args.reference)
    fieldnames = [
        "atom",
        "alpha0_eft",
        "alpha0_ref",
        "err_alpha_pct",
        "C6_eft",
        "C6_ref",
        "err_C6_pct",
    ]
    print(",".join(fieldnames))
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['alpha0_eft']:.8f},"
            f"{row['alpha0_ref']:.8f},"
            f"{row['err_alpha_pct']:.6f},"
            f"{row['C6_eft']:.8f},"
            f"{row['C6_ref']:.8f},"
            f"{row['err_C6_pct']:.6f}"
        )


if __name__ == "__main__":
    main()
