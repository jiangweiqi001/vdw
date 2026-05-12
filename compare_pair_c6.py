import argparse
import csv


def _pair_key(a, b):
    return tuple(sorted((a.strip(), b.strip())))


def _read_model_pairs(path):
    pairs = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            pairs[_pair_key(row["A"], row["B"])] = float(row["C6_au"])
    return pairs


def _read_reference_pairs(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def compare_pair_tables(model_path="c6_table.csv", reference_path="reference_pair_c6.csv"):
    model_pairs = _read_model_pairs(model_path)
    rows = []
    for ref in _read_reference_pairs(reference_path):
        a = ref["A"].strip()
        b = ref["B"].strip()
        c6_ref = ref.get("C6_ref", "").strip()
        if not a or not b or not c6_ref:
            continue

        key = _pair_key(a, b)
        if key not in model_pairs:
            continue

        c6_model = model_pairs[key]
        c6_ref = float(c6_ref)
        rows.append(
            {
                "pair": f"{a}-{b}",
                "A": a,
                "B": b,
                "C6_model": c6_model,
                "C6_ref": c6_ref,
                "err_pct": 100.0 * (c6_model - c6_ref) / c6_ref,
                "source": ref.get("source", ""),
            }
        )
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare model heteronuclear C6 values to reference data.")
    parser.add_argument("--model", default="c6_table.csv")
    parser.add_argument("--reference", default="reference_pair_c6.csv")
    args = parser.parse_args(argv)

    print("A-B,C6_model,C6_ref,err_pct")
    for row in compare_pair_tables(args.model, args.reference):
        print(
            f"{row['pair']},"
            f"{row['C6_model']:.8f},"
            f"{row['C6_ref']:.8f},"
            f"{row['err_pct']:.6f}"
        )


if __name__ == "__main__":
    main()
