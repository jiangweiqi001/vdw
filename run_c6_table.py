import csv
import itertools
import argparse
from eft_alpha import load_channels_csv, c6_from_alpha

def build_c6_rows(input_path="atomic_channels.csv"):
    data = load_channels_csv(input_path)
    rows = []
    for a, b in itertools.combinations_with_replacement(sorted(data.keys()), 2):
        c6 = float(c6_from_alpha(
            data[a]["delta"],
            data[a]["osc"],
            data[b]["delta"],
            data[b]["osc"],
        ))
        rows.append({"A": a, "B": b, "C6_au": c6})
    return rows


def write_c6_rows(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["A", "B", "C6_au"])
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build pair C6 table from atomic channels.")
    parser.add_argument("--input", default="atomic_channels.csv")
    parser.add_argument("--output", default="c6_table.csv")
    args = parser.parse_args(argv)

    rows = build_c6_rows(args.input)
    write_c6_rows(args.output, rows)

    for row in rows:
        print(f"{row['A']}-{row['B']}: C6_au={row['C6_au']:.6f}")

if __name__ == "__main__":
    main()
