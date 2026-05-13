import csv
import argparse
from eft_alpha import alpha0_from_osc, self_c6_from_osc, load_channels_csv


def build_alpha_rows(input_path="atomic_channels.csv"):
    data = load_channels_csv(input_path)
    rows = []
    for atom in sorted(data.keys()):
        delta = data[atom]["delta"]
        osc = data[atom]["osc"]
        alpha0 = alpha0_from_osc(delta, osc)
        c6_self = self_c6_from_osc(delta, osc)
        rows.append(
            {
                "atom": atom,
                "alpha0_au": f"{alpha0:.8f}",
                "C6_self_au": f"{c6_self:.8f}",
                "n_channels": len(delta),
            }
        )
    return rows


def write_alpha_rows(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "alpha0_au", "C6_self_au", "n_channels"])
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build alpha0 and self C6 table from atomic channels.")
    parser.add_argument("--input", default="atomic_channels.csv")
    parser.add_argument("--output", default="alpha_c6_table.csv")
    args = parser.parse_args(argv)

    rows = build_alpha_rows(args.input)
    write_alpha_rows(args.output, rows)

    print("atom,alpha0_au,C6_self_au,n_channels")
    for row in rows:
        print(f"{row['atom']},{row['alpha0_au']},{row['C6_self_au']},{row['n_channels']}")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
