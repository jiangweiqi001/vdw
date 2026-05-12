import csv
import itertools
from eft_alpha import load_channels_csv, c6_from_alpha

def main():
    data = load_channels_csv("atomic_channels.csv")
    rows = []
    for a, b in itertools.combinations_with_replacement(sorted(data.keys()), 2):
        c6 = float(c6_from_alpha(
            data[a]["delta"],
            data[a]["osc"],
            data[b]["delta"],
            data[b]["osc"],
        ))
        rows.append({"A": a, "B": b, "C6_au": c6})

    with open("c6_table.csv", "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["A", "B", "C6_au"])
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(f"{row['A']}-{row['B']}: C6_au={row['C6_au']:.6f}")

if __name__ == "__main__":
    main()
