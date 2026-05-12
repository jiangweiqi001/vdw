import csv
from eft_alpha import alpha0_from_osc, self_c6_from_osc, load_channels_csv


def main():
    data = load_channels_csv("atomic_channels.csv")
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

    with open("alpha_c6_table.csv", "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "alpha0_au", "C6_self_au", "n_channels"])
        writer.writeheader()
        writer.writerows(rows)

    print("atom,alpha0_au,C6_self_au,n_channels")
    for row in rows:
        print(f"{row['atom']},{row['alpha0_au']},{row['C6_self_au']},{row['n_channels']}")
    print("\nWrote alpha_c6_table.csv")


if __name__ == "__main__":
    main()
