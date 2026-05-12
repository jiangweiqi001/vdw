import argparse
import csv

from eft_alpha import alpha0_from_osc, load_channels_csv, self_c6_from_osc


def _core_electron_counts(orbital_path):
    counts = {}
    with open(orbital_path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            if row["type"].strip().lower() != "core":
                continue
            atom = row["atom"].strip()
            counts[atom] = counts.get(atom, 0.0) + float(row["occupation"])
    return counts


def summarize_sum_rules(
    orbital_path="atomic_spectral_input.csv",
    channel_path="atomic_channels.csv",
):
    n_core_by_atom = _core_electron_counts(orbital_path)
    channels = load_channels_csv(channel_path)
    rows = []
    for atom in sorted(set(n_core_by_atom) | set(channels)):
        n_core = n_core_by_atom.get(atom, 0.0)
        if atom in channels:
            delta = channels[atom]["delta"]
            osc = channels[atom]["osc"]
            sum_osc = float(osc.sum())
            alpha0 = float(alpha0_from_osc(delta, osc))
            c6_self = float(self_c6_from_osc(delta, osc))
        else:
            sum_osc = 0.0
            alpha0 = 0.0
            c6_self = 0.0

        rows.append(
            {
                "atom": atom,
                "N_core": n_core,
                "sum_osc": sum_osc,
                "sum_osc_over_N_core": sum_osc / n_core if n_core else 0.0,
                "alpha0": alpha0,
                "C6_self": c6_self,
            }
        )
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Check core oscillator-strength sum rules.")
    parser.add_argument("--orbitals", default="atomic_spectral_input.csv")
    parser.add_argument("--channels", default="atomic_channels.csv")
    args = parser.parse_args(argv)

    print("atom,N_core,sum_osc,sum_osc/N_core,alpha0,C6_self")
    for row in summarize_sum_rules(args.orbitals, args.channels):
        print(
            f"{row['atom']},"
            f"{row['N_core']:.8f},"
            f"{row['sum_osc']:.8f},"
            f"{row['sum_osc_over_N_core']:.8f},"
            f"{row['alpha0']:.8f},"
            f"{row['C6_self']:.8f}"
        )


if __name__ == "__main__":
    main()
