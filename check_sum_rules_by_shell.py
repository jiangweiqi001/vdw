import argparse
import csv


def _load_core_occupations(path):
    occupations = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            if row["type"].strip().lower() != "core":
                continue
            occupations[(row["atom"].strip(), row["orbital"].strip())] = float(row["occupation"])
    return occupations


def _from_orbital(channel):
    if channel == "missing_core_continuum":
        return None
    if "_to_" not in channel:
        return None
    return channel.split("_to_", 1)[0]


def summarize_by_shell(
    channel_path="atomic_channels.csv",
    orbital_path="atomic_spectral_input.csv",
    warn_ratio=1.1,
):
    occupations = _load_core_occupations(orbital_path)
    sums = {key: 0.0 for key in occupations}
    with open(channel_path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            atom = row["atom"].strip()
            source_orbital = _from_orbital(row["channel"].strip())
            if source_orbital is None:
                continue
            key = (atom, source_orbital)
            if key not in sums:
                continue
            osc_text = row.get("osc", "").strip()
            if not osc_text:
                continue
            sums[key] += float(osc_text)

    rows = []
    for (atom, orbital), occupation in sorted(occupations.items()):
        sum_osc = sums[(atom, orbital)]
        ratio = sum_osc / occupation if occupation else 0.0
        rows.append(
            {
                "atom": atom,
                "from_orbital": orbital,
                "occupation": occupation,
                "sum_osc": sum_osc,
                "ratio_to_occupation": ratio,
                "warn": ratio > warn_ratio,
            }
        )
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Check oscillator strength sums by occupied shell.")
    parser.add_argument("--channels", default="atomic_channels.csv")
    parser.add_argument("--orbitals", default="atomic_spectral_input.csv")
    parser.add_argument("--warn-ratio", type=float, default=1.1)
    args = parser.parse_args(argv)

    print("atom,from_orbital,occupation,sum_osc,ratio_to_occupation,warn")
    for row in summarize_by_shell(args.channels, args.orbitals, args.warn_ratio):
        print(
            f"{row['atom']},"
            f"{row['from_orbital']},"
            f"{row['occupation']:.12g},"
            f"{row['sum_osc']:.12g},"
            f"{row['ratio_to_occupation']:.8f},"
            f"{str(row['warn']).lower()}"
        )


if __name__ == "__main__":
    main()
