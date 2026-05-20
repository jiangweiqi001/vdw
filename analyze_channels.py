import argparse
import csv

import numpy as np

from eft_alpha import c6_from_alpha


def _load_channel_rows(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            delta = float(row["delta_Ha"])
            osc_text = row.get("osc", "").strip()
            d2_text = row.get("d2", "").strip()
            if osc_text:
                osc = float(osc_text)
            elif d2_text:
                osc = (2.0 / 3.0) * delta * float(d2_text)
            else:
                raise ValueError(f"Channel {row['atom']} {row['channel']} has no osc or d2.")
            rows.append(
                {
                    "atom": row["atom"].strip(),
                    "channel": row["channel"].strip(),
                    "delta_Ha": delta,
                    "osc": osc,
                    "is_residual": row["channel"].strip() == "missing_core_continuum",
                }
            )
    return rows


def _pair_c6_contributions(deltas, oscs):
    n = len(deltas)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            c6 = float(c6_from_alpha([deltas[i]], [oscs[i]], [deltas[j]], [oscs[j]]))
            matrix[i, j] = c6
            matrix[j, i] = c6
    total = 0.0
    channel_totals = np.zeros(n, dtype=float)
    for i in range(n):
        total += matrix[i, i]
        channel_totals[i] += matrix[i, i]
        for j in range(i + 1, n):
            contribution = 2.0 * matrix[i, j]
            total += contribution
            channel_totals[i] += matrix[i, j]
            channel_totals[j] += matrix[i, j]
    return matrix.diagonal(), channel_totals, total


def analyze_channels(path="atomic_channels.csv"):
    input_rows = _load_channel_rows(path)
    output_rows = []
    for atom in sorted({row["atom"] for row in input_rows}):
        atom_rows = [row for row in input_rows if row["atom"] == atom]
        deltas = np.asarray([row["delta_Ha"] for row in atom_rows], dtype=float)
        oscs = np.asarray([row["osc"] for row in atom_rows], dtype=float)

        alpha_contrib = oscs / deltas**2
        alpha_total = float(alpha_contrib.sum())
        single_c6, cross_c6, c6_total = _pair_c6_contributions(deltas, oscs)

        for row, alpha_i, single_i, cross_i in zip(atom_rows, alpha_contrib, single_c6, cross_c6):
            output_rows.append(
                {
                    "atom": atom,
                    "channel": row["channel"],
                    "delta_Ha": row["delta_Ha"],
                    "osc": row["osc"],
                    "alpha0_contribution": float(alpha_i),
                    "alpha0_fraction": float(alpha_i / alpha_total) if alpha_total else 0.0,
                    "single_channel_c6": float(single_i),
                    "single_channel_c6_fraction": float(single_i / c6_total) if c6_total else 0.0,
                    "cross_inclusive_c6": float(cross_i),
                    "cross_inclusive_c6_fraction": float(cross_i / c6_total) if c6_total else 0.0,
                    "is_residual": row["is_residual"],
                }
            )

    return sorted(output_rows, key=lambda row: (row["atom"], -row["alpha0_contribution"]))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Analyze per-channel alpha0 and C6 contributions.")
    parser.add_argument("--input", default="atomic_channels.csv")
    args = parser.parse_args(argv)

    fieldnames = [
        "atom",
        "channel",
        "delta_Ha",
        "osc",
        "alpha0_contribution",
        "alpha0_fraction",
        "single_channel_c6",
        "single_channel_c6_fraction",
        "cross_inclusive_c6",
        "cross_inclusive_c6_fraction",
        "is_residual",
    ]
    print(",".join(fieldnames))
    for row in analyze_channels(args.input):
        print(
            f"{row['atom']},"
            f"{row['channel']},"
            f"{row['delta_Ha']:.12f},"
            f"{row['osc']:.12f},"
            f"{row['alpha0_contribution']:.12f},"
            f"{row['alpha0_fraction']:.8f},"
            f"{row['single_channel_c6']:.12f},"
            f"{row['single_channel_c6_fraction']:.8f},"
            f"{row['cross_inclusive_c6']:.12f},"
            f"{row['cross_inclusive_c6_fraction']:.8f},"
            f"{str(row['is_residual']).lower()}"
        )


if __name__ == "__main__":
    main()
