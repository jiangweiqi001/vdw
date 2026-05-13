import argparse
import csv

import numpy as np


def _value(row, names):
    for name in names:
        if name in row and str(row[name]).strip() != "":
            return float(row[name])
    raise KeyError(f"Missing any of columns: {', '.join(names)}")


def _normalize_converted(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault((row["atom"], row["orbital"]), []).append(row)

    normalized = []
    for key in sorted(grouped):
        group = sorted(grouped[key], key=lambda row: row["r_bohr"])
        r = np.asarray([row["r_bohr"] for row in group], dtype=float)
        u = np.asarray([row["u"] for row in group], dtype=float)
        norm = float(np.trapezoid(u * u, r))
        if norm <= 0.0:
            raise ValueError(f"Cannot normalize {key[0]} {key[1]}: non-positive norm.")
        scale = 1.0 / np.sqrt(norm)
        for row in group:
            normalized.append({**row, "u": row["u"] * scale})
    return normalized


def convert_rows(rows, input_kind="auto", normalize=False):
    converted = []
    for row in rows:
        atom = row["atom"].strip()
        orbital = row["orbital"].strip()
        r = _value(row, ["r_bohr", "r"])

        kind = input_kind.lower()
        if kind == "auto":
            kind = "u" if row.get("u", "").strip() else "r"

        if kind == "u":
            u = _value(row, ["u"])
        elif kind == "r":
            radial_R = _value(row, ["R", "R_of_r", "radial_R"])
            u = r * radial_R
        else:
            raise ValueError("input_kind must be one of: auto, R, u")

        converted.append({"atom": atom, "orbital": orbital, "r_bohr": r, "u": u})

    if normalize:
        converted = _normalize_converted(converted)
    return converted


def read_input(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_radial_orbitals(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "orbital", "r_bohr", "u"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "atom": row["atom"],
                    "orbital": row["orbital"],
                    "r_bohr": f"{row['r_bohr']:.12g}",
                    "u": f"{row['u']:.12g}",
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Convert atomic solver radial CSV output to atom,orbital,r_bohr,u. "
            "Use --input-kind R when the solver outputs R(r); this writes u(r)=rR(r)."
        )
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="radial_orbitals.csv")
    parser.add_argument("--input-kind", choices=["auto", "R", "u"], default="auto")
    parser.add_argument("--normalize", action="store_true")
    args = parser.parse_args(argv)

    rows = convert_rows(read_input(args.input), args.input_kind, normalize=args.normalize)
    write_radial_orbitals(args.output, rows)
    print(f"Wrote {args.output} with {len(rows)} radial grid rows.")


if __name__ == "__main__":
    main()
