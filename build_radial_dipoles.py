import argparse
import csv

import numpy as np


def _load_orbitals(path):
    orbitals = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            atom = row["atom"].strip()
            orbital = {
                "atom": atom,
                "orbital": row["orbital"].strip(),
                "type": row["type"].strip().lower(),
                "n": int(row["n"]),
                "l": int(row["l"]),
                "occupation": float(row["occupation"]),
                "energy_Ha": float(row["energy_Ha"]),
            }
            orbitals.setdefault(atom, {})[orbital["orbital"]] = orbital
    return orbitals


def _load_radial_orbitals(path):
    radial = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            key = (row["atom"].strip(), row["orbital"].strip())
            radial.setdefault(key, {"r": [], "u": []})
            radial[key]["r"].append(float(row["r_bohr"]))
            radial[key]["u"].append(float(row["u"]))

    for key, values in radial.items():
        order = np.argsort(values["r"])
        values["r"] = np.asarray(values["r"], dtype=float)[order]
        values["u"] = np.asarray(values["u"], dtype=float)[order]
    return radial


def _normalize_u(r, u):
    norm = np.trapezoid(u * u, r)
    if norm <= 0.0:
        raise ValueError("Radial orbital norm must be positive.")
    return u / np.sqrt(norm)


def radial_integral(orbital_i, orbital_a, normalize=True):
    r_i = orbital_i["r"]
    u_i = orbital_i["u"]
    r_a = orbital_a["r"]
    u_a = orbital_a["u"]
    if normalize:
        u_i = _normalize_u(r_i, u_i)
        u_a = _normalize_u(r_a, u_a)

    r_min = max(float(r_i[0]), float(r_a[0]))
    r_max = min(float(r_i[-1]), float(r_a[-1]))
    if r_max <= r_min:
        raise ValueError("Radial orbital grids do not overlap.")

    grid = np.unique(np.concatenate([r_i[(r_i >= r_min) & (r_i <= r_max)], r_a[(r_a >= r_min) & (r_a <= r_max)]]))
    u_i_grid = np.interp(grid, r_i, u_i)
    u_a_grid = np.interp(grid, r_a, u_a)
    return float(np.trapezoid(u_i_grid * u_a_grid * grid, grid))


def angular_factor(l_i, l_a):
    if abs(l_a - l_i) != 1:
        return 0.0
    return max(l_i, l_a) / (2.0 * l_i + 1.0)


def build_radial_dipoles_from_orbitals(
    orbital_path="atomic_spectral_input.csv",
    radial_path="radial_orbitals.csv",
    normalize=True,
):
    orbitals_by_atom = _load_orbitals(orbital_path)
    radial = _load_radial_orbitals(radial_path)
    rows = []
    for atom, orbitals in sorted(orbitals_by_atom.items()):
        core_orbitals = [orb for orb in orbitals.values() if orb["type"] == "core"]
        virtual_orbitals = [orb for orb in orbitals.values() if orb["type"] == "virtual"]
        for occupied in core_orbitals:
            for virtual in virtual_orbitals:
                if virtual["energy_Ha"] <= occupied["energy_Ha"]:
                    continue
                factor = angular_factor(occupied["l"], virtual["l"])
                if factor == 0.0:
                    continue

                key_i = (atom, occupied["orbital"])
                key_a = (atom, virtual["orbital"])
                if key_i not in radial or key_a not in radial:
                    continue

                integral = radial_integral(radial[key_i], radial[key_a], normalize=normalize)
                d2 = factor * integral**2
                rows.append(
                    {
                        "atom": atom,
                        "from_orbital": occupied["orbital"],
                        "to_orbital": virtual["orbital"],
                        "d2": d2,
                        "radial_integral": integral,
                        "angular_factor": factor,
                    }
                )
    return rows


def write_radial_dipoles(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "from_orbital", "to_orbital", "d2"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "atom": row["atom"],
                    "from_orbital": row["from_orbital"],
                    "to_orbital": row["to_orbital"],
                    "d2": f"{row['d2']:.12f}",
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute radial_dipoles.csv from reduced radial orbitals u(r).")
    parser.add_argument("--orbitals", default="atomic_spectral_input.csv")
    parser.add_argument("--radial-orbitals", default="radial_orbitals.csv")
    parser.add_argument("--output", default="radial_dipoles.csv")
    parser.add_argument("--no-normalize", action="store_true")
    args = parser.parse_args(argv)

    rows = build_radial_dipoles_from_orbitals(
        args.orbitals,
        args.radial_orbitals,
        normalize=not args.no_normalize,
    )
    write_radial_dipoles(args.output, rows)

    print("atom,from_orbital,to_orbital,d2,radial_integral,angular_factor")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['from_orbital']},"
            f"{row['to_orbital']},"
            f"{row['d2']:.12f},"
            f"{row['radial_integral']:.12f},"
            f"{row['angular_factor']:.12f}"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
