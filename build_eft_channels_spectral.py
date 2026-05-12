import argparse
import csv
import json

import numpy as np


FIELDNAMES = ["atom", "channel", "delta_Ha", "d2", "osc", "is_core"]


def build_spectral_channels(atom, occupied_core_orbitals, virtual_orbitals):
    """Convert frozen-core orbital transitions into dipole-channel CSV rows.

    Orbital dictionaries use energies in Hartree. Occupied orbitals marked with
    `is_frozen_core=False` are skipped to avoid double counting explicit valence.
    """
    rows = []
    for occupied in occupied_core_orbitals:
        if not occupied.get("is_frozen_core", True):
            continue

        i_label = occupied["label"]
        i_energy = float(occupied["energy_Ha"])
        occupation = float(occupied.get("occupation", 1.0))

        for virtual in virtual_orbitals:
            a_label = virtual["label"]
            a_energy = float(virtual["energy_Ha"])
            delta = a_energy - i_energy
            if delta <= 0.0:
                raise ValueError(f"Transition {i_label}_to_{a_label} has non-positive delta.")

            dipole = np.asarray(virtual["dipole"], dtype=float)
            if dipole.shape != (3,):
                raise ValueError(f"Transition {i_label}_to_{a_label} dipole must have 3 components.")

            d2 = float(np.dot(dipole, dipole))
            osc = (2.0 / 3.0) * delta * occupation * d2
            rows.append(
                {
                    "atom": atom,
                    "channel": f"{i_label}_to_{a_label}",
                    "delta_Ha": f"{delta:.12f}",
                    "d2": f"{d2:.12f}",
                    "osc": f"{osc:.12f}",
                    "is_core": "true",
                }
            )
    return rows


def build_spectral_channels_from_json(path):
    """Read spectral channel specs from JSON and return atomic channel rows."""
    with open(path, encoding="utf-8") as fp:
        spec = json.load(fp)

    atoms = spec.get("atoms", [spec])
    rows = []
    for atom_spec in atoms:
        rows.extend(
            build_spectral_channels(
                atom=atom_spec["atom"],
                occupied_core_orbitals=atom_spec.get("occupied_core_orbitals", []),
                virtual_orbitals=atom_spec.get("virtual_orbitals", []),
            )
        )
    return rows


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


def _load_dipoles(path):
    dipoles = []
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            dipoles.append(
                {
                    "atom": row["atom"].strip(),
                    "from_orbital": row["from_orbital"].strip(),
                    "to_orbital": row["to_orbital"].strip(),
                    "d2": float(row["d2"]),
                }
            )
    return dipoles


def _load_residual_deltas(path):
    deltas = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            deltas[row["atom"].strip()] = float(row["delta_missing_Ha"])
    return deltas


def _core_electron_counts(orbitals_by_atom):
    counts = {}
    for atom, orbitals in orbitals_by_atom.items():
        counts[atom] = sum(
            orbital["occupation"]
            for orbital in orbitals.values()
            if orbital["type"] == "core"
        )
    return counts


def add_residual_oscillators(rows, orbitals_by_atom, residual_path):
    residual_deltas = _load_residual_deltas(residual_path)
    n_core_by_atom = _core_electron_counts(orbitals_by_atom)
    osc_by_atom = {}
    for row in rows:
        osc_by_atom[row["atom"]] = osc_by_atom.get(row["atom"], 0.0) + float(row["osc"])

    rows = list(rows)
    for atom in sorted(n_core_by_atom):
        missing = n_core_by_atom[atom] - osc_by_atom.get(atom, 0.0)
        if missing <= 0.0:
            continue
        if atom not in residual_deltas:
            raise ValueError(f"Missing residual oscillator delta for {atom}.")
        rows.append(
            {
                "atom": atom,
                "channel": "missing_core_continuum",
                "delta_Ha": f"{residual_deltas[atom]:.12f}",
                "d2": "",
                "osc": f"{missing:.12f}",
                "is_core": "true",
            }
        )
    return rows


def build_spectral_channels_from_csv(
    orbital_path="atomic_spectral_input.csv",
    dipole_path="radial_dipoles.csv",
    add_residual_oscillator=False,
    residual_path="residual_oscillators.csv",
):
    """Build spectral channels from discrete orbital energies and radial d2 inputs."""
    orbitals_by_atom = _load_orbitals(orbital_path)
    rows = []
    for dipole in _load_dipoles(dipole_path):
        atom = dipole["atom"]
        atom_orbitals = orbitals_by_atom.get(atom, {})
        occupied = atom_orbitals.get(dipole["from_orbital"])
        virtual = atom_orbitals.get(dipole["to_orbital"])
        if occupied is None or virtual is None:
            continue
        if occupied["type"] != "core" or virtual["type"] != "virtual":
            continue
        if abs(occupied["l"] - virtual["l"]) != 1:
            continue

        delta = virtual["energy_Ha"] - occupied["energy_Ha"]
        if delta <= 0.0:
            raise ValueError(
                f"Transition {atom}:{occupied['orbital']}_to_{virtual['orbital']} has non-positive delta."
            )

        d2 = dipole["d2"]
        osc = (2.0 / 3.0) * delta * occupied["occupation"] * d2
        rows.append(
            {
                "atom": atom,
                "channel": f"{occupied['orbital']}_to_{virtual['orbital']}",
                "delta_Ha": f"{delta:.12f}",
                "d2": f"{d2:.12f}",
                "osc": f"{osc:.12f}",
                "is_core": "true",
            }
        )
    if add_residual_oscillator:
        rows = add_residual_oscillators(rows, orbitals_by_atom, residual_path)
    return rows


def write_channels(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build EFT dipole channels from spectral orbital input.")
    parser.add_argument("--input-json", help="JSON file with atom, occupied_core_orbitals, virtual_orbitals.")
    parser.add_argument("--orbitals", default="atomic_spectral_input.csv")
    parser.add_argument("--dipoles", default="radial_dipoles.csv")
    parser.add_argument("--residuals", default="residual_oscillators.csv")
    parser.add_argument("--add-residual-oscillator", action="store_true")
    parser.add_argument("--output", default="atomic_channels.csv")
    args = parser.parse_args(argv)

    if args.input_json:
        rows = build_spectral_channels_from_json(args.input_json)
    else:
        rows = build_spectral_channels_from_csv(
            args.orbitals,
            args.dipoles,
            add_residual_oscillator=args.add_residual_oscillator,
            residual_path=args.residuals,
        )
    write_channels(args.output, rows)
    print(f"Wrote {args.output} with {len(rows)} spectral dipole-channel rows.")


if __name__ == "__main__":
    main()
