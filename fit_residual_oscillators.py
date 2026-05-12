import argparse
import csv
import math

from build_eft_channels_spectral import build_spectral_channels_from_csv


def _core_electron_counts(orbital_path):
    counts = {}
    with open(orbital_path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            if row["type"].strip().lower() == "core":
                atom = row["atom"].strip()
                counts[atom] = counts.get(atom, 0.0) + float(row["occupation"])
    return counts


def _reference_alpha0(reference_path):
    refs = {}
    with open(reference_path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            alpha0 = row.get("alpha0_ref", "").strip()
            if alpha0:
                refs[row["atom"].strip()] = float(alpha0)
    return refs


def _discrete_response_by_atom(orbital_path, dipole_path):
    rows = build_spectral_channels_from_csv(
        orbital_path,
        dipole_path,
        add_residual_oscillator=False,
    )
    response = {}
    for row in rows:
        atom = row["atom"]
        delta = float(row["delta_Ha"])
        osc = float(row["osc"])
        atom_response = response.setdefault(atom, {"sum_osc": 0.0, "alpha0": 0.0})
        atom_response["sum_osc"] += osc
        atom_response["alpha0"] += osc / delta**2
    return response


def fit_residual_deltas(
    orbital_path="atomic_spectral_input.csv",
    dipole_path="radial_dipoles.csv",
    reference_path="reference_alpha_c6.csv",
):
    n_core = _core_electron_counts(orbital_path)
    alpha_ref = _reference_alpha0(reference_path)
    response = _discrete_response_by_atom(orbital_path, dipole_path)
    rows = []
    for atom in sorted(n_core):
        if atom not in alpha_ref:
            continue

        discrete = response.get(atom, {"sum_osc": 0.0, "alpha0": 0.0})
        f_missing = n_core[atom] - discrete["sum_osc"]
        alpha_missing = alpha_ref[atom] - discrete["alpha0"]
        if f_missing <= 0.0:
            continue
        if alpha_missing <= 0.0:
            raise ValueError(
                f"Cannot fit residual delta for {atom}: discrete alpha0 exceeds reference."
            )

        rows.append(
            {
                "atom": atom,
                "delta_missing_Ha": math.sqrt(f_missing / alpha_missing),
                "f_missing": f_missing,
                "alpha_discrete": discrete["alpha0"],
                "alpha_ref": alpha_ref[atom],
            }
        )
    return rows


def write_residual_deltas(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "delta_missing_Ha"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "atom": row["atom"],
                    "delta_missing_Ha": f"{row['delta_missing_Ha']:.12f}",
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Diagnostic/empirical helper: fit spectral residual oscillator energies "
            "to reference alpha0. Do not treat fitted residuals as EFT predictions."
        )
    )
    parser.add_argument("--orbitals", default="atomic_spectral_input.csv")
    parser.add_argument("--dipoles", default="radial_dipoles.csv")
    parser.add_argument("--reference", default="reference_alpha_c6.csv")
    parser.add_argument("--output", default="residual_oscillators.csv")
    args = parser.parse_args(argv)

    rows = fit_residual_deltas(args.orbitals, args.dipoles, args.reference)
    write_residual_deltas(args.output, rows)

    print("atom,delta_missing_Ha,f_missing,alpha_discrete,alpha_ref")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['delta_missing_Ha']:.12f},"
            f"{row['f_missing']:.12f},"
            f"{row['alpha_discrete']:.12f},"
            f"{row['alpha_ref']:.12f}"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
