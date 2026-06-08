import argparse
import csv
from pathlib import Path

import numpy as np


SOURCE = "CORE_STERNHEIMER_FINITE_BASIS"

DEFAULT_TARGETS = {
    "Sr": {"basis": "def2-TZVPP", "shells": {"4s", "4p"}},
    "Zn": {"basis": "def2-TZVPP", "shells": {"3d"}},
    "Cd": {"basis": "def2-TZVPP", "shells": {"4d"}},
}


def _selected_indices(mo_occ, mo_to_shell, selected_shells):
    selected_shells = set(selected_shells)
    return [
        idx
        for idx, occ in enumerate(np.asarray(mo_occ, dtype=float))
        if occ > 1e-12 and mo_to_shell.get(int(idx)) in selected_shells
    ]


def sternheimer_alpha_iw_from_arrays(xi, mo_energy, mo_occ, dipole_mo, mo_to_shell, selected_shells):
    """Compute frozen-core isotropic alpha(i xi) from a finite-basis Sternheimer solve.

    In the canonical MO basis the projected Sternheimer equation for one occupied
    core orbital c and Cartesian dipole component mu is diagonal in virtual
    orbitals a.  The summed +/- imaginary-frequency response is

        delta_c^mu(a; i xi) = -2 Delta_ac d_ac^mu / (Delta_ac^2 + xi^2).

    Contracting -<c|r_mu|delta_c^mu>/3 over x/y/z gives the usual finite-basis
    oscillator expression, while keeping the response-equation interface needed
    for continuum/discretized-continuum basis extensions.
    """
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    mo_energy = np.asarray(mo_energy, dtype=float)
    mo_occ = np.asarray(mo_occ, dtype=float)
    dipole_mo = np.asarray(dipole_mo, dtype=float)

    occupied = _selected_indices(mo_occ, mo_to_shell, selected_shells)
    virtual = np.where(mo_occ <= 1e-12)[0]
    alpha = np.zeros_like(xi, dtype=float)

    for i in occupied:
        occ = float(mo_occ[i])
        for a in virtual:
            delta = float(mo_energy[a] - mo_energy[i])
            if delta <= 0.0:
                continue
            d2 = float(np.sum(np.abs(dipole_mo[:, i, a]) ** 2))
            alpha += (2.0 / 3.0) * occ * delta * d2 / (delta**2 + xi**2)
    return alpha


def sternheimer_channels_from_arrays(atom, mo_energy, mo_occ, dipole_mo, mo_to_shell, selected_shells, min_osc=1e-14):
    """Return equivalent oscillator rows for the finite-basis Sternheimer response."""
    mo_energy = np.asarray(mo_energy, dtype=float)
    mo_occ = np.asarray(mo_occ, dtype=float)
    dipole_mo = np.asarray(dipole_mo, dtype=float)
    occupied = _selected_indices(mo_occ, mo_to_shell, selected_shells)
    virtual = np.where(mo_occ <= 1e-12)[0]
    rows = []

    for i in occupied:
        shell = mo_to_shell[int(i)]
        occ = float(mo_occ[i])
        for a in virtual:
            delta = float(mo_energy[a] - mo_energy[i])
            if delta <= 0.0:
                continue
            d2 = float(np.sum(np.abs(dipole_mo[:, i, a]) ** 2))
            osc = (2.0 / 3.0) * occ * delta * d2
            if abs(osc) <= min_osc:
                continue
            rows.append(
                {
                    "atom": atom,
                    "channel": f"sternheimer_{shell}_occ_{int(i):03d}_to_resp_{int(a):03d}",
                    "delta_Ha": delta,
                    "osc": osc,
                    "is_core": "true",
                    "source": SOURCE,
                    "from_shell": shell,
                    "occ_mo": int(i),
                    "response_basis_index": int(a),
                    "occ_energy_Ha": float(mo_energy[i]),
                    "response_energy_Ha": float(mo_energy[a]),
                    "d2_mo": d2,
                }
            )
    return rows


def compute_core_sternheimer_channels(atom, basis, selected_shells):
    from pyscf.scf import atom_hf
    from pyscf_export_ar_radials import grouped_shells, make_atom_molecule

    mol = make_atom_molecule(atom, basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} atom HF did not converge for {basis}.")

    dipole_ao = mol.intor("int1e_r", comp=3)
    mo_coeff = mf.mo_coeff
    dipole_mo = np.asarray([mo_coeff.T @ dipole_ao[mu] @ mo_coeff for mu in range(3)])
    mo_to_shell = {}
    for shell in grouped_shells(mol, mf.mo_energy, mf.mo_occ, mf.mo_coeff):
        for column in shell["columns"]:
            mo_to_shell[int(column)] = shell["label"]
    return sternheimer_channels_from_arrays(
        atom=atom,
        mo_energy=mf.mo_energy,
        mo_occ=mf.mo_occ,
        dipole_mo=dipole_mo,
        mo_to_shell=mo_to_shell,
        selected_shells=set(selected_shells),
    )


def write_channels(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "atom",
        "channel",
        "delta_Ha",
        "osc",
        "is_core",
        "source",
        "from_shell",
        "occ_mo",
        "response_basis_index",
        "occ_energy_Ha",
        "response_energy_Ha",
        "d2_mo",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "atom": row["atom"],
                    "channel": row["channel"],
                    "delta_Ha": f"{row['delta_Ha']:.12f}",
                    "osc": f"{row['osc']:.12f}",
                    "is_core": row["is_core"],
                    "source": row["source"],
                    "from_shell": row["from_shell"],
                    "occ_mo": row["occ_mo"],
                    "response_basis_index": row["response_basis_index"],
                    "occ_energy_Ha": f"{row['occ_energy_Ha']:.12f}",
                    "response_energy_Ha": f"{row['response_energy_Ha']:.12f}",
                    "d2_mo": f"{row['d2_mo']:.12f}",
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute frozen-core l=1 Sternheimer response channels.")
    parser.add_argument("--atom", action="append")
    parser.add_argument("--basis", action="append")
    parser.add_argument("--shell", action="append")
    parser.add_argument("--output", default="results/core_sternheimer_channels.csv")
    args = parser.parse_args(argv)

    all_rows = []
    atoms = args.atom or sorted(DEFAULT_TARGETS)
    basis_values = args.basis or [DEFAULT_TARGETS[atom]["basis"] for atom in atoms]
    if len(basis_values) == 1 and len(atoms) > 1:
        basis_values = basis_values * len(atoms)

    for atom, basis in zip(atoms, basis_values):
        shells = set(args.shell or DEFAULT_TARGETS.get(atom, {}).get("shells", set()))
        if not shells:
            raise ValueError(f"No target shells supplied for {atom}.")
        all_rows.extend(compute_core_sternheimer_channels(atom, basis, shells))

    write_channels(args.output, all_rows)
    print("atom,from_shell,n_channels,sum_osc,output")
    grouped = {}
    for row in all_rows:
        grouped.setdefault((row["atom"], row["from_shell"]), {"count": 0, "sum_osc": 0.0})
        grouped[(row["atom"], row["from_shell"])]["count"] += 1
        grouped[(row["atom"], row["from_shell"])]["sum_osc"] += row["osc"]
    for (atom, shell), summary in sorted(grouped.items()):
        print(f"{atom},{shell},{summary['count']},{summary['sum_osc']:.12f},{args.output}")


if __name__ == "__main__":
    main()
