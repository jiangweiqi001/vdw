import argparse
import csv

import numpy as np
from pyscf import gto
from pyscf.gto.basis import BasisNotFoundError
from pyscf.scf import atom_hf
from pyscf.scf.atom_hf import _angular_momentum_for_each_ao


L_LETTER = {0: "s", 1: "p", 2: "d", 3: "f"}


def shell_label(l_value, shell_index):
    return f"{l_value + 1 + shell_index}{L_LETTER[l_value]}"


def shell_occupation(l_value, component_occ):
    return component_occ * (2 * l_value + 1)


def infer_mo_l_values(mol, mo_coeff):
    ao_l = _angular_momentum_for_each_ao(mol)
    l_values = []
    for column in range(mo_coeff.shape[1]):
        weights = [
            np.linalg.norm(mo_coeff[ao_l == l_value, column])
            for l_value in range(max(L_LETTER) + 1)
        ]
        l_values.append(int(np.argmax(weights)))
    return np.asarray(l_values, dtype=int)


def grouped_shells(mol, mo_energy, mo_occ, mo_coeff):
    l_by_mo = infer_mo_l_values(mol, mo_coeff)
    shells = []
    for l_value in sorted(set(l_by_mo)):
        shell_columns = np.where(l_by_mo == l_value)[0]
        shell_columns = sorted(shell_columns, key=lambda column: mo_energy[column])
        used = set()
        shell_index = 0
        for column in shell_columns:
            if column in used:
                continue
            same_energy = [
                other
                for other in shell_columns
                if other not in used and abs(mo_energy[other] - mo_energy[column]) < 1e-8
            ]
            used.update(same_energy)
            shells.append(
                {
                    "l": l_value,
                    "label": shell_label(l_value, shell_index),
                    "columns": same_energy,
                    "energy_Ha": float(mo_energy[column]),
                    "occupation": float(shell_occupation(l_value, mo_occ[column])),
                    "type": "core" if mo_occ[column] > 0.0 else "virtual",
                    "n": l_value + 1 + shell_index,
                }
            )
            shell_index += 1
    return shells


def double_augmented_basis(symbol="Ar", parent_basis="aug-cc-pVTZ"):
    """Construct an even-tempered second augmentation shell from aug-cc-pVTZ.

    This is a local fallback for environments without Basis Set Exchange data.
    It adds one more diffuse primitive per angular momentum using the ratio of
    the two most diffuse existing exponents for that angular momentum.
    """
    basis = gto.basis.load(parent_basis, symbol)
    augmented = [shell.copy() for shell in basis]
    exponents_by_l = {}
    for shell in basis:
        l_value = shell[0]
        for primitive in shell[1:]:
            exponents_by_l.setdefault(l_value, []).append(float(primitive[0]))

    for l_value, exponents in sorted(exponents_by_l.items()):
        unique = sorted(set(exponents))
        if len(unique) < 2:
            continue
        most_diffuse, second_diffuse = unique[0], unique[1]
        new_exponent = most_diffuse * most_diffuse / second_diffuse
        augmented.append([l_value, [new_exponent, 1.0]])
    return augmented


def make_atom_molecule(atom, basis):
    try:
        return gto.M(atom=f"{atom} 0 0 0", basis=basis, spin=0, charge=0, cart=False, verbose=0)
    except BasisNotFoundError:
        if atom == "Ar" and basis in {"d-aug-cc-pVTZ", "d-aug-cc-pVTZ-local"}:
            return gto.M(
                atom="Ar 0 0 0",
                basis={"Ar": double_augmented_basis("Ar", "aug-cc-pVTZ")},
                spin=0,
                charge=0,
                cart=False,
                verbose=0,
            )
        if basis.endswith("-aug-local"):
            parent_basis = basis[: -len("-aug-local")]
            return gto.M(
                atom=f"{atom} 0 0 0",
                basis={atom: double_augmented_basis(atom, parent_basis)},
                spin=0,
                charge=0,
                cart=False,
                verbose=0,
            )
        raise


def make_ar_molecule(basis):
    return make_atom_molecule("Ar", basis)


def export_ar_radials(
    solver_output_path="atomic_solver_output.csv",
    spectral_input_path="atomic_spectral_input.csv",
    basis="cc-pvtz",
    r_min=1e-5,
    r_max=40.0,
    n_grid=2400,
):
    mol = make_ar_molecule(basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("PySCF spherically averaged Ar HF did not converge.")

    shells = grouped_shells(mol, mf.mo_energy, mf.mo_occ, mf.mo_coeff)
    r_grid = np.geomspace(r_min, r_max, n_grid)
    coords = np.column_stack([np.zeros_like(r_grid), np.zeros_like(r_grid), r_grid])
    ao_values = mol.eval_gto("GTOval_sph", coords)
    mo_values = ao_values @ mf.mo_coeff

    with open(solver_output_path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "orbital", "r_bohr", "R"])
        writer.writeheader()
        for shell in shells:
            l_value = shell["l"]
            y_l0_at_z = np.sqrt((2 * l_value + 1) / (4.0 * np.pi))
            # Along the z-axis only the m=0 real spherical component is non-zero.
            column = max(shell["columns"], key=lambda idx: np.max(np.abs(mo_values[:, idx])))
            radial_R = mo_values[:, column] / y_l0_at_z
            for r_value, R_value in zip(r_grid, radial_R):
                writer.writerow(
                    {
                        "atom": "Ar",
                        "orbital": shell["label"],
                        "r_bohr": f"{r_value:.12g}",
                        "R": f"{R_value:.12g}",
                    }
                )

    with open(spectral_input_path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
        )
        writer.writeheader()
        for shell in shells:
            writer.writerow(
                {
                    "atom": "Ar",
                    "orbital": shell["label"],
                    "type": shell["type"],
                    "n": shell["n"],
                    "l": shell["l"],
                    "occupation": f"{shell['occupation']:.12g}",
                    "energy_Ha": f"{shell['energy_Ha']:.12g}",
                }
            )

    print(f"Ar HF converged with E = {mf.e_tot:.12f} Ha using {basis}.")
    print(f"Wrote {solver_output_path} and {spectral_input_path}.")
    for shell in shells:
        print(
            f"Ar {shell['label']} {shell['type']} "
            f"l={shell['l']} occ={shell['occupation']:.6g} energy={shell['energy_Ha']:.8f}"
        )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export Ar spherically averaged PySCF HF radial orbitals.")
    parser.add_argument("--basis", default="cc-pvtz")
    parser.add_argument("--solver-output", default="atomic_solver_output.csv")
    parser.add_argument("--spectral-input", default="atomic_spectral_input.csv")
    parser.add_argument("--r-min", type=float, default=1e-5)
    parser.add_argument("--r-max", type=float, default=40.0)
    parser.add_argument("--n-grid", type=int, default=2400)
    args = parser.parse_args(argv)

    export_ar_radials(
        solver_output_path=args.solver_output,
        spectral_input_path=args.spectral_input,
        basis=args.basis,
        r_min=args.r_min,
        r_max=args.r_max,
        n_grid=args.n_grid,
    )


if __name__ == "__main__":
    main()
