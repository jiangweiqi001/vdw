import argparse
import csv

import numpy as np


BASIS_LIST = ["cc-pVTZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]


def compute_mo_trk_for_basis(basis):
    from pyscf.scf import atom_hf
    from pyscf_export_ar_radials import make_ar_molecule

    mol = make_ar_molecule(basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"Ar atom HF did not converge for {basis}.")

    dipole_ao = mol.intor("int1e_r", comp=3)
    mo_coeff = mf.mo_coeff
    dipole_mo = np.asarray([mo_coeff.T @ dipole_ao[mu] @ mo_coeff for mu in range(3)])
    mo_energy = mf.mo_energy
    mo_occ = mf.mo_occ

    occupied = np.where(mo_occ > 1e-12)[0]
    virtual = np.where(mo_occ <= 1e-12)[0]
    sum_f = 0.0
    n_transitions = 0
    for i in occupied:
        for a in virtual:
            delta = mo_energy[a] - mo_energy[i]
            if delta <= 0.0:
                continue
            d2 = float(np.sum(np.abs(dipole_mo[:, i, a]) ** 2))
            sum_f += (2.0 / 3.0) * delta * mo_occ[i] * d2
            n_transitions += 1

    n_electrons = float(np.sum(mo_occ[occupied]))
    return {
        "basis": basis,
        "N_electrons": n_electrons,
        "sum_f_mo": sum_f,
        "sum_f_mo_over_N_electrons": sum_f / n_electrons if n_electrons else 0.0,
        "n_occ_mo": int(len(occupied)),
        "n_virt_mo": int(len(virtual)),
        "n_transitions": int(n_transitions),
    }


def run_mo_trk(bases=None):
    return [compute_mo_trk_for_basis(basis) for basis in (bases or BASIS_LIST)]


def write_rows(path, rows):
    fieldnames = [
        "basis",
        "N_electrons",
        "sum_f_mo",
        "sum_f_mo_over_N_electrons",
        "n_occ_mo",
        "n_virt_mo",
        "n_transitions",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Check TRK sum rule from PySCF 3D AO dipole integrals and MO coefficients.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST)
    parser.add_argument("--output", default="results/ar/ar_mo_trk.csv")
    args = parser.parse_args(argv)

    rows = run_mo_trk(args.basis or BASIS_LIST)
    write_rows(args.output, rows)
    print("basis,N_electrons,sum_f_mo,sum_f_mo/N_electrons,n_occ_mo,n_virt_mo,n_transitions")
    for row in rows:
        print(
            f"{row['basis']},"
            f"{row['N_electrons']:.8f},"
            f"{row['sum_f_mo']:.8f},"
            f"{row['sum_f_mo_over_N_electrons']:.8f},"
            f"{row['n_occ_mo']},"
            f"{row['n_virt_mo']},"
            f"{row['n_transitions']}"
        )


if __name__ == "__main__":
    main()
