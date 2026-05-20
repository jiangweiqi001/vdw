import argparse
import csv
from pathlib import Path

import numpy as np


BASIS_LIST = ["cc-pVTZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]


def mo_oscillator_channels_from_arrays(mo_energy, mo_occ, dipole_mo, atom="Ar"):
    mo_energy = np.asarray(mo_energy, dtype=float)
    mo_occ = np.asarray(mo_occ, dtype=float)
    dipole_mo = np.asarray(dipole_mo, dtype=float)
    occupied = np.where(mo_occ > 1e-12)[0]
    virtual = np.where(mo_occ <= 1e-12)[0]
    rows = []
    for i in occupied:
        for a in virtual:
            delta = mo_energy[a] - mo_energy[i]
            if delta <= 0.0:
                continue
            d2 = float(np.sum(np.abs(dipole_mo[:, i, a]) ** 2))
            osc = (2.0 / 3.0) * delta * mo_occ[i] * d2
            rows.append(
                {
                    "atom": atom,
                    "channel": f"mo_occ_{i:03d}_to_virt_{a:03d}",
                    "delta_Ha": float(delta),
                    "osc": float(osc),
                    "is_core": "true",
                    "source": "PySCF_MO",
                    "occ_mo": int(i),
                    "virt_mo": int(a),
                    "mo_occ": float(mo_occ[i]),
                    "d2_mo": d2,
                }
            )
    return rows


def export_ar_mo_oscillators(basis="cc-pVTZ"):
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
    rows = mo_oscillator_channels_from_arrays(mf.mo_energy, mf.mo_occ, dipole_mo, atom="Ar")
    n_electrons = float(np.sum(mf.mo_occ[mf.mo_occ > 1e-12]))
    return rows, {
        "basis": basis,
        "N_electrons": n_electrons,
        "sum_f_mo": sum(row["osc"] for row in rows),
        "sum_f_mo_over_N_electrons": sum(row["osc"] for row in rows) / n_electrons,
        "n_occ": int(np.count_nonzero(mf.mo_occ > 1e-12)),
        "n_virt": int(np.count_nonzero(mf.mo_occ <= 1e-12)),
    }


def write_mo_channels(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["atom", "channel", "delta_Ha", "osc", "is_core", "source"],
        )
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
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export Ar PySCF 3D MO oscillator-strength channels.")
    parser.add_argument("--basis", choices=BASIS_LIST, default="cc-pVTZ")
    parser.add_argument("--output", default="ar_mo_channels.csv")
    args = parser.parse_args(argv)

    rows, summary = export_ar_mo_oscillators(args.basis)
    write_mo_channels(args.output, rows)
    print(
        "basis,N_electrons,sum_f_mo,sum_f_mo/N_electrons,n_occ,n_virt,n_channels"
    )
    print(
        f"{summary['basis']},"
        f"{summary['N_electrons']:.8f},"
        f"{summary['sum_f_mo']:.8f},"
        f"{summary['sum_f_mo_over_N_electrons']:.8f},"
        f"{summary['n_occ']},"
        f"{summary['n_virt']},"
        f"{len(rows)}"
    )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
