import argparse
import csv
from pathlib import Path

import numpy as np


DEFAULT_TARGETS = {
    "Mg": {"basis": "cc-pVQZ", "shells": {"2s", "2p"}},
    "Ca": {"basis": "cc-pVQZ", "shells": {"3s", "3p"}},
}


def dipole_wilson_channels_from_arrays(atom, mo_energy, mo_occ, dipole_mo, mo_to_shell, selected_shells):
    mo_energy = np.asarray(mo_energy, dtype=float)
    mo_occ = np.asarray(mo_occ, dtype=float)
    dipole_mo = np.asarray(dipole_mo, dtype=float)
    occupied = np.where(mo_occ > 1e-12)[0]
    virtual = np.where(mo_occ <= 1e-12)[0]
    rows = []
    for i in occupied:
        shell = mo_to_shell.get(int(i))
        if shell not in selected_shells:
            continue
        for a in virtual:
            delta = float(mo_energy[a] - mo_energy[i])
            if delta <= 0.0:
                continue
            d2 = float(np.sum(np.abs(dipole_mo[:, i, a]) ** 2))
            osc = (2.0 / 3.0) * delta * float(mo_occ[i]) * d2
            if abs(osc) <= 1e-14:
                continue
            rows.append(
                {
                    "atom": atom,
                    "channel": f"eft_dipole_{shell}_occ_{int(i):03d}_to_virt_{int(a):03d}",
                    "delta_Ha": delta,
                    "osc": osc,
                    "is_core": "true",
                    "source": "EFT_CORE_DIPOLE_WILSON_MO_APPROX",
                    "from_shell": shell,
                    "occ_mo": int(i),
                    "virt_mo": int(a),
                    "occ_energy_Ha": float(mo_energy[i]),
                    "virt_energy_Ha": float(mo_energy[a]),
                    "d2_mo": d2,
                }
            )
    return rows


def _mo_shell_map(mol, mf):
    from pyscf_export_ar_radials import grouped_shells

    mapping = {}
    for shell in grouped_shells(mol, mf.mo_energy, mf.mo_occ, mf.mo_coeff):
        for column in shell["columns"]:
            mapping[int(column)] = shell["label"]
    return mapping


def compute_dipole_wilson_channels(atom, basis, selected_shells):
    from pyscf.scf import atom_hf
    from pyscf_export_ar_radials import make_atom_molecule

    mol = make_atom_molecule(atom, basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} atom HF did not converge for {basis}.")

    dipole_ao = mol.intor("int1e_r", comp=3)
    mo_coeff = mf.mo_coeff
    dipole_mo = np.asarray([mo_coeff.T @ dipole_ao[mu] @ mo_coeff for mu in range(3)])
    mo_to_shell = _mo_shell_map(mol, mf)
    return dipole_wilson_channels_from_arrays(
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
        "virt_mo",
        "occ_energy_Ha",
        "virt_energy_Ha",
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
                    "virt_mo": row["virt_mo"],
                    "occ_energy_Ha": f"{row['occ_energy_Ha']:.12f}",
                    "virt_energy_Ha": f"{row['virt_energy_Ha']:.12f}",
                    "d2_mo": f"{row['d2_mo']:.12f}",
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute l=1 dipole Wilson channels from selected core/semicore shells.")
    parser.add_argument("--atom", action="append")
    parser.add_argument("--basis", action="append")
    parser.add_argument("--shell", action="append")
    parser.add_argument("--output", default="results/eft_core_dipole_wilson_channels.csv")
    args = parser.parse_args(argv)

    all_rows = []
    if args.atom:
        basis_values = args.basis or ["cc-pVQZ"] * len(args.atom)
        if len(basis_values) == 1 and len(args.atom) > 1:
            basis_values = basis_values * len(args.atom)
        shells = set(args.shell or [])
        for atom, basis in zip(args.atom, basis_values):
            target_shells = shells or DEFAULT_TARGETS[atom]["shells"]
            all_rows.extend(compute_dipole_wilson_channels(atom, basis, target_shells))
    else:
        for atom, spec in DEFAULT_TARGETS.items():
            all_rows.extend(compute_dipole_wilson_channels(atom, spec["basis"], spec["shells"]))

    write_channels(args.output, all_rows)
    print("atom,from_shell,n_channels,sum_osc,output")
    grouped = {}
    for row in all_rows:
        grouped.setdefault((row["atom"], row["from_shell"]), 0.0)
        grouped[(row["atom"], row["from_shell"])] += row["osc"]
    for (atom, shell), sum_osc in sorted(grouped.items()):
        count = sum(1 for row in all_rows if row["atom"] == atom and row["from_shell"] == shell)
        print(f"{atom},{shell},{count},{sum_osc:.12f},{args.output}")


if __name__ == "__main__":
    main()
