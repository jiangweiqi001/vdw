import argparse
import csv
from pathlib import Path

import numpy as np


CORE_ION_DEFAULTS = {
    "Mg": {"charge": 2, "basis": "aug-cc-pVQZ", "nstates": 100},
    "Ca": {"charge": 2, "basis": "cc-pVQZ", "nstates": 100},
}


def transition_density_dipole(dipole_ov, amplitudes):
    """Contract TDHF transition density amplitudes with <occ|r|virt>.

    PySCF RHF TDHF transition dipoles use 2 * sum_ia (X+Y)_ia <i|r|a>.
    The factor of 2 accounts for spin in a closed-shell singlet response.
    """
    return 2.0 * np.einsum("xij,ij->x", np.asarray(dipole_ov), np.asarray(amplitudes))


def oscillator_from_transition_dipole(atom, state_index, energy, dipole_vector, source):
    dipole_vector = np.asarray(dipole_vector, dtype=float)
    d2 = float(np.dot(dipole_vector, dipole_vector))
    energy = float(energy)
    return {
        "atom": atom,
        "channel": f"multipole_tdhf_{int(state_index):03d}",
        "delta_Ha": energy,
        "d_x": float(dipole_vector[0]),
        "d_y": float(dipole_vector[1]),
        "d_z": float(dipole_vector[2]),
        "d2": d2,
        "osc": float((2.0 / 3.0) * energy * d2),
        "is_core": "true",
        "source": source,
    }


def compute_multipole_core_channels(atom, charge, basis, nstates=100, min_osc=1e-10):
    from pyscf import gto, scf, tdscf

    mol = gto.M(atom=f"{atom} 0 0 0", basis=basis, charge=int(charge), spin=0, cart=False, verbose=0)
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom}^{charge}+ core ion RHF did not converge for {basis}.")

    td = tdscf.TDHF(mf)
    td.nstates = int(nstates)
    td.verbose = 0
    energies = np.asarray(td.kernel()[0], dtype=float)

    dipole_ao = mol.intor_symmetric("int1e_r", comp=3)
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    occ = np.where(mo_occ == 2)[0]
    virt = np.where(mo_occ == 0)[0]
    dipole_ov = np.asarray([mo_coeff[:, occ].T @ dipole_ao[mu] @ mo_coeff[:, virt] for mu in range(3)])

    rows = []
    dipoles = []
    for idx, (energy, (x, y)) in enumerate(zip(energies, td.xy), start=1):
        dipole = transition_density_dipole(dipole_ov, x + y)
        row = oscillator_from_transition_dipole(
            atom=atom,
            state_index=idx,
            energy=energy,
            dipole_vector=dipole,
            source="EFT_CORE_MULTIPOLE_TDENSITY_TDHF",
        )
        if row["osc"] > min_osc:
            rows.append(row)
        dipoles.append(dipole)

    return rows, {
        "atom": atom,
        "charge": int(charge),
        "basis": basis,
        "nstates": int(nstates),
        "n_channels": len(rows),
        "sum_osc": sum(row["osc"] for row in rows),
    }


def write_channels(path, rows):
    fieldnames = ["atom", "channel", "delta_Ha", "d_x", "d_y", "d_z", "d2", "osc", "is_core", "source"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "atom": row["atom"],
                    "channel": row["channel"],
                    "delta_Ha": f"{row['delta_Ha']:.12f}",
                    "d_x": f"{row['d_x']:.12f}",
                    "d_y": f"{row['d_y']:.12f}",
                    "d_z": f"{row['d_z']:.12f}",
                    "d2": f"{row['d2']:.12f}",
                    "osc": f"{row['osc']:.12f}",
                    "is_core": row["is_core"],
                    "source": row["source"],
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute explicit multipole transition-density core Wilson channels.")
    parser.add_argument("--atom", action="append")
    parser.add_argument("--charge", action="append", type=int)
    parser.add_argument("--basis", action="append")
    parser.add_argument("--nstates", type=int, default=100)
    parser.add_argument("--output", default="results/eft_core_multipole_wilson_channels.csv")
    args = parser.parse_args(argv)

    all_rows = []
    summaries = []
    atoms = args.atom or list(CORE_ION_DEFAULTS)
    charges = args.charge or [CORE_ION_DEFAULTS[atom]["charge"] for atom in atoms]
    bases = args.basis or [CORE_ION_DEFAULTS[atom]["basis"] for atom in atoms]
    if len(charges) == 1 and len(atoms) > 1:
        charges *= len(atoms)
    if len(bases) == 1 and len(atoms) > 1:
        bases *= len(atoms)
    for atom, charge, basis in zip(atoms, charges, bases):
        rows, summary = compute_multipole_core_channels(atom, charge, basis, args.nstates)
        all_rows.extend(rows)
        summaries.append(summary)
    write_channels(args.output, all_rows)
    print("atom,charge,basis,nstates,n_channels,sum_osc,output")
    for summary in summaries:
        print(
            f"{summary['atom']},{summary['charge']},{summary['basis']},"
            f"{summary['nstates']},{summary['n_channels']},{summary['sum_osc']:.12f},{args.output}"
        )


if __name__ == "__main__":
    main()
