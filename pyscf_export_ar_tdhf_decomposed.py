import argparse
import csv
from pathlib import Path

import numpy as np

AR_VALENCE_SHELLS = {"3s", "3p"}
AR_CORE_SHELLS = {"1s", "2s", "2p"}


def oscillator_from_dipole(delta, dipole):
    dipole = np.asarray(dipole, dtype=float)
    return float((2.0 / 3.0) * float(delta) * np.dot(dipole, dipole))


def tdhf_component_rows(energies, all_dipoles, valence_dipoles, core_dipoles, atom="Ar", min_osc=1e-10):
    rows = {"all": [], "valence": [], "core": [], "cross": []}
    for idx, (energy, d_all, d_val, d_core) in enumerate(
        zip(energies, all_dipoles, valence_dipoles, core_dipoles),
        start=1,
    ):
        f_all = oscillator_from_dipole(energy, d_all)
        f_val = oscillator_from_dipole(energy, d_val)
        f_core = oscillator_from_dipole(energy, d_core)
        f_cross = f_all - f_val - f_core
        values = {
            "all": f_all,
            "valence": f_val,
            "core": f_core,
            "cross": f_cross,
        }
        for component, osc in values.items():
            if component != "cross" and abs(osc) <= min_osc:
                continue
            if component == "cross" and abs(osc) <= min_osc:
                continue
            rows[component].append(
                {
                    "atom": atom,
                    "channel": f"tdhf_{idx:03d}",
                    "delta_Ha": float(energy),
                    "osc": float(osc),
                    "is_core": "true" if component in {"all", "core", "cross"} else "false",
                    "source": "PySCF_TDHF_PROJECTED",
                    "component": component,
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


def export_decomposed_tdhf(
    atom="Ar",
    basis="aug-cc-pVQZ",
    nstates=200,
    valence_shells=AR_VALENCE_SHELLS,
    core_shells=AR_CORE_SHELLS,
    min_osc=1e-10,
):
    from pyscf import tdscf
    from pyscf.scf import atom_hf
    from pyscf_export_ar_radials import make_atom_molecule

    mol = make_atom_molecule(atom, basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} atom HF did not converge for {basis}.")

    td = tdscf.TDHF(mf)
    td.nstates = nstates
    td.verbose = 0
    energies = np.asarray(td.kernel()[0], dtype=float)

    dipole_ao = mol.intor_symmetric("int1e_r", comp=3)
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    occupied = np.where(mo_occ == 2)[0]
    virtual = np.where(mo_occ == 0)[0]
    dipole_ov = np.asarray(
        [
            mo_coeff[:, occupied].T @ dipole_ao[mu] @ mo_coeff[:, virtual]
            for mu in range(3)
        ]
    )
    mo_to_shell = _mo_shell_map(mol, mf)
    valence_mask = np.asarray([mo_to_shell[int(i)] in valence_shells for i in occupied])
    core_mask = np.asarray([mo_to_shell[int(i)] in core_shells for i in occupied])

    all_dipoles = []
    valence_dipoles = []
    core_dipoles = []
    for x, y in td.xy:
        amplitudes = x + y
        all_dipoles.append(2.0 * np.einsum("xij,ij->x", dipole_ov, amplitudes))
        valence_dipoles.append(2.0 * np.einsum("xij,ij->x", dipole_ov[:, valence_mask, :], amplitudes[valence_mask, :]))
        core_dipoles.append(2.0 * np.einsum("xij,ij->x", dipole_ov[:, core_mask, :], amplitudes[core_mask, :]))

    rows = tdhf_component_rows(
        energies,
        np.asarray(all_dipoles),
        np.asarray(valence_dipoles),
        np.asarray(core_dipoles),
        atom=atom,
        min_osc=min_osc,
    )
    return rows, {
        "atom": atom,
        "basis": basis,
        "nstates": int(nstates),
        "n_all": len(rows["all"]),
        "n_valence": len(rows["valence"]),
        "n_core": len(rows["core"]),
        "n_cross": len(rows["cross"]),
    }


def export_decomposed_ar_tdhf(basis="aug-cc-pVQZ", nstates=200, min_osc=1e-10):
    return export_decomposed_tdhf(
        atom="Ar",
        basis=basis,
        nstates=nstates,
        valence_shells=AR_VALENCE_SHELLS,
        core_shells=AR_CORE_SHELLS,
        min_osc=min_osc,
    )


def write_component_channels(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["atom", "channel", "delta_Ha", "osc", "is_core", "source", "component"],
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
                    "component": row["component"],
                }
            )


def export_to_dir(output_dir, basis="aug-cc-pVQZ", nstates=200, min_osc=1e-10):
    rows, summary = export_decomposed_ar_tdhf(basis, nstates, min_osc)
    output_dir = Path(output_dir)
    write_component_channels(output_dir / "ar_tdhf_all_channels.csv", rows["all"])
    write_component_channels(output_dir / "ar_tdhf_valence_projected_channels.csv", rows["valence"])
    write_component_channels(output_dir / "ar_tdhf_core_projected_channels.csv", rows["core"])
    write_component_channels(output_dir / "ar_tdhf_cross_projected_channels.csv", rows["cross"])
    return rows, summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export projected Ar TDHF valence/core/cross oscillator channels.")
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--nstates", type=int, default=200)
    parser.add_argument("--output-dir", default="results/ar/core_valence_tdhf")
    parser.add_argument("--min-osc", type=float, default=1e-10)
    args = parser.parse_args(argv)

    _rows, summary = export_to_dir(args.output_dir, args.basis, args.nstates, args.min_osc)
    print("basis,nstates,n_all,n_valence,n_core,n_cross")
    print(
        f"{summary['basis']},"
        f"{summary['nstates']},"
        f"{summary['n_all']},"
        f"{summary['n_valence']},"
        f"{summary['n_core']},"
        f"{summary['n_cross']}"
    )


if __name__ == "__main__":
    main()
