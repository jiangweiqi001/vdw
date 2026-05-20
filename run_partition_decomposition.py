import argparse
import csv
from pathlib import Path

import numpy as np

from pyscf_export_ar_tdhf_decomposed import oscillator_from_dipole, write_component_channels
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from run_noble_gas_tdhf import load_references


def load_partition_definitions(path, atom):
    partitions = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            if row["atom"].strip() != atom:
                continue
            shells = {shell.strip() for shell in row["shells"].split(";") if shell.strip()}
            partitions[row["partition"].strip()] = shells
    if not partitions:
        raise ValueError(f"No partition definitions found for {atom}.")
    return partitions


def partition_combinations(partitions):
    combos = {name: set(shells) for name, shells in partitions.items()}
    if "valence" in partitions and "semicore" in partitions:
        combos["valence_plus_semicore"] = set(partitions["valence"]) | set(partitions["semicore"])
    if {"valence", "semicore", "deep_core"} <= set(partitions):
        combos["valence_plus_semicore_plus_deep_core"] = (
            set(partitions["valence"]) | set(partitions["semicore"]) | set(partitions["deep_core"])
        )
    return combos


def _mo_shell_map(mol, mf):
    from pyscf_export_ar_radials import grouped_shells

    mapping = {}
    for shell in grouped_shells(mol, mf.mo_energy, mf.mo_occ, mf.mo_coeff):
        for column in shell["columns"]:
            mapping[int(column)] = shell["label"]
    return mapping


def projected_dipoles_by_component(atom, basis, nstates, component_shells):
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
        [mo_coeff[:, occupied].T @ dipole_ao[mu] @ mo_coeff[:, virtual] for mu in range(3)]
    )
    mo_to_shell = _mo_shell_map(mol, mf)
    masks = {
        name: np.asarray([mo_to_shell[int(i)] in shells for i in occupied])
        for name, shells in component_shells.items()
    }

    dipoles = {name: [] for name in component_shells}
    dipoles["all"] = []
    for x, y in td.xy:
        amplitudes = x + y
        dipoles["all"].append(2.0 * np.einsum("xij,ij->x", dipole_ov, amplitudes))
        for name, mask in masks.items():
            dipoles[name].append(2.0 * np.einsum("xij,ij->x", dipole_ov[:, mask, :], amplitudes[mask, :]))
    return energies, {name: np.asarray(values) for name, values in dipoles.items()}


def rows_from_dipoles(atom, component, energies, dipoles, min_osc=1e-10):
    rows = []
    for idx, (energy, dipole) in enumerate(zip(energies, dipoles), start=1):
        osc = oscillator_from_dipole(energy, dipole)
        if abs(osc) <= min_osc:
            continue
        rows.append(
            {
                "atom": atom,
                "channel": f"tdhf_{idx:03d}",
                "delta_Ha": float(energy),
                "osc": float(osc),
                "is_core": "false" if component == "valence" else "true",
                "source": "PySCF_TDHF_PROJECTED",
                "component": component,
            }
        )
    return rows


def analyze_component(path, prefix, atom):
    alpha_rows = build_alpha_rows(path)
    write_alpha_rows(f"{prefix}_alpha_c6_table.csv", alpha_rows)
    write_c6_rows(f"{prefix}_c6_table.csv", build_c6_rows(path))
    return next(
        (row for row in alpha_rows if row["atom"] == atom),
        {"atom": atom, "alpha0_au": "0.00000000", "C6_self_au": "0.00000000", "n_channels": 0},
    )


def percent_error(c6, atom):
    refs = load_references()
    return 100.0 * (float(c6) - refs[atom]["C6_ref"]) / refs[atom]["C6_ref"]


def run_partition_decomposition(
    atom,
    basis="aug-cc-pVQZ",
    nstates=200,
    partition_file="partition_definitions.csv",
    output_root="results",
    min_osc=1e-10,
):
    partitions = load_partition_definitions(partition_file, atom)
    components = partition_combinations(partitions)
    energies, dipoles = projected_dipoles_by_component(atom, basis, nstates, components)

    output_dir = Path(output_root) / atom.lower() / "partition_decomposition"
    output_dir.mkdir(parents=True, exist_ok=True)
    all_components = {"all": None, **components}
    summaries = {}
    for component in all_components:
        rows = rows_from_dipoles(atom, component, energies, dipoles[component], min_osc=min_osc)
        path = output_dir / f"{atom.lower()}_{component}_channels.csv"
        write_component_channels(path, rows)
        summaries[component] = analyze_component(path, output_dir / component, atom)

    results = []
    for component, summary in summaries.items():
        c6 = float(summary["C6_self_au"])
        results.append(
            {
                "component": component,
                "alpha0": float(summary["alpha0_au"]),
                "C6": c6,
                "C6_error_pct": percent_error(c6, atom),
                "n_channels": int(summary["n_channels"]),
                "note": "TDHF transition-dipole projection",
            }
        )

    if "valence" in summaries and "valence_plus_semicore" in summaries:
        c6_val = float(summaries["valence"]["C6_self_au"])
        c6_val_sc = float(summaries["valence_plus_semicore"]["C6_self_au"])
        c6_all = float(summaries["all"]["C6_self_au"])
        delta = c6_val_sc - c6_val
        results.append(
            {
                "component": "Delta_C6_semicore",
                "alpha0": "",
                "C6": delta,
                "C6_error_pct": "",
                "n_channels": "",
                "note": f"relative_semicore_contribution={delta / c6_all:.8f}",
            }
        )

    output_path = Path(output_root) / atom.lower() / f"{atom.lower()}_core_valence_decomposition.csv"
    with open(output_path, "w", newline="", encoding="utf-8") as fp:
        fieldnames = ["component", "alpha0", "C6", "C6_error_pct", "n_channels", "note"]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run TDHF partition decomposition for an atom.")
    parser.add_argument("--atom", required=True)
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--nstates", type=int, default=200)
    parser.add_argument("--partition-file", default="partition_definitions.csv")
    parser.add_argument("--output-root", default="results")
    args = parser.parse_args(argv)

    rows = run_partition_decomposition(args.atom, args.basis, args.nstates, args.partition_file, args.output_root)
    print("component,alpha0,C6,C6_error_pct,n_channels,note")
    for row in rows:
        print(f"{row['component']},{row['alpha0']},{row['C6']},{row['C6_error_pct']},{row['n_channels']},{row['note']}")


if __name__ == "__main__":
    main()
