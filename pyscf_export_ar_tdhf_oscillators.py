import argparse
import csv
from pathlib import Path

import numpy as np


BASIS_LIST = ["cc-pVTZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]


def tdhf_channel_rows_from_arrays(energies, oscillator_strengths, atom="Ar", min_osc=1e-10):
    rows = []
    for idx, (energy, osc) in enumerate(zip(energies, oscillator_strengths), start=1):
        energy = float(energy)
        osc = float(osc)
        if energy <= 0.0 or osc <= min_osc:
            continue
        rows.append(
            {
                "atom": atom,
                "channel": f"tdhf_{idx:03d}",
                "delta_Ha": energy,
                "osc": osc,
                "is_core": "true",
                "source": "PySCF_TDHF",
            }
        )
    return rows


def export_ar_tdhf_oscillators(basis="cc-pVTZ", nstates=100, min_osc=1e-10):
    from pyscf import tdscf
    from pyscf.scf import atom_hf
    from pyscf_export_ar_radials import make_ar_molecule

    mol = make_ar_molecule(basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"Ar atom HF did not converge for {basis}.")

    td = tdscf.TDHF(mf)
    td.nstates = nstates
    td.verbose = 0
    energies = np.asarray(td.kernel()[0], dtype=float)
    oscillator_strengths = np.asarray(td.oscillator_strength(), dtype=float)
    rows = tdhf_channel_rows_from_arrays(energies, oscillator_strengths, min_osc=min_osc)
    return rows, {
        "basis": basis,
        "nstates_requested": int(nstates),
        "nstates_returned": int(len(energies)),
        "n_channels": int(len(rows)),
        "sum_osc": float(sum(row["osc"] for row in rows)),
    }


def write_tdhf_channels(path, rows):
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
    parser = argparse.ArgumentParser(description="Export Ar PySCF TDHF oscillator-strength channels.")
    parser.add_argument("--basis", choices=BASIS_LIST, default="cc-pVTZ")
    parser.add_argument("--nstates", type=int, default=100)
    parser.add_argument("--min-osc", type=float, default=1e-10)
    parser.add_argument("--output", default="ar_tdhf_channels.csv")
    args = parser.parse_args(argv)

    rows, summary = export_ar_tdhf_oscillators(args.basis, args.nstates, args.min_osc)
    write_tdhf_channels(args.output, rows)
    print("basis,nstates_requested,nstates_returned,n_channels,sum_osc")
    print(
        f"{summary['basis']},"
        f"{summary['nstates_requested']},"
        f"{summary['nstates_returned']},"
        f"{summary['n_channels']},"
        f"{summary['sum_osc']:.8f}"
    )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
