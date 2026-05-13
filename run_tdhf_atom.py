import argparse
import csv
from pathlib import Path


def tdhf_rows_from_arrays(atom, energies, oscillator_strengths, min_osc=1e-10):
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


def export_tdhf_atom(atom, basis, nstates=200, min_osc=1e-10):
    from pyscf import gto, tdscf
    from pyscf.scf import atom_hf

    mol = gto.M(atom=f"{atom} 0 0 0", basis=basis, spin=0, charge=0, cart=False, verbose=0)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} atom HF did not converge for {basis}.")

    td = tdscf.TDHF(mf)
    td.nstates = nstates
    td.verbose = 0
    energies = td.kernel()[0]
    oscillator_strengths = td.oscillator_strength()
    rows = tdhf_rows_from_arrays(atom, energies, oscillator_strengths, min_osc=min_osc)
    return rows, {
        "atom": atom,
        "basis": basis,
        "nstates_requested": int(nstates),
        "nstates_returned": int(len(energies)),
        "n_channels": int(len(rows)),
        "sum_osc": float(sum(row["osc"] for row in rows)),
    }


def write_channels(path, rows):
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
    parser = argparse.ArgumentParser(description="Export TDHF oscillator-strength channels for a closed-shell atom.")
    parser.add_argument("--atom", required=True)
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--nstates", type=int, default=200)
    parser.add_argument("--min-osc", type=float, default=1e-10)
    parser.add_argument("--output", default=None)
    args = parser.parse_args(argv)

    output = args.output or f"{args.atom.lower()}_tdhf_channels.csv"
    rows, summary = export_tdhf_atom(args.atom, args.basis, args.nstates, args.min_osc)
    write_channels(output, rows)
    print("atom,basis,nstates_requested,nstates_returned,n_channels,sum_osc")
    print(
        f"{summary['atom']},"
        f"{summary['basis']},"
        f"{summary['nstates_requested']},"
        f"{summary['nstates_returned']},"
        f"{summary['n_channels']},"
        f"{summary['sum_osc']:.8f}"
    )
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()
