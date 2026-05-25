import argparse
import csv
from pathlib import Path


CORE_ION_DEFAULTS = {
    "Mg": {"charge": 2, "basis": "aug-cc-pVQZ", "nstates": 100},
    "Ca": {"charge": 2, "basis": "cc-pVQZ", "nstates": 100},
}


def core_tdhf_rows_from_arrays(atom, energies, oscillator_strengths, min_osc=1e-10):
    rows = []
    for idx, (energy, osc) in enumerate(zip(energies, oscillator_strengths), start=1):
        energy = float(energy)
        osc = float(osc)
        if energy <= 0.0 or osc <= min_osc:
            continue
        rows.append(
            {
                "atom": atom,
                "channel": f"core_tdhf_{idx:03d}",
                "delta_Ha": energy,
                "osc": osc,
                "is_core": "true",
                "source": "EFT_CORE_DIPOLE_WILSON_CORE_TDHF",
            }
        )
    return rows


def compute_core_tdhf_channels(atom, charge, basis, nstates=100, min_osc=1e-10):
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
    energies = td.kernel()[0]
    oscillator_strengths = td.oscillator_strength()
    rows = core_tdhf_rows_from_arrays(atom, energies, oscillator_strengths, min_osc=min_osc)
    return rows, {
        "atom": atom,
        "charge": int(charge),
        "basis": basis,
        "nstates": int(nstates),
        "n_channels": len(rows),
        "sum_osc": sum(row["osc"] for row in rows),
    }


def write_channels(path, rows):
    fieldnames = ["atom", "channel", "delta_Ha", "osc", "is_core", "source"]
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
                    "osc": f"{row['osc']:.12f}",
                    "is_core": row["is_core"],
                    "source": row["source"],
                }
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute core-ion TDHF dipole Wilson channels.")
    parser.add_argument("--atom", action="append")
    parser.add_argument("--charge", action="append", type=int)
    parser.add_argument("--basis", action="append")
    parser.add_argument("--nstates", type=int, default=100)
    parser.add_argument("--output", default="results/eft_core_tdhf_wilson_channels.csv")
    args = parser.parse_args(argv)

    all_rows = []
    summaries = []
    if args.atom:
        charges = args.charge or [CORE_ION_DEFAULTS[atom]["charge"] for atom in args.atom]
        bases = args.basis or [CORE_ION_DEFAULTS[atom]["basis"] for atom in args.atom]
        if len(charges) == 1 and len(args.atom) > 1:
            charges = charges * len(args.atom)
        if len(bases) == 1 and len(args.atom) > 1:
            bases = bases * len(args.atom)
        for atom, charge, basis in zip(args.atom, charges, bases):
            rows, summary = compute_core_tdhf_channels(atom, charge, basis, args.nstates)
            all_rows.extend(rows)
            summaries.append(summary)
    else:
        for atom, spec in CORE_ION_DEFAULTS.items():
            rows, summary = compute_core_tdhf_channels(atom, spec["charge"], spec["basis"], args.nstates)
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
