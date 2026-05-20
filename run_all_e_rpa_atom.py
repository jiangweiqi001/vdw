import argparse
import csv
from pathlib import Path

from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from run_noble_gas_tdhf import load_references, percent_error


XC_MAP = {
    "lda": "lda,vwn",
    "pbe": "pbe,pbe",
}


def response_rows_from_arrays(atom, energies, oscillator_strengths, source, min_osc=1e-10):
    rows = []
    for idx, (energy, osc) in enumerate(zip(energies, oscillator_strengths), start=1):
        energy = float(energy)
        osc = float(osc)
        if energy <= 0.0 or osc <= min_osc:
            continue
        rows.append(
            {
                "atom": atom,
                "channel": f"rpa_{idx:03d}",
                "delta_Ha": energy,
                "osc": osc,
                "is_core": "true",
                "source": source,
            }
        )
    return rows


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


def export_all_e_response(atom, xc="lda", basis="aug-cc-pVTZ", nstates=100, method="TDDFT", min_osc=1e-10):
    from pyscf import gto, scf, tdscf

    mol = gto.M(atom=f"{atom} 0 0 0", basis=basis, spin=0, charge=0, cart=False, verbose=0)
    mf = scf.RKS(mol)
    mf.xc = XC_MAP[xc]
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} RKS-{xc} did not converge for {basis}.")

    td_cls = tdscf.TDDFT if method.upper() == "TDDFT" else tdscf.TDA
    td = td_cls(mf)
    td.nstates = nstates
    td.verbose = 0
    energies = td.kernel()[0]
    oscillator_strengths = td.oscillator_strength()
    source = f"PySCF_{method.upper()}_{xc.upper()}"
    rows = response_rows_from_arrays(atom, energies, oscillator_strengths, source=source, min_osc=min_osc)
    return rows, {
        "atom": atom,
        "xc": xc,
        "basis": basis,
        "nstates": int(nstates),
        "method": method.upper(),
        "n_channels": len(rows),
        "sum_osc": float(sum(row["osc"] for row in rows)),
    }


def summarize_response(atom, xc, basis, nstates, method, alpha_row, alpha0_ref, c6_ref):
    alpha0 = float(alpha_row["alpha0_au"])
    c6 = float(alpha_row["C6_self_au"])
    return {
        "atom": atom,
        "xc": xc,
        "basis": basis,
        "nstates": int(nstates),
        "method": method,
        "alpha0": alpha0,
        "C6": c6,
        "alpha0_ref": float(alpha0_ref),
        "alpha0_err": percent_error(alpha0, alpha0_ref),
        "C6_ref": float(c6_ref),
        "C6_err": percent_error(c6, c6_ref),
        "n_channels": int(alpha_row["n_channels"]),
    }


def run_atom(atom, xc, basis, nstates, method, output_root="results/all_e_rpa", min_osc=1e-10):
    output_dir = Path(output_root) / atom.lower() / f"{xc}_{method.lower()}_{basis.lower()}"
    output_dir.mkdir(parents=True, exist_ok=True)
    channels_path = output_dir / f"{atom.lower()}_all_e_channels.csv"
    alpha_path = output_dir / "alpha_c6_table.csv"
    c6_path = output_dir / "c6_table.csv"

    rows, meta = export_all_e_response(atom, xc, basis, nstates, method, min_osc=min_osc)
    write_channels(channels_path, rows)
    alpha_rows = build_alpha_rows(channels_path)
    write_alpha_rows(alpha_path, alpha_rows)
    write_c6_rows(c6_path, build_c6_rows(channels_path))
    alpha_row = next(row for row in alpha_rows if row["atom"] == atom)
    refs = load_references()
    return summarize_response(atom, xc, basis, nstates, meta["method"], alpha_row, refs[atom]["alpha0_ref"], refs[atom]["C6_ref"])


def write_summary(path, rows):
    fieldnames = [
        "atom",
        "xc",
        "basis",
        "nstates",
        "method",
        "alpha0",
        "C6",
        "alpha0_ref",
        "alpha0_err",
        "C6_ref",
        "C6_err",
        "n_channels",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run all-electron KS TDDFT/TDA response for an atom.")
    parser.add_argument("--atom", required=True)
    parser.add_argument("--xc", choices=sorted(XC_MAP), default="lda")
    parser.add_argument("--basis", default="aug-cc-pVTZ")
    parser.add_argument("--nstates", type=int, default=100)
    parser.add_argument("--method", choices=["TDDFT", "TDA"], default="TDDFT")
    parser.add_argument("--output-root", default="results/all_e_rpa")
    parser.add_argument("--summary", default="results/all_e_rpa_summary.csv")
    args = parser.parse_args(argv)

    row = run_atom(args.atom, args.xc, args.basis, args.nstates, args.method, args.output_root)
    summary_path = Path(args.summary)
    existing = []
    if summary_path.exists():
        with open(summary_path, newline="", encoding="utf-8") as fp:
            existing = list(csv.DictReader(fp))
        existing = [
            old for old in existing
            if not (
                old["atom"] == row["atom"]
                and old["xc"] == row["xc"]
                and old["basis"] == row["basis"]
                and int(old["nstates"]) == row["nstates"]
                and old["method"] == row["method"]
            )
        ]
    existing.append(row)
    write_summary(summary_path, existing)

    print("atom,xc,basis,nstates,method,alpha0,C6,alpha0_ref,alpha0_err,C6_ref,C6_err,n_channels")
    print(
        f"{row['atom']},"
        f"{row['xc']},"
        f"{row['basis']},"
        f"{row['nstates']},"
        f"{row['method']},"
        f"{row['alpha0']:.8f},"
        f"{row['C6']:.8f},"
        f"{row['alpha0_ref']:.8f},"
        f"{row['alpha0_err']:.6f},"
        f"{row['C6_ref']:.8f},"
        f"{row['C6_err']:.6f},"
        f"{row['n_channels']}"
    )


if __name__ == "__main__":
    main()
