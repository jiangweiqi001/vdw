import argparse
import csv
import re
from pathlib import Path

from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_all_e_rpa_atom import response_rows_from_arrays, write_channels, XC_MAP
from run_c6_table import build_c6_rows, write_c6_rows


ACTIVE_SHELLS_BY_ATOM_AND_NELEC = {
    ("Be", 2): "2s",
    ("Be", 4): "1s;2s",
    ("Ne", 8): "2s;2p",
    ("Ar", 8): "3s;3p",
    ("Kr", 8): "4s;4p",
    ("Mg", 2): "3s",
    ("Mg", 10): "2s;2p;3s",
    ("Ca", 2): "4s",
    ("Ca", 10): "3s;3p;4s",
}


def infer_active_shells(atom, nelec):
    return ACTIVE_SHELLS_BY_ATOM_AND_NELEC.get((atom, int(nelec)), "")


def summarize_psp_response(atom, psp, basis, xc, nstates, method, nelec, alpha_row):
    return {
        "atom": atom,
        "psp": psp,
        "basis": basis,
        "xc": xc,
        "nstates": int(nstates),
        "method": method,
        "active_electrons": int(nelec),
        "active_shells": infer_active_shells(atom, nelec),
        "alpha0_psp": float(alpha_row["alpha0_au"]),
        "C6_psp": float(alpha_row["C6_self_au"]),
        "n_channels": int(alpha_row["n_channels"]),
        "status": "ok",
        "note": "",
    }


def unavailable_row(atom, psp, basis, xc, nstates, reason, method="TDDFT"):
    return {
        "atom": atom,
        "psp": psp,
        "basis": basis,
        "xc": xc,
        "nstates": int(nstates),
        "method": method,
        "active_electrons": "",
        "active_shells": "",
        "alpha0_psp": "",
        "C6_psp": "",
        "n_channels": "",
        "status": "unavailable",
        "note": str(reason),
    }


def extract_cp2k_named_block(path, atom, name):
    text = Path(path).read_text(encoding="utf-8")
    lines = text.splitlines()
    start = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if parts and parts[0] == atom and name in parts[1:]:
            start = idx
            break
    if start is None:
        raise ValueError(f"No CP2K block for {atom} {name} in {path}.")

    end = len(lines)
    header_pattern = re.compile(r"^[A-Z][a-z]?\s+\S+")
    for idx in range(start + 1, len(lines)):
        stripped = lines[idx].strip()
        if stripped and not stripped.startswith("#") and header_pattern.match(stripped):
            end = idx
            break
    return "\n".join(lines[start:end])


def load_cp2k_basis(atom, basis_file, basis_name):
    from pyscf.gto.basis import parse_cp2k

    return parse_cp2k.parse(extract_cp2k_named_block(basis_file, atom, basis_name))


def load_cp2k_pseudo(atom, pseudo_file, pseudo_name):
    from pyscf.gto.basis import parse_cp2k_pp

    return parse_cp2k_pp.parse(extract_cp2k_named_block(pseudo_file, atom, pseudo_name))


def export_psp_response(atom, psp="gth-lda", basis="gth-dzvp", xc="lda", nstates=50, method="TDDFT", min_osc=1e-10):
    from pyscf import gto, scf, tdscf

    mol = gto.M(atom=f"{atom} 0 0 0", basis=basis, pseudo=psp, spin=0, charge=0, cart=False, verbose=0)
    mf = scf.RKS(mol)
    mf.xc = XC_MAP[xc]
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} PSP RKS-{xc} did not converge for {basis}/{psp}.")

    td_cls = tdscf.TDDFT if method.upper() == "TDDFT" else tdscf.TDA
    td = td_cls(mf)
    td.nstates = nstates
    td.verbose = 0
    energies = td.kernel()[0]
    oscillator_strengths = td.oscillator_strength()
    source = f"PySCF_PSP_{psp}_{method.upper()}_{xc.upper()}"
    rows = response_rows_from_arrays(atom, energies, oscillator_strengths, source=source, min_osc=min_osc)
    return rows, {"active_electrons": mol.nelectron}


def run_atom(
    atom,
    psp,
    basis,
    xc,
    nstates,
    method,
    output_root="results/psp_rpa",
    min_osc=1e-10,
    basis_file=None,
    basis_name=None,
    pseudo_file=None,
    pseudo_name=None,
):
    output_dir = Path(output_root) / atom.lower() / f"{psp}_{basis}_{xc}_{method.lower()}".replace("/", "_")
    output_dir.mkdir(parents=True, exist_ok=True)
    channels_path = output_dir / f"{atom.lower()}_psp_channels.csv"
    alpha_path = output_dir / "alpha_c6_table.csv"
    c6_path = output_dir / "c6_table.csv"

    if basis_file and basis_name:
        basis_obj = {atom: load_cp2k_basis(atom, basis_file, basis_name)}
        basis_label = basis_name
    else:
        basis_obj = basis
        basis_label = basis

    if pseudo_file and pseudo_name:
        pseudo_obj = {atom: load_cp2k_pseudo(atom, pseudo_file, pseudo_name)}
        psp_label = pseudo_name
    else:
        pseudo_obj = psp
        psp_label = psp

    rows, meta = export_psp_response(atom, pseudo_obj, basis_obj, xc, nstates, method, min_osc=min_osc)
    write_channels(channels_path, rows)
    alpha_rows = build_alpha_rows(channels_path)
    write_alpha_rows(alpha_path, alpha_rows)
    write_c6_rows(c6_path, build_c6_rows(channels_path))
    alpha_row = next(row for row in alpha_rows if row["atom"] == atom)
    return summarize_psp_response(atom, psp_label, basis_label, xc, nstates, method, meta["active_electrons"], alpha_row)


def write_summary(path, rows):
    fieldnames = [
        "atom",
        "psp",
        "basis",
        "xc",
        "nstates",
        "method",
        "active_electrons",
        "active_shells",
        "alpha0_psp",
        "C6_psp",
        "n_channels",
        "status",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run pseudopotential valence-only RKS TDDFT response for an atom.")
    parser.add_argument("--atom", required=True)
    parser.add_argument("--psp", default="gth-lda")
    parser.add_argument("--basis", default="gth-dzvp")
    parser.add_argument("--basis-file", default=None)
    parser.add_argument("--basis-name", default=None)
    parser.add_argument("--xc", choices=sorted(XC_MAP), default="lda")
    parser.add_argument("--pseudo-file", default=None)
    parser.add_argument("--pseudo-name", default=None)
    parser.add_argument("--nstates", type=int, default=50)
    parser.add_argument("--method", choices=["TDDFT", "TDA"], default="TDDFT")
    parser.add_argument("--output-root", default="results/psp_rpa")
    parser.add_argument("--summary", default="results/psp_rpa_summary.csv")
    args = parser.parse_args(argv)

    try:
        row = run_atom(
            args.atom,
            args.psp,
            args.basis,
            args.xc,
            args.nstates,
            args.method,
            args.output_root,
            basis_file=args.basis_file,
            basis_name=args.basis_name,
            pseudo_file=args.pseudo_file,
            pseudo_name=args.pseudo_name,
        )
    except Exception as exc:
        row = unavailable_row(args.atom, args.psp, args.basis, args.xc, args.nstates, f"{type(exc).__name__}: {exc}", args.method)

    summary_path = Path(args.summary)
    existing = []
    if summary_path.exists():
        with open(summary_path, newline="", encoding="utf-8") as fp:
            existing = list(csv.DictReader(fp))
        existing = [
            old for old in existing
            if not (
                old["atom"] == row["atom"]
                and old["psp"] == row["psp"]
                and old["basis"] == row["basis"]
                and old["xc"] == row["xc"]
                and int(old["nstates"]) == row["nstates"]
                and old["method"] == row["method"]
            )
        ]
    existing.append(row)
    write_summary(summary_path, existing)

    print("atom,psp,basis,xc,nstates,method,active_electrons,alpha0_psp,C6_psp,n_channels,status,note")
    print(
        f"{row['atom']},"
        f"{row['psp']},"
        f"{row['basis']},"
        f"{row['xc']},"
        f"{row['nstates']},"
        f"{row['method']},"
        f"{row['active_electrons']},"
        f"{row['alpha0_psp']},"
        f"{row['C6_psp']},"
        f"{row['n_channels']},"
        f"{row['status']},"
        f"{row['note']}"
    )


if __name__ == "__main__":
    main()
