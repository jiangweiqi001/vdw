import argparse
import csv
import re
from pathlib import Path

from run_psp_rpa_atom import load_cp2k_basis, load_cp2k_pseudo
from run_all_e_rpa_atom import XC_MAP


DEFAULT_ATOMS = ["Ca", "Sr", "Ba", "Zn", "Cd", "Hg"]
DEFAULT_PSEUDO_FILE = "external_data/cp2k/GTH_POTENTIALS"
DEFAULT_BASIS_FILES = [
    "external_data/cp2k/BASIS_MOLOPT_UCL",
    "external_data/cp2k/BASIS_MOLOPT_UZH",
    "external_data/cp2k/BASIS_MOLOPT",
]

Q2_ACTIVE_SHELLS = {
    "Mg": "3s",
    "Ca": "4s",
    "Sr": "5s",
    "Ba": "6s",
    "Zn": "4s",
    "Cd": "5s",
    "Hg": "6s",
}


def has_q_token(name, q_token):
    return re.search(rf"(^|[-_]){re.escape(q_token)}($|[-_])", name, flags=re.IGNORECASE) is not None


def discover_headers(lines, atom, required_token="q2"):
    names = []
    header_pattern = re.compile(r"^[A-Z][a-z]?\s+\S+")
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or not header_pattern.match(stripped):
            continue
        parts = stripped.split()
        if parts[0] != atom:
            continue
        aliases = parts[1:]
        for alias in aliases:
            if has_q_token(alias, required_token) and alias not in names:
                names.append(alias)
    return names


def discover_file_headers(path, atom, required_token="q2"):
    return discover_headers(Path(path).read_text(encoding="utf-8").splitlines(), atom, required_token)


def candidate_status(
    atom,
    pseudo_name,
    basis_name,
    can_build,
    can_run_rks,
    can_run_tddft,
    note,
    basis_provenance="library_native",
    basis_label=None,
    basis_note="",
):
    if not basis_name:
        status = "no_matched_q2_basis"
    elif not pseudo_name:
        status = "no_q2_pseudo"
    elif can_run_tddft:
        status = "tddft_smoke_ok"
    elif can_run_rks:
        status = "rks_ok"
    elif can_build:
        status = "build_ok"
    else:
        status = "build_failed"
    return {
        "atom": atom,
        "pseudo_name": pseudo_name,
        "basis_name": basis_name,
        "basis_label": basis_label or basis_name,
        "basis_provenance": basis_provenance,
        "active_electrons": 2,
        "active_shells": Q2_ACTIVE_SHELLS.get(atom, ""),
        "can_build_pyscf_mol": str(bool(can_build)).lower(),
        "can_run_rks": str(bool(can_run_rks)).lower(),
        "can_run_tddft_smoke": str(bool(can_run_tddft)).lower(),
        "candidate_status": status,
        "note": "; ".join(part for part in [note, basis_note] if part),
    }


def preferred_pseudos(pseudos):
    order = ["GTH-PBE-q2", "GTH-LDA-q2", "GTH-PADE-q2", "GTH-BLYP-q2", "GTH-BP-q2"]
    ranked = []
    for name in pseudos:
        rank = order.index(name) if name in order else len(order)
        ranked.append((rank, name))
    return [name for _rank, name in sorted(ranked)]


def preferred_bases(bases):
    order_tokens = ["QZ", "TZV2P", "TZVP", "DZVP", "SZV"]
    ranked = []
    for name in bases:
        rank = len(order_tokens)
        upper = name.upper()
        for idx, token in enumerate(order_tokens):
            if token in upper:
                rank = idx
                break
        ranked.append((rank, name))
    return [name for _rank, name in sorted(ranked)]


def basis_rank(name):
    ordered = preferred_bases([name])
    if not ordered:
        return (999, name)
    order_tokens = ["QZ", "TZV2P", "TZVP", "DZVP", "SZV"]
    upper = name.upper()
    for idx, token in enumerate(order_tokens):
        if token in upper:
            return (idx, name)
    return (len(order_tokens), name)


def smoke_candidate(atom, pseudo_file, pseudo_name, basis_file, basis_name, xc="pbe", nstates=3, run_tddft=True):
    from pyscf import gto, scf, tdscf

    pseudo = {atom: load_cp2k_pseudo(atom, pseudo_file, pseudo_name)}
    basis = {atom: load_cp2k_basis(atom, basis_file, basis_name)}
    mol = gto.M(atom=f"{atom} 0 0 0", basis=basis, pseudo=pseudo, spin=0, charge=0, cart=False, verbose=0)
    can_build = True
    mf = scf.RKS(mol)
    mf.xc = XC_MAP[xc]
    mf.verbose = 0
    mf.kernel()
    can_run_rks = bool(mf.converged)
    can_run_tddft = False
    note = f"nelectron={mol.nelectron}; nao={mol.nao_nr()}"
    if can_run_rks and run_tddft:
        td = tdscf.TDDFT(mf)
        td.nstates = nstates
        td.verbose = 0
        energies = td.kernel()[0]
        can_run_tddft = len(energies) > 0
        note += f"; tddft_states={len(energies)}"
    return can_build, can_run_rks, can_run_tddft, note


def load_imported_basis_candidates(path, atom=None):
    rows = []
    if not path:
        return rows
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            if atom and row.get("atom") != atom:
                continue
            rows.append(
                {
                    "atom": row["atom"],
                    "basis_file": row["basis_file"],
                    "basis_name": row["basis_name"],
                    "basis_label": row.get("basis_label") or row["basis_name"],
                    "basis_provenance": row.get("basis_provenance") or "imported_external",
                    "basis_note": row.get("basis_note", ""),
                }
            )
    return rows


def scan_candidates(atoms=None, pseudo_file=DEFAULT_PSEUDO_FILE, basis_files=None, run_tddft=True, imported_basis_candidates=None):
    rows = []
    basis_files = basis_files or DEFAULT_BASIS_FILES
    for atom in atoms or DEFAULT_ATOMS:
        pseudos = preferred_pseudos(discover_file_headers(pseudo_file, atom, "q2"))
        basis_options = []
        for basis_file in basis_files:
            for basis_name in discover_file_headers(basis_file, atom, "q2"):
                basis_options.append(
                    {
                        "basis_file": basis_file,
                        "basis_name": basis_name,
                        "basis_label": basis_name,
                        "basis_provenance": "library_native",
                        "basis_note": "",
                    }
                )
        basis_options = [row for row in basis_options if has_q_token(row["basis_name"], "q2")]
        for row in imported_basis_candidates or []:
            if row.get("atom") == atom:
                basis_options.append(row)
        basis_options = sorted(basis_options, key=lambda row: basis_rank(row["basis_label"]))
        if not pseudos:
            rows.append(candidate_status(atom, "", "", False, False, False, "no q2 pseudo found"))
            continue
        if not basis_options:
            rows.append(candidate_status(atom, pseudos[0], "", False, False, False, "no matched q2 basis found"))
            continue
        for pseudo_name in pseudos:
            for basis_candidate in basis_options:
                try:
                    can_build, can_run_rks, can_run_tddft, note = smoke_candidate(
                        atom,
                        pseudo_file,
                        pseudo_name,
                        basis_candidate["basis_file"],
                        basis_candidate["basis_name"],
                        run_tddft=run_tddft,
                    )
                except Exception as exc:
                    rows.append(
                        candidate_status(
                            atom,
                            pseudo_name,
                            basis_candidate["basis_label"],
                            False,
                            False,
                            False,
                            f"{type(exc).__name__}: {exc}",
                            basis_provenance=basis_candidate["basis_provenance"],
                            basis_note=basis_candidate.get("basis_note", ""),
                        )
                    )
                    continue
                row = candidate_status(
                    atom,
                    pseudo_name,
                    basis_candidate["basis_label"],
                    can_build,
                    can_run_rks,
                    can_run_tddft,
                    note,
                    basis_provenance=basis_candidate["basis_provenance"],
                    basis_note=basis_candidate.get("basis_note", ""),
                )
                row["basis_file"] = basis_candidate["basis_file"]
                row["basis_block_name"] = basis_candidate["basis_name"]
                row["pseudo_file"] = pseudo_file
                rows.append(row)
    return rows


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "pseudo_name",
        "basis_name",
        "basis_label",
        "basis_block_name",
        "basis_provenance",
        "basis_file",
        "pseudo_file",
        "active_electrons",
        "active_shells",
        "can_build_pyscf_mol",
        "can_run_rks",
        "can_run_tddft_smoke",
        "candidate_status",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main(argv=None):
    parser = argparse.ArgumentParser(description="Probe clean large-core q2 pseudo/basis candidates for a second benchmark.")
    parser.add_argument("--atom", action="append")
    parser.add_argument("--pseudo-file", default=DEFAULT_PSEUDO_FILE)
    parser.add_argument("--basis-file", action="append")
    parser.add_argument("--candidate-basis-csv", action="append")
    parser.add_argument("--no-tddft", action="store_true")
    parser.add_argument("--output", default="results/large_core_q2_candidate_scan.csv")
    args = parser.parse_args(argv)
    imported = []
    for path in args.candidate_basis_csv or []:
        imported.extend(load_imported_basis_candidates(path))
    rows = scan_candidates(
        args.atom or DEFAULT_ATOMS,
        args.pseudo_file,
        args.basis_file or DEFAULT_BASIS_FILES,
        run_tddft=not args.no_tddft,
        imported_basis_candidates=imported,
    )
    write_rows(args.output, rows)
    print("atom,pseudo_name,basis_name,candidate_status,note")
    for row in rows:
        print(f"{row['atom']},{row['pseudo_name']},{row['basis_name']},{row['candidate_status']},{row['note']}")


if __name__ == "__main__":
    main()
