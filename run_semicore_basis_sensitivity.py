import argparse
import csv
from pathlib import Path

from compute_core_sternheimer import compute_core_sternheimer_channels, write_channels as write_core_channels
from run_psp_rpa_atom import run_atom as run_psp_atom
from run_semicore_c6_validation import run_validation
from semicore_c6_targets import parse_shells


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def infer_xc_from_pseudo(pseudo_name):
    upper = str(pseudo_name).upper()
    if "PBE" in upper or "GGA" in upper:
        return "pbe"
    if "LDA" in upper or "PADE" in upper:
        return "lda"
    return "pbe"


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "pseudo_name",
        "basis_label",
        "basis_block_name",
        "basis_provenance",
        "active_shells",
        "correction_shells",
        "C6_PSP",
        "C6_PSP_plus_sternheimer",
        "Delta_C6_core",
        "relative_delta_pct",
        "target_audit_status",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _candidate_rank(row):
    basis = row.get("basis_label") or row.get("basis_name", "")
    pseudo = row.get("pseudo_name", "")
    pseudo_rank = 0 if pseudo == "GTH-LDA-q2" else 1
    basis_tokens = ["QZVPP", "TZV2P", "TZVPP", "TZVP", "DZVP", "SVP"]
    basis_rank = next((idx for idx, token in enumerate(basis_tokens) if token in basis.upper()), 99)
    return (pseudo_rank, basis_rank, basis)


def select_candidates(candidate_rows, atoms=None, pseudo_name="GTH-LDA-q2"):
    atoms = set(atoms or [])
    rows = []
    for row in candidate_rows:
        if atoms and row.get("atom") not in atoms:
            continue
        if row.get("candidate_status") != "tddft_smoke_ok":
            continue
        if pseudo_name and row.get("pseudo_name") != pseudo_name:
            continue
        rows.append(row)
    return sorted(rows, key=_candidate_rank)


def _channels_path(output_root, atom, pseudo_name, basis_label, xc, method):
    leaf = f"{pseudo_name}_{basis_label}_{xc}_{method.lower()}".replace("/", "_")
    return Path(output_root) / "psp_rpa" / atom.lower() / leaf / f"{atom.lower()}_psp_channels.csv"


def run_sensitivity(
    candidate_scan,
    atoms,
    output_root,
    pseudo_name="GTH-LDA-q2",
    psp_nstates=100,
    method="TDDFT",
    core_basis="def2-TZVPP",
):
    from semicore_c6_targets import SEMICORE_C6_TARGETS

    output_root = Path(output_root)
    candidates = select_candidates(read_rows(candidate_scan), atoms=atoms, pseudo_name=pseudo_name)
    core_paths = {}
    summaries = []

    for atom in sorted({row["atom"] for row in candidates}):
        target = SEMICORE_C6_TARGETS[atom]
        core_path = output_root / "core_sternheimer" / f"{atom.lower()}_core_sternheimer_channels.csv"
        core_rows = compute_core_sternheimer_channels(atom, core_basis, target.correction_shells)
        write_core_channels(core_path, core_rows)
        core_paths[atom] = core_path

    for candidate in candidates:
        atom = candidate["atom"]
        target = SEMICORE_C6_TARGETS[atom]
        basis_label = candidate.get("basis_label") or candidate["basis_name"]
        basis_name = candidate.get("basis_block_name") or candidate.get("basis_name") or basis_label
        xc = infer_xc_from_pseudo(candidate["pseudo_name"])
        run_psp_atom(
            atom=atom,
            psp=candidate["pseudo_name"],
            basis=basis_label,
            xc=xc,
            nstates=psp_nstates,
            method=method,
            output_root=output_root / "psp_rpa",
            basis_file=candidate.get("basis_file", ""),
            basis_name=basis_name,
            pseudo_file=candidate.get("pseudo_file", ""),
            pseudo_name=candidate["pseudo_name"],
        )
        summary = run_validation(
            atom=atom,
            psp_channels=_channels_path(output_root, atom, candidate["pseudo_name"], basis_label, xc, method),
            core_channels=core_paths[atom],
            output_root=output_root / "validation" / atom.lower() / basis_label,
            active_electrons=target.active_electrons,
            active_shells=";".join(sorted(target.explicit_shells)),
        )
        c6_psp = float(summary["C6_PSP"])
        delta = float(summary["Delta_C6_core"])
        summaries.append(
            {
                "atom": atom,
                "pseudo_name": candidate["pseudo_name"],
                "basis_label": basis_label,
                "basis_block_name": basis_name,
                "basis_provenance": candidate.get("basis_provenance", ""),
                "active_shells": ";".join(sorted(target.explicit_shells)),
                "correction_shells": ";".join(sorted(parse_shells(summary["correction_shells"]))),
                "C6_PSP": summary["C6_PSP"],
                "C6_PSP_plus_sternheimer": summary["C6_PSP_plus_sternheimer"],
                "Delta_C6_core": summary["Delta_C6_core"],
                "relative_delta_pct": 100.0 * delta / c6_psp if c6_psp else "",
                "target_audit_status": summary["target_audit_status"],
                "note": candidate.get("note", ""),
            }
        )
    return summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Zn/Cd adapted-basis semicore C6 sensitivity.")
    parser.add_argument("--candidate-scan", required=True)
    parser.add_argument("--atom", action="append", choices=["Zn", "Cd"])
    parser.add_argument("--pseudo-name", default="GTH-LDA-q2")
    parser.add_argument("--output-root", default="results/semicore_zn_cd_basis_sensitivity")
    parser.add_argument("--psp-nstates", type=int, default=100)
    parser.add_argument("--method", choices=["TDDFT", "TDA"], default="TDDFT")
    parser.add_argument("--core-basis", default="def2-TZVPP")
    args = parser.parse_args(argv)

    summaries = run_sensitivity(
        candidate_scan=args.candidate_scan,
        atoms=args.atom or ["Zn", "Cd"],
        output_root=args.output_root,
        pseudo_name=args.pseudo_name,
        psp_nstates=args.psp_nstates,
        method=args.method,
        core_basis=args.core_basis,
    )
    output = Path(args.output_root) / "summary.csv"
    write_rows(output, summaries)
    print("atom,pseudo_name,basis_label,C6_PSP,C6_PSP_plus_sternheimer,Delta_C6_core,relative_delta_pct")
    for row in summaries:
        print(
            f"{row['atom']},"
            f"{row['pseudo_name']},"
            f"{row['basis_label']},"
            f"{row['C6_PSP']},"
            f"{row['C6_PSP_plus_sternheimer']},"
            f"{row['Delta_C6_core']},"
            f"{row['relative_delta_pct']}"
        )
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()
