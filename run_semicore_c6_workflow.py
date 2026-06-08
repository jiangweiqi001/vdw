import argparse
import csv
from pathlib import Path

from semicore_c6_targets import SEMICORE_C6_TARGETS, audit_semicore_target


TARGET_ATOMS = ["Sr", "Zn", "Cd"]


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "dimer",
        "workflow_status",
        "pseudo_name",
        "basis_label",
        "basis_block_name",
        "basis_file",
        "pseudo_file",
        "candidate_status",
        "active_electrons",
        "active_shells",
        "expected_active_shells",
        "correction_shells",
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
    pseudo = row.get("pseudo_name", "")
    basis = row.get("basis_label") or row.get("basis_name", "")
    pseudo_order = ["GTH-PBE-q2", "GTH-LDA-q2", "GTH-PADE-q2", "GTH-BLYP-q2", "GTH-BP-q2"]
    pseudo_rank = pseudo_order.index(pseudo) if pseudo in pseudo_order else len(pseudo_order)
    basis_rank = 0
    upper_basis = basis.upper()
    for idx, token in enumerate(["QZ", "TZV2P", "TZVP", "DZVP", "SZV"]):
        if token in upper_basis:
            basis_rank = idx
            break
    else:
        basis_rank = 99
    return (pseudo_rank, basis_rank, basis)


def _audit_candidate(atom, row):
    return audit_semicore_target(
        atom,
        active_electrons=row.get("active_electrons", ""),
        active_shells=row.get("active_shells", ""),
    )


def select_runnable_candidate(atom, candidate_rows):
    """Select a TDDFT-runnable large-core candidate with no semicore double counting."""
    matches = []
    for row in candidate_rows:
        if row.get("atom") != atom:
            continue
        if row.get("candidate_status") != "tddft_smoke_ok":
            continue
        try:
            audit = _audit_candidate(atom, row)
        except Exception:
            continue
        if audit["audit_status"] != "pass":
            continue
        matches.append(row)
    if not matches:
        return None
    return sorted(matches, key=_candidate_rank)[0]


def _candidate_notes(atom, candidate_rows):
    notes = []
    for row in candidate_rows:
        if row.get("atom") != atom:
            continue
        status = row.get("candidate_status", "")
        note = row.get("note", "")
        label = row.get("basis_label") or row.get("basis_name", "")
        parts = [part for part in [status, row.get("pseudo_name", ""), label, note] if part]
        if parts:
            notes.append(" / ".join(parts))
    return " | ".join(notes)


def _blocked_by_environment(atom, candidate_rows):
    atom_rows = [row for row in candidate_rows if row.get("atom") == atom]
    if not atom_rows:
        return False
    return any("No module named 'pyscf'" in row.get("note", "") for row in atom_rows)


def unavailable_target_row(atom, candidate_rows):
    target = SEMICORE_C6_TARGETS[atom]
    return {
        "atom": atom,
        "dimer": target.dimer,
        "workflow_status": "environment_blocked" if _blocked_by_environment(atom, candidate_rows) else "unavailable",
        "pseudo_name": "",
        "basis_label": "",
        "basis_block_name": "",
        "basis_file": "",
        "pseudo_file": "",
        "candidate_status": "",
        "active_electrons": "",
        "active_shells": "",
        "expected_active_shells": ";".join(sorted(target.explicit_shells)),
        "correction_shells": ";".join(sorted(target.correction_shells)),
        "target_audit_status": "not_run",
        "note": _candidate_notes(atom, candidate_rows) or "no candidate rows for target",
    }


def candidate_ready_row(atom, candidate):
    audit = _audit_candidate(atom, candidate)
    return {
        "atom": atom,
        "dimer": audit["dimer"],
        "workflow_status": "candidate_ready",
        "pseudo_name": candidate.get("pseudo_name", ""),
        "basis_label": candidate.get("basis_label") or candidate.get("basis_name", ""),
        "basis_block_name": candidate.get("basis_block_name", ""),
        "basis_file": candidate.get("basis_file", ""),
        "pseudo_file": candidate.get("pseudo_file", ""),
        "candidate_status": candidate.get("candidate_status", ""),
        "active_electrons": audit["active_electrons"],
        "active_shells": audit["active_shells"],
        "expected_active_shells": audit["expected_active_shells"],
        "correction_shells": audit["correction_shells"],
        "target_audit_status": audit["audit_status"],
        "note": candidate.get("note", ""),
    }


def build_workflow_rows(candidate_rows, atoms=None):
    rows = []
    for atom in atoms or TARGET_ATOMS:
        candidate = select_runnable_candidate(atom, candidate_rows)
        if candidate is None:
            rows.append(unavailable_target_row(atom, candidate_rows))
        else:
            rows.append(candidate_ready_row(atom, candidate))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Prepare clean Sr/Zn/Cd large-core semicore C6 workflow targets.")
    parser.add_argument("--candidate-scan", required=True)
    parser.add_argument("--atom", action="append", choices=TARGET_ATOMS)
    parser.add_argument("--output", default="results/semicore_c6_workflow_targets.csv")
    args = parser.parse_args(argv)

    rows = build_workflow_rows(read_rows(args.candidate_scan), atoms=args.atom or TARGET_ATOMS)
    write_rows(args.output, rows)
    print("atom,workflow_status,pseudo_name,basis_label,target_audit_status,note")
    for row in rows:
        print(
            f"{row['atom']},{row['workflow_status']},{row['pseudo_name']},"
            f"{row['basis_label']},{row['target_audit_status']},{row['note']}"
        )
    return rows


if __name__ == "__main__":
    main()
