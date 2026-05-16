import argparse
import csv
from pathlib import Path

from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from run_tdhf_atom import export_tdhf_atom, write_channels


DEFAULT_ATOMS = ["Ne", "Ar", "Kr"]


def percent_error(value, reference):
    return 100.0 * (float(value) - float(reference)) / float(reference)


def load_references(path="reference_alpha_c6.csv"):
    refs = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            if row.get("alpha0_ref", "").strip() and row.get("C6_ref", "").strip():
                refs[row["atom"].strip()] = {
                    "alpha0_ref": float(row["alpha0_ref"]),
                    "C6_ref": float(row["C6_ref"]),
                    "source": row.get("source", ""),
                }
    return refs


def summarize_atom_result(atom, alpha0_ref, c6_ref, alpha_row):
    alpha0 = float(alpha_row["alpha0_au"])
    c6 = float(alpha_row["C6_self_au"])
    return {
        "atom": atom,
        "alpha0_ref": float(alpha0_ref),
        "alpha0_tdhf": alpha0,
        "alpha0_err": percent_error(alpha0, alpha0_ref),
        "C6_ref": float(c6_ref),
        "C6_tdhf": c6,
        "C6_err": percent_error(c6, c6_ref),
        "n_channels": int(alpha_row["n_channels"]),
    }


def run_atom(atom, basis, nstates, output_root):
    atom_dir = Path(output_root) / atom.lower()
    atom_dir.mkdir(parents=True, exist_ok=True)
    channels_path = atom_dir / f"{atom.lower()}_tdhf_channels.csv"
    alpha_path = atom_dir / "alpha_c6_table.csv"
    c6_path = atom_dir / "c6_table.csv"

    rows, _summary = export_tdhf_atom(atom, basis, nstates=nstates)
    write_channels(channels_path, rows)
    alpha_rows = build_alpha_rows(channels_path)
    write_alpha_rows(alpha_path, alpha_rows)
    write_c6_rows(c6_path, build_c6_rows(channels_path))
    return next(row for row in alpha_rows if row["atom"] == atom)


def write_summary(path, rows):
    fieldnames = [
        "atom",
        "alpha0_ref",
        "alpha0_tdhf",
        "alpha0_err",
        "C6_ref",
        "C6_tdhf",
        "C6_err",
        "n_channels",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_noble_gas_tdhf(atoms=None, basis="aug-cc-pVQZ", nstates=200, output_root="results", reference_path="reference_alpha_c6.csv"):
    refs = load_references(reference_path)
    rows = []
    for atom in atoms or DEFAULT_ATOMS:
        alpha_row = run_atom(atom, basis, nstates, output_root)
        rows.append(summarize_atom_result(atom, refs[atom]["alpha0_ref"], refs[atom]["C6_ref"], alpha_row))
    write_summary(Path(output_root) / "noble_gas_tdhf_summary.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run TDHF alpha0/C6 benchmark for noble gases.")
    parser.add_argument("--atom", action="append", choices=DEFAULT_ATOMS)
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--nstates", type=int, default=200)
    parser.add_argument("--output-root", default="results")
    args = parser.parse_args(argv)

    rows = run_noble_gas_tdhf(args.atom or DEFAULT_ATOMS, args.basis, args.nstates, args.output_root)
    print("atom,alpha0_ref,alpha0_tdhf,alpha0_err,C6_ref,C6_tdhf,C6_err")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['alpha0_ref']:.8f},"
            f"{row['alpha0_tdhf']:.8f},"
            f"{row['alpha0_err']:.6f},"
            f"{row['C6_ref']:.8f},"
            f"{row['C6_tdhf']:.8f},"
            f"{row['C6_err']:.6f}"
        )


if __name__ == "__main__":
    main()
