import argparse
import csv
from pathlib import Path


DEFAULT_ATOMS = ["Ne", "Ar", "Kr", "Mg", "Ca"]
DEFAULT_PSEUDOS = ["gth-lda", "gth-pbe", "gth-hf"]
DEFAULT_BASES = [
    "gth-szv",
    "gth-dzv",
    "gth-dzvp",
    "gth-tzv2p",
    "gth-tzvp",
    "gth-qzv2p",
    "gth-qzvp",
    "gth-szv-molopt-sr",
    "gth-dzvp-molopt-sr",
    "gth-tzvp-molopt-sr",
    "gth-tzv2p-molopt-sr",
    "gth-qzv2p-molopt-sr",
]


def probe_case(atom, pseudo, basis):
    from pyscf import gto

    try:
        mol = gto.M(
            atom=f"{atom} 0 0 0",
            basis=basis,
            pseudo=pseudo,
            spin=0,
            charge=0,
            cart=False,
            verbose=0,
        )
        return {
            "atom": atom,
            "pseudo": pseudo,
            "basis": basis,
            "available": True,
            "nelectron": mol.nelectron,
            "nao": mol.nao_nr(),
            "error": "",
        }
    except Exception as exc:
        return {
            "atom": atom,
            "pseudo": pseudo,
            "basis": basis,
            "available": False,
            "nelectron": "",
            "nao": "",
            "error": f"{type(exc).__name__}: {str(exc).splitlines()[0]}",
        }


def run_probe(atoms=None, pseudos=None, bases=None):
    rows = []
    for atom in atoms or DEFAULT_ATOMS:
        for pseudo in pseudos or DEFAULT_PSEUDOS:
            for basis in bases or DEFAULT_BASES:
                rows.append(probe_case(atom, pseudo, basis))
    return rows


def write_rows(path, rows):
    fieldnames = ["atom", "pseudo", "basis", "available", "nelectron", "nao", "error"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Probe PySCF GTH pseudopotential basis availability.")
    parser.add_argument("--atom", action="append")
    parser.add_argument("--pseudo", action="append")
    parser.add_argument("--basis", action="append")
    parser.add_argument("--output", default="results/psp_basis_availability.csv")
    args = parser.parse_args(argv)

    rows = run_probe(args.atom or DEFAULT_ATOMS, args.pseudo or DEFAULT_PSEUDOS, args.basis or DEFAULT_BASES)
    write_rows(args.output, rows)
    print("atom,pseudo,basis,available,nelectron,nao,error")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['pseudo']},"
            f"{row['basis']},"
            f"{str(row['available']).lower()},"
            f"{row['nelectron']},"
            f"{row['nao']},"
            f"{row['error']}"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
