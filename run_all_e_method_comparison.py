import argparse
import csv
from pathlib import Path

from run_noble_gas_tdhf import load_references, percent_error


ATOM_BASIS = {
    "Ne": "aug-cc-pVQZ",
    "Ar": "aug-cc-pVQZ",
    "Kr": "aug-cc-pVQZ",
    "Mg": "aug-cc-pVQZ",
    "Ca": "cc-pVQZ",
}


def _read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def _all_e_by_atom_xc(path="results/all_e_rpa_summary.csv"):
    rows = {}
    for row in _read_rows(path):
        if int(row["nstates"]) != 200 or row["method"] != "TDDFT":
            continue
        rows[(row["atom"], row["xc"])] = row
    return rows


def _tdhf_reference_rows():
    rows = {}
    for row in _read_rows("results/noble_gas_tdhf_summary.csv"):
        rows[row["atom"]] = {
            "alpha0": float(row["alpha0_tdhf"]),
            "C6": float(row["C6_tdhf"]),
            "source": "HF_TDHF",
        }
    for atom, path in {
        "Mg": "results/mg/mg_core_valence_decomposition.csv",
        "Ca": "results/ca/ca_core_valence_decomposition.csv",
    }.items():
        if not Path(path).exists():
            continue
        all_row = next(row for row in _read_rows(path) if row["component"] == "all")
        rows[atom] = {
            "alpha0": float(all_row["alpha0"]),
            "C6": float(all_row["C6"]),
            "source": "HF_TDHF_projected_all",
        }
    return rows


def comparison_row(atom, xc, ks_alpha0, ks_c6, hf_alpha0, hf_c6, ref_alpha0, ref_c6, basis):
    return {
        "atom": atom,
        "xc": xc,
        "basis": basis,
        "nstates": 200,
        "KS_alpha0": float(ks_alpha0),
        "KS_C6": float(ks_c6),
        "KS_alpha0_err": percent_error(ks_alpha0, ref_alpha0),
        "KS_C6_err": percent_error(ks_c6, ref_c6),
        "HF_TDHF_alpha0": float(hf_alpha0),
        "HF_TDHF_C6": float(hf_c6),
        "HF_TDHF_alpha0_err": percent_error(hf_alpha0, ref_alpha0),
        "HF_TDHF_C6_err": percent_error(hf_c6, ref_c6),
        "delta_alpha0_KS_minus_HF": float(ks_alpha0) - float(hf_alpha0),
        "delta_C6_KS_minus_HF": float(ks_c6) - float(hf_c6),
    }


def build_comparison_rows():
    refs = load_references()
    all_e = _all_e_by_atom_xc()
    tdhf = _tdhf_reference_rows()
    rows = []
    for atom in ["Ne", "Ar", "Kr", "Mg", "Ca"]:
        if atom not in tdhf:
            continue
        for xc in ["lda", "pbe"]:
            ks = all_e.get((atom, xc))
            if not ks:
                continue
            ref = refs[atom]
            rows.append(
                comparison_row(
                    atom=atom,
                    xc=xc,
                    basis=ks["basis"],
                    ks_alpha0=float(ks["alpha0"]),
                    ks_c6=float(ks["C6"]),
                    hf_alpha0=tdhf[atom]["alpha0"],
                    hf_c6=tdhf[atom]["C6"],
                    ref_alpha0=ref["alpha0_ref"],
                    ref_c6=ref["C6_ref"],
                )
            )
    return rows


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "xc",
        "basis",
        "nstates",
        "KS_alpha0",
        "KS_C6",
        "KS_alpha0_err",
        "KS_C6_err",
        "HF_TDHF_alpha0",
        "HF_TDHF_C6",
        "HF_TDHF_alpha0_err",
        "HF_TDHF_C6_err",
        "delta_alpha0_KS_minus_HF",
        "delta_C6_KS_minus_HF",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare all-electron KS TDDFT against HF/TDHF route.")
    parser.add_argument("--output", default="results/all_e_rpa_vs_hf_tdhf_summary.csv")
    args = parser.parse_args(argv)

    rows = build_comparison_rows()
    write_rows(args.output, rows)
    print("atom,xc,basis,KS_C6,KS_C6_err,HF_TDHF_C6,HF_TDHF_C6_err,delta_C6_KS_minus_HF")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['xc']},"
            f"{row['basis']},"
            f"{row['KS_C6']:.8f},"
            f"{row['KS_C6_err']:.6f},"
            f"{row['HF_TDHF_C6']:.8f},"
            f"{row['HF_TDHF_C6_err']:.6f},"
            f"{row['delta_C6_KS_minus_HF']:.8f}"
        )


if __name__ == "__main__":
    main()
