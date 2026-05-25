import argparse
import csv
from pathlib import Path

from compute_multipole_core_wilson import compute_multipole_core_channels, write_channels as write_core_channels
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_all_e_rpa_atom import response_rows_from_arrays, write_channels as write_response_channels
from run_c6_table import build_c6_rows, write_c6_rows
from run_eft_core_dipole_validation import write_combined
from run_psp_rpa_atom import run_atom as run_psp_atom
from run_psp_rpa_atom import write_summary as write_psp_summary


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def c6_for(path, prefix, atom="Na"):
    alpha_rows = build_alpha_rows(path)
    write_alpha_rows(f"{prefix}_alpha_c6_table.csv", alpha_rows)
    write_c6_rows(f"{prefix}_c6_table.csv", build_c6_rows(path))
    return next(row for row in alpha_rows if row["atom"] == atom)


def run_all_e_na(output_dir, basis="aug-cc-pVQZ", nstates=100):
    from pyscf import gto, dft, tdscf

    mol = gto.M(atom="Na 0 0 0", basis=basis, spin=1, charge=0, cart=False, verbose=0)
    mf = dft.UKS(mol)
    mf.xc = "pbe,pbe"
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("Na all-electron UKS-PBE did not converge.")
    td = tdscf.TDDFT(mf)
    td.nstates = nstates
    td.verbose = 0
    energies = td.kernel()[0]
    oscillator_strengths = td.oscillator_strength()
    rows = response_rows_from_arrays("Na", energies, oscillator_strengths, "PySCF_TDDFT_PBE", min_osc=1e-10)
    path = Path(output_dir) / "na_all_e_channels.csv"
    write_response_channels(path, rows)
    return c6_for(path, Path(output_dir) / "all_e", atom="Na")


def summary_row(psp, eft, all_e):
    c6_psp = float(psp["C6_self_au"])
    c6_eft = float(eft["C6_self_au"])
    c6_all = float(all_e["C6_self_au"])
    missing = c6_all - c6_psp
    correction = c6_eft - c6_psp
    return {
        "case": "Na_q1_adapted",
        "C6_PSP": c6_psp,
        "C6_PSP_EFT": c6_eft,
        "C6_all_e": c6_all,
        "Delta_C6_missing": missing,
        "Delta_C6_EFT": correction,
        "closure_fraction": correction / missing if missing else 0.0,
        "closure_pct": 100.0 * correction / missing if missing else 0.0,
        "double_counting_status": "clean",
        "note": "adapted diagnostic: q1 pseudo with q1 UZH basis; open-shell UKS TDDFT",
    }


def write_summary(path, rows):
    fieldnames = [
        "case",
        "C6_PSP",
        "C6_PSP_EFT",
        "C6_all_e",
        "Delta_C6_missing",
        "Delta_C6_EFT",
        "closure_fraction",
        "closure_pct",
        "double_counting_status",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_benchmark(output_dir="results/na_q1"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    psp_summary = run_psp_atom(
        atom="Na",
        psp="GTH-PBE-q1",
        basis="TZV2P-MOLOPT-PBE-GTH-q1",
        xc="pbe",
        nstates=100,
        method="TDDFT",
        output_root="results/psp_rpa",
        basis_file="external_data/cp2k/BASIS_MOLOPT_UZH",
        basis_name="TZV2P-MOLOPT-PBE-GTH-q1",
        pseudo_file="external_data/cp2k/GTH_POTENTIALS",
        pseudo_name="GTH-PBE-q1",
        spin=1,
    )
    write_psp_summary(output_dir / "psp_summary.csv", [psp_summary])
    psp_path = "results/psp_rpa/na/GTH-PBE-q1_TZV2P-MOLOPT-PBE-GTH-q1_pbe_tddft/na_psp_channels.csv"
    psp_c6 = c6_for(psp_path, output_dir / "psp", atom="Na")

    core_rows, _core_summary = compute_multipole_core_channels("Na", 1, "aug-cc-pVQZ", 100)
    core_path = output_dir / "na_core_multipole_channels.csv"
    write_core_channels(core_path, core_rows)
    combined_path = output_dir / "na_psp_plus_eft_channels.csv"
    write_combined(combined_path, read_rows(psp_path) + core_rows)
    eft_c6 = c6_for(combined_path, output_dir / "psp_plus_eft", atom="Na")
    all_e_c6 = run_all_e_na(output_dir)
    row = summary_row(psp_c6, eft_c6, all_e_c6)
    write_summary(output_dir / "summary.csv", [row])
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Na q1 adapted clean PSP+EFT benchmark.")
    parser.add_argument("--output-dir", default="results/na_q1")
    args = parser.parse_args(argv)
    row = run_benchmark(args.output_dir)
    print("case,C6_PSP,C6_PSP_EFT,C6_all_e,closure_pct,double_counting_status,note")
    print(
        f"{row['case']},{row['C6_PSP']:.8f},{row['C6_PSP_EFT']:.8f},"
        f"{row['C6_all_e']:.8f},{row['closure_pct']:.6f},"
        f"{row['double_counting_status']},{row['note']}"
    )


if __name__ == "__main__":
    main()
