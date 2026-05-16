import argparse
import csv
import re
import shutil
from pathlib import Path

from pyscf_export_ar_radials import grouped_shells, make_ar_molecule
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows


AR_VALENCE_SHELLS = {"3s", "3p"}
AR_CORE_SHELLS = {"1s", "2s", "2p"}


def mo_shell_map_for_ar(basis="aug-cc-pVQZ"):
    from pyscf.scf import atom_hf

    mol = make_ar_molecule(basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"Ar atom HF did not converge for {basis}.")
    shells = grouped_shells(mol, mf.mo_energy, mf.mo_occ, mf.mo_coeff)
    mapping = {}
    for shell in shells:
        for column in shell["columns"]:
            mapping[int(column)] = shell["label"]
    return mapping


def occupied_mo_index(channel):
    match = re.search(r"mo_occ_(\d+)_to_virt_", channel)
    if not match:
        return None
    return int(match.group(1))


def split_rows(rows, mo_to_shell, valence_shells=AR_VALENCE_SHELLS, core_shells=AR_CORE_SHELLS):
    valence = []
    core = []
    for row in rows:
        occ_idx = occupied_mo_index(row["channel"])
        shell = mo_to_shell.get(occ_idx)
        enriched = {**row, "from_orbital": shell or ""}
        if shell in valence_shells:
            valence.append(enriched)
        elif shell in core_shells:
            core.append(enriched)
    return valence, core


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_channel_rows(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "channel", "delta_Ha", "osc", "is_core", "source"])
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in writer.fieldnames})


def alpha_c6_for_channels(path, output_prefix):
    alpha_path = Path(f"{output_prefix}_alpha_c6_table.csv")
    c6_path = Path(f"{output_prefix}_c6_table.csv")
    alpha_rows = build_alpha_rows(path)
    write_alpha_rows(alpha_path, alpha_rows)
    write_c6_rows(c6_path, build_c6_rows(path))
    return next(row for row in alpha_rows if row["atom"] == "Ar")


def decompose_mo_channels(input_path, output_dir, basis="aug-cc-pVQZ"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_path = output_dir / "ar_mo_all_channels.csv"
    valence_path = output_dir / "ar_mo_valence_channels.csv"
    core_path = output_dir / "ar_mo_core_channels.csv"

    rows = read_rows(input_path)
    mo_to_shell = mo_shell_map_for_ar(basis)
    valence_rows, core_rows = split_rows(rows, mo_to_shell)
    if Path(input_path).resolve() != all_path.resolve():
        shutil.copyfile(input_path, all_path)
    write_channel_rows(valence_path, valence_rows)
    write_channel_rows(core_path, core_rows)

    summaries = {
        "all": alpha_c6_for_channels(all_path, output_dir / "all"),
        "valence": alpha_c6_for_channels(valence_path, output_dir / "valence"),
        "core": alpha_c6_for_channels(core_path, output_dir / "core"),
    }
    c6_all = float(summaries["all"]["C6_self_au"])
    c6_val = float(summaries["valence"]["C6_self_au"])
    delta_core = c6_all - c6_val
    decomp_path = output_dir / "ar_mo_core_valence_decomposition.csv"
    with open(decomp_path, "w", newline="", encoding="utf-8") as fp:
        fieldnames = ["method", "component", "alpha0", "C6", "C6_error_pct", "n_channels", "note"]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for component, row in summaries.items():
            writer.writerow(
                {
                    "method": "MO",
                    "component": component,
                    "alpha0": row["alpha0_au"],
                    "C6": row["C6_self_au"],
                    "C6_error_pct": "",
                    "n_channels": row["n_channels"],
                    "note": "independent-particle channel split",
                }
            )
        writer.writerow(
            {
                "method": "MO",
                "component": "Delta_C6_core",
                "alpha0": "",
                "C6": f"{delta_core:.8f}",
                "C6_error_pct": "",
                "n_channels": "",
                "note": f"relative_core_contribution={delta_core / c6_all:.8f}",
            }
        )
    return summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description="Split Ar MO oscillator channels by occupied shell.")
    parser.add_argument("--input", default="results/ar/mo/aug-cc-pvqz/ar_mo_channels.csv")
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--output-dir", default="results/ar/core_valence_mo")
    args = parser.parse_args(argv)

    summaries = decompose_mo_channels(args.input, args.output_dir, args.basis)
    print("component,alpha0,C6,n_channels")
    for component, row in summaries.items():
        print(f"{component},{row['alpha0_au']},{row['C6_self_au']},{row['n_channels']}")


if __name__ == "__main__":
    main()
