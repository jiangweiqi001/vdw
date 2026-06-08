import argparse
import csv
from pathlib import Path

from compute_core_sternheimer import compute_core_sternheimer_channels, write_channels as write_core_channels
from run_psp_rpa_atom import run_atom as run_psp_atom
from run_semicore_c6_validation import run_validation


REFERENCE_PAIR_C6 = {
    ("Sr", "Sr"): 3170.0,
}


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def parse_shells(value):
    if isinstance(value, set):
        return set(value)
    return {shell.strip() for shell in str(value).replace(",", ";").split(";") if shell.strip()}


def infer_xc_from_pseudo(pseudo_name):
    upper = str(pseudo_name).upper()
    if "PBE" in upper or "GGA" in upper:
        return "pbe"
    if "LDA" in upper or "PADE" in upper:
        return "lda"
    return "pbe"


def select_ready_target(rows, atom="Sr"):
    matches = [
        row for row in rows
        if row.get("atom") == atom and row.get("workflow_status") == "candidate_ready"
    ]
    if not matches:
        raise ValueError(f"No candidate_ready target found for {atom}.")
    if len(matches) > 1:
        raise ValueError(f"Multiple candidate_ready targets found for {atom}.")
    return matches[0]


def _psp_channels_path(output_root, atom, pseudo_name, basis_label, xc, method):
    leaf = f"{pseudo_name}_{basis_label}_{xc}_{method.lower()}".replace("/", "_")
    return Path(output_root) / "psp_rpa" / atom.lower() / leaf / f"{atom.lower()}_psp_channels.csv"


def validation_config_from_target(target, output_root="results/semicore_sr_validation", psp_nstates=100, method="TDDFT", core_basis="def2-TZVPP"):
    if target.get("workflow_status") != "candidate_ready":
        raise ValueError("Sr validation requires a candidate_ready workflow target.")
    atom = target["atom"]
    pseudo_name = target["pseudo_name"]
    basis_label = target["basis_label"]
    basis_name = target.get("basis_block_name") or basis_label
    xc = infer_xc_from_pseudo(pseudo_name)
    output_root = Path(output_root)
    return {
        "atom": atom,
        "pseudo_name": pseudo_name,
        "basis_label": basis_label,
        "basis_name": basis_name,
        "basis_file": target.get("basis_file", ""),
        "pseudo_file": target.get("pseudo_file", ""),
        "xc": xc,
        "method": method,
        "psp_nstates": int(psp_nstates),
        "active_electrons": int(target["active_electrons"]),
        "active_shells": target["active_shells"],
        "correction_shells": parse_shells(target["correction_shells"]),
        "core_basis": core_basis,
        "psp_output_root": output_root / "psp_rpa",
        "psp_channels_path": _psp_channels_path(output_root, atom, pseudo_name, basis_label, xc, method),
        "core_channels_path": output_root / "core_sternheimer" / f"{atom.lower()}_core_sternheimer_channels.csv",
        "validation_output_root": output_root / "validation",
        "reference_c6": REFERENCE_PAIR_C6.get((atom, atom)),
    }


def run_from_config(config):
    run_psp_atom(
        atom=config["atom"],
        psp=config["pseudo_name"],
        basis=config["basis_label"],
        xc=config["xc"],
        nstates=config["psp_nstates"],
        method=config["method"],
        output_root=config["psp_output_root"],
        basis_file=config["basis_file"],
        basis_name=config["basis_name"],
        pseudo_file=config["pseudo_file"],
        pseudo_name=config["pseudo_name"],
    )

    core_rows = compute_core_sternheimer_channels(
        atom=config["atom"],
        basis=config["core_basis"],
        selected_shells=config["correction_shells"],
    )
    write_core_channels(config["core_channels_path"], core_rows)

    return run_validation(
        atom=config["atom"],
        psp_channels=config["psp_channels_path"],
        core_channels=config["core_channels_path"],
        output_root=config["validation_output_root"],
        active_electrons=config["active_electrons"],
        active_shells=config["active_shells"],
        reference_c6=config["reference_c6"],
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Sr large-core PSP + semicore Sternheimer C6 validation.")
    parser.add_argument("--targets", required=True, help="Workflow target CSV from run_semicore_c6_workflow.py")
    parser.add_argument("--output-root", default="results/semicore_sr_validation")
    parser.add_argument("--psp-nstates", type=int, default=100)
    parser.add_argument("--method", choices=["TDDFT", "TDA"], default="TDDFT")
    parser.add_argument("--core-basis", default="def2-TZVPP")
    args = parser.parse_args(argv)

    target = select_ready_target(read_rows(args.targets), atom="Sr")
    config = validation_config_from_target(
        target,
        output_root=args.output_root,
        psp_nstates=args.psp_nstates,
        method=args.method,
        core_basis=args.core_basis,
    )
    summary = run_from_config(config)
    print(",".join(summary.keys()))
    print(",".join(str(value) for value in summary.values()))
    return summary


if __name__ == "__main__":
    main()
