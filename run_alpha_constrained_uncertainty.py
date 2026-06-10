import argparse
import csv
from pathlib import Path

import numpy as np

from eft_alpha import alpha0_from_osc, c6_from_alpha, load_channels_csv
from run_radial_cphf_alpha_constrained_validation import local_field_scale


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_rows(path, rows):
    fieldnames = list(rows[0].keys()) if rows else []
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def alpha_reference_row(rows, atom):
    for row in rows:
        if row["atom"] == atom:
            return row
    raise ValueError(f"No static polarizability reference for {atom}.")


def pair_reference_rows(rows, atom):
    return [row for row in rows if row["A"] == atom and row["B"] == atom]


def corrected_c6_for_alpha_target(atom, psp_channels, core_channels, alpha_target):
    psp = load_channels_csv(psp_channels)[atom]
    core = load_channels_csv(core_channels)[atom]
    alpha_psp = alpha0_from_osc(psp["delta"], psp["osc"])
    alpha_core = alpha0_from_osc(core["delta"], core["osc"])
    scale = local_field_scale(alpha_target, alpha_psp, alpha_core)
    delta = np.concatenate([psp["delta"], core["delta"]])
    osc = np.concatenate([psp["osc"], core["osc"] * scale])
    return {
        "alpha0_psp": float(alpha_psp),
        "alpha0_core_unconstrained": float(alpha_core),
        "alpha0_corrected": float(alpha_psp + alpha_core * scale),
        "local_field_scale": float(scale),
        "C6_model": float(c6_from_alpha(delta, osc, delta, osc)),
    }


def alpha_targets(alpha_ref, uncertainty):
    alpha_ref = float(alpha_ref)
    uncertainty = float(uncertainty)
    return [
        ("low", alpha_ref - uncertainty),
        ("center", alpha_ref),
        ("high", alpha_ref + uncertainty),
    ]


def build_uncertainty_rows(
    atom,
    psp_channels,
    core_channels,
    alpha_references,
    pair_references,
    tolerance_pct=10.0,
):
    alpha_ref = alpha_reference_row(read_rows(alpha_references), atom)
    pair_refs = pair_reference_rows(read_rows(pair_references), atom)
    rows = []
    for alpha_case, alpha_target in alpha_targets(alpha_ref["alpha0_ref"], alpha_ref["uncertainty"]):
        model = corrected_c6_for_alpha_target(atom, psp_channels, core_channels, alpha_target)
        for ref in pair_refs:
            c6_ref = float(ref["C6_ref"])
            error_pct = 100.0 * (model["C6_model"] - c6_ref) / c6_ref
            rows.append(
                {
                    "atom": atom,
                    "alpha_case": alpha_case,
                    "alpha0_target": f"{alpha_target:.8f}",
                    "alpha0_psp": f"{model['alpha0_psp']:.8f}",
                    "alpha0_core_unconstrained": f"{model['alpha0_core_unconstrained']:.8f}",
                    "alpha0_corrected": f"{model['alpha0_corrected']:.8f}",
                    "local_field_scale": f"{model['local_field_scale']:.8f}",
                    "C6_model": f"{model['C6_model']:.8f}",
                    "C6_ref": f"{c6_ref:.8f}",
                    "reference_label": ref["reference_label"],
                    "error_pct": f"{error_pct:.8f}",
                    "within_tolerance": str(abs(error_pct) <= float(tolerance_pct)).lower(),
                    "tolerance_pct": f"{float(tolerance_pct):.8f}",
                    "alpha_reference_source": alpha_ref.get("source", ""),
                    "c6_reference_source": ref.get("source", ""),
                }
            )
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Propagate static-alpha uncertainty through alpha-constrained radial CPHF C6.")
    parser.add_argument("--atom", action="append", required=True, choices=["Zn", "Cd"])
    parser.add_argument("--psp-channels", action="append", required=True)
    parser.add_argument("--core-channels", action="append", required=True)
    parser.add_argument("--alpha-references", default="reference_static_polarizability_group12.csv")
    parser.add_argument("--pair-references", default="reference_pair_c6_alternates.csv")
    parser.add_argument("--tolerance-pct", type=float, default=10.0)
    parser.add_argument("--output", default="results/radial_grid_sternheimer/alpha_constrained/uncertainty_summary.csv")
    args = parser.parse_args(argv)

    if not (len(args.atom) == len(args.psp_channels) == len(args.core_channels)):
        raise ValueError("--atom, --psp-channels, and --core-channels must have the same length.")

    rows = []
    for atom, psp_channels, core_channels in zip(args.atom, args.psp_channels, args.core_channels):
        rows.extend(
            build_uncertainty_rows(
                atom=atom,
                psp_channels=psp_channels,
                core_channels=core_channels,
                alpha_references=args.alpha_references,
                pair_references=args.pair_references,
                tolerance_pct=args.tolerance_pct,
            )
        )
    write_rows(args.output, rows)
    print("atom,alpha_case,C6_model,C6_ref,reference_label,error_pct,within_tolerance")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['alpha_case']},"
            f"{row['C6_model']},"
            f"{row['C6_ref']},"
            f"{row['reference_label']},"
            f"{row['error_pct']},"
            f"{row['within_tolerance']}"
        )
    print(f"\nWrote {args.output}")
    return rows


if __name__ == "__main__":
    main()
