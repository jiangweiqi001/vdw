import argparse
import csv
from pathlib import Path

from run_alpha_table import build_alpha_rows
from run_semicore_c6_validation import run_validation


SOURCE_SUFFIX = "STATIC_ALPHA_LOCAL_FIELD_CONSTRAINED"


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_rows(path, rows, fieldnames=None):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_alpha_reference(path, atom):
    rows = read_rows(path)
    for row in rows:
        if row["atom"] == atom:
            return row
    raise ValueError(f"No static polarizability reference for {atom} in {path}.")


def alpha_row(path, atom):
    return next(row for row in build_alpha_rows(path) if row["atom"] == atom)


def table_row(path, atom):
    return next(row for row in read_rows(path) if row["atom"] == atom)


def local_field_scale(alpha_target, alpha_psp, alpha_core):
    missing = float(alpha_target) - float(alpha_psp)
    if float(alpha_core) <= 0.0:
        raise ValueError("Core static polarizability must be positive.")
    return max(0.0, missing / float(alpha_core))


def scale_core_rows(core_rows, scale):
    scaled = []
    for row in core_rows:
        updated = dict(row)
        updated["osc"] = f"{float(row['osc']) * float(scale):.12f}"
        updated["local_field_scale"] = f"{float(scale):.12f}"
        updated["source"] = f"{row.get('source', '').strip()}_{SOURCE_SUFFIX}".strip("_")
        scaled.append(updated)
    return scaled


def run_alpha_constrained_validation(
    atom,
    psp_channels,
    core_channels,
    alpha_reference,
    output_root,
    active_electrons,
    active_shells,
    c6_reference=None,
):
    output_root = Path(output_root)
    psp_alpha = alpha_row(psp_channels, atom)
    core_alpha = alpha_row(core_channels, atom)
    alpha_ref = read_alpha_reference(alpha_reference, atom)
    scale = local_field_scale(alpha_ref["alpha0_ref"], psp_alpha["alpha0_au"], core_alpha["alpha0_au"])

    core_rows = read_rows(core_channels)
    scaled_rows = scale_core_rows(core_rows, scale)
    fieldnames = sorted(set().union(*(row.keys() for row in scaled_rows)))
    scaled_core_path = output_root / f"{atom.lower()}_alpha_constrained_cphf_channels.csv"
    write_rows(scaled_core_path, scaled_rows, fieldnames)

    validation = run_validation(
        atom=atom,
        psp_channels=psp_channels,
        core_channels=scaled_core_path,
        output_root=output_root / "validation",
        active_electrons=active_electrons,
        active_shells=active_shells,
        reference_c6=c6_reference,
    )
    corrected_alpha = table_row(output_root / "validation" / "psp_plus_sternheimer_alpha_c6_table.csv", atom)
    summary = {
        "atom": atom,
        "alpha0_psp": psp_alpha["alpha0_au"],
        "alpha0_core_unconstrained": core_alpha["alpha0_au"],
        "alpha0_reference": alpha_ref["alpha0_ref"],
        "alpha0_uncertainty": alpha_ref.get("uncertainty", ""),
        "local_field_scale": scale,
        "alpha0_corrected": corrected_alpha["alpha0_au"],
        "C6_PSP": validation["C6_PSP"],
        "C6_PSP_plus_alpha_constrained_CPHF": validation["C6_PSP_plus_sternheimer"],
        "C6_reference": validation["C6_reference"],
        "reference_error_pct": validation["reference_error_pct"],
        "go_no_go": validation["go_no_go"],
        "scaled_core_channels": str(scaled_core_path),
        "alpha_reference_source": alpha_ref.get("source", ""),
    }
    write_rows(output_root / "kernel_summary.csv", [summary], list(summary.keys()))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Apply static-alpha constrained local-field scaling to radial CPHF core channels.")
    parser.add_argument("--atom", required=True, choices=["Zn", "Cd"])
    parser.add_argument("--psp-channels", required=True)
    parser.add_argument("--core-channels", required=True)
    parser.add_argument("--alpha-reference", default="reference_static_polarizability_group12.csv")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--active-electrons", required=True, type=int)
    parser.add_argument("--active-shells", required=True)
    parser.add_argument("--c6-reference")
    args = parser.parse_args(argv)

    summary = run_alpha_constrained_validation(
        atom=args.atom,
        psp_channels=args.psp_channels,
        core_channels=args.core_channels,
        alpha_reference=args.alpha_reference,
        output_root=args.output_root,
        active_electrons=args.active_electrons,
        active_shells=args.active_shells,
        c6_reference=args.c6_reference,
    )
    print(",".join(summary.keys()))
    print(",".join(str(value) for value in summary.values()))
    return summary


if __name__ == "__main__":
    main()
