import argparse
import csv
from pathlib import Path

from analyze_channels import analyze_channels
from build_eft_channels_spectral import build_spectral_channels_from_csv, write_channels
from build_radial_dipoles import build_radial_dipoles_from_orbitals, write_radial_dipoles
from check_radial_orbitals import check_radial_orbitals
from compare_alpha_c6 import compare_tables
from convert_atomic_solver_output import convert_rows, read_input, write_radial_orbitals
from eft_alpha import alpha0_from_osc, c6_from_alpha, load_channels_csv, self_c6_from_osc


BASIS_LIST = ["cc-pVTZ", "aug-cc-pVTZ", "d-aug-cc-pVTZ-local", "aug-cc-pVQZ"]


def safe_basis_name(basis):
    return basis.lower().replace("+", "p").replace("*", "s")


def write_alpha_table(path, channels_path):
    data = load_channels_csv(channels_path)
    rows = []
    for atom in sorted(data):
        delta = data[atom]["delta"]
        osc = data[atom]["osc"]
        rows.append(
            {
                "atom": atom,
                "alpha0_au": f"{alpha0_from_osc(delta, osc):.8f}",
                "C6_self_au": f"{self_c6_from_osc(delta, osc):.8f}",
                "n_channels": len(delta),
            }
        )
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["atom", "alpha0_au", "C6_self_au", "n_channels"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_c6_table(path, channels_path):
    data = load_channels_csv(channels_path)
    rows = []
    for atom in sorted(data):
        c6 = float(c6_from_alpha(data[atom]["delta"], data[atom]["osc"], data[atom]["delta"], data[atom]["osc"]))
        rows.append({"A": atom, "B": atom, "C6_au": f"{c6:.8f}"})
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["A", "B", "C6_au"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_channel_analysis(path, rows):
    fieldnames = [
        "atom",
        "channel",
        "delta_Ha",
        "osc",
        "alpha0_contribution",
        "alpha0_fraction",
        "single_channel_c6",
        "single_channel_c6_fraction",
        "cross_inclusive_c6",
        "cross_inclusive_c6_fraction",
        "is_residual",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def core_electron_count(spectral_input_path):
    total = 0.0
    with open(spectral_input_path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            if row["atom"] == "Ar" and row["type"].strip().lower() == "core":
                total += float(row["occupation"])
    return total


def oscillator_sums(channel_rows):
    discrete = sum(float(row["osc"]) for row in channel_rows if not row["is_residual"])
    residual = sum(float(row["osc"]) for row in channel_rows if row["is_residual"])
    return discrete, residual, discrete + residual


def summarize_basis_result(basis, alpha_rows, channel_rows, n_core=0.0):
    ar_alpha = next(row for row in alpha_rows if row["atom"] == "Ar")
    channel_3p_3d = next((row for row in channel_rows if row["atom"] == "Ar" and row["channel"] == "3p_to_3d"), None)
    residual = next((row for row in channel_rows if row["atom"] == "Ar" and row["channel"] == "missing_core_continuum"), None)
    sum_osc_discrete, sum_osc_residual, sum_osc_total = oscillator_sums(channel_rows)
    discrete_ratio = sum_osc_discrete / n_core if n_core else 0.0
    total_ratio = sum_osc_total / n_core if n_core else 0.0
    invalid_for_prediction = discrete_ratio > 1.1
    invalid_reason = (
        f"sum_osc_discrete/N_core={discrete_ratio:.6f}"
        if invalid_for_prediction
        else ""
    )
    return {
        "basis": basis,
        "alpha0": float(ar_alpha["alpha0_au"]),
        "C6": float(ar_alpha["C6_self_au"]),
        "n_channels": int(ar_alpha["n_channels"]),
        "N_core": float(n_core),
        "sum_osc_discrete": sum_osc_discrete,
        "sum_osc_residual": sum_osc_residual,
        "sum_osc_total": sum_osc_total,
        "sum_osc_discrete_over_N_core": discrete_ratio,
        "sum_osc_total_over_N_core": total_ratio,
        "invalid_for_prediction": invalid_for_prediction,
        "invalid_reason": invalid_reason,
        "3p_to_3d_delta": float(channel_3p_3d["delta_Ha"]) if channel_3p_3d else 0.0,
        "3p_to_3d_osc": float(channel_3p_3d["osc"]) if channel_3p_3d else 0.0,
        "3p_to_3d_alpha_fraction": float(channel_3p_3d["alpha0_fraction"]) if channel_3p_3d else 0.0,
        "3p_to_3d_c6_fraction": float(channel_3p_3d["cross_inclusive_c6_fraction"]) if channel_3p_3d else 0.0,
        "residual_fraction": float(residual["alpha0_fraction"]) if residual else 0.0,
        "residual_c6_fraction": float(residual["cross_inclusive_c6_fraction"]) if residual else 0.0,
    }


def write_basis_summary(path, rows):
    fieldnames = [
        "basis",
        "alpha0",
        "C6",
        "n_channels",
        "N_core",
        "sum_osc_discrete",
        "sum_osc_residual",
        "sum_osc_total",
        "sum_osc_discrete_over_N_core",
        "sum_osc_total_over_N_core",
        "invalid_for_prediction",
        "invalid_reason",
        "3p_to_3d_delta",
        "3p_to_3d_osc",
        "3p_to_3d_alpha_fraction",
        "3p_to_3d_c6_fraction",
        "residual_fraction",
        "residual_c6_fraction",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_one_basis(basis, output_root, n_grid=2400):
    from pyscf_export_ar_radials import export_ar_radials

    result_dir = output_root / safe_basis_name(basis)
    result_dir.mkdir(parents=True, exist_ok=True)

    solver_output = result_dir / "atomic_solver_output.csv"
    spectral_input = result_dir / "atomic_spectral_input.csv"
    radial_orbitals = result_dir / "radial_orbitals.csv"
    radial_dipoles = result_dir / "radial_dipoles.csv"
    atomic_channels = result_dir / "atomic_channels.csv"
    alpha_table = result_dir / "alpha_c6_table.csv"
    c6_table = result_dir / "c6_table.csv"
    channel_table = result_dir / "channel_analysis.csv"

    export_ar_radials(
        solver_output_path=solver_output,
        spectral_input_path=spectral_input,
        basis=basis,
        n_grid=n_grid,
    )
    write_radial_orbitals(radial_orbitals, convert_rows(read_input(solver_output), input_kind="R", normalize=True))
    radial_report = check_radial_orbitals(radial_orbitals)
    if any(not row["norm_ok"] or not row["monotonic_grid"] for row in radial_report["orbitals"]):
        raise RuntimeError(f"Radial orbital validation failed for {basis}.")

    write_radial_dipoles(radial_dipoles, build_radial_dipoles_from_orbitals(spectral_input, radial_orbitals))
    write_channels(
        atomic_channels,
        build_spectral_channels_from_csv(
            spectral_input,
            radial_dipoles,
            add_residual_oscillator=True,
            residual_path="residual_oscillators.csv",
        ),
    )
    alpha_rows = write_alpha_table(alpha_table, atomic_channels)
    write_c6_table(c6_table, atomic_channels)
    channel_rows = analyze_channels(atomic_channels)
    write_channel_analysis(channel_table, channel_rows)
    compare_rows = compare_tables(alpha_table, "reference_alpha_c6.csv")
    with open(result_dir / "compare_alpha_c6.csv", "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["atom", "alpha0_eft", "alpha0_ref", "err_alpha_pct", "C6_eft", "C6_ref", "err_C6_pct", "source"],
        )
        writer.writeheader()
        writer.writerows(compare_rows)
    return summarize_basis_result(basis, alpha_rows, channel_rows, core_electron_count(spectral_input)), channel_rows


def run_basis_convergence(bases=None, output_root=Path("results/ar"), n_grid=2400):
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    summaries = []
    all_channels = []
    for basis in bases or BASIS_LIST:
        summary, channel_rows = run_one_basis(basis, output_root, n_grid=n_grid)
        summaries.append(summary)
        for row in channel_rows:
            all_channels.append({"basis": basis, **row})
    write_basis_summary(output_root / "ar_basis_convergence.csv", summaries)
    with open(output_root / "ar_channel_convergence.csv", "w", newline="", encoding="utf-8") as fp:
        fieldnames = ["basis"] + [
            "atom",
            "channel",
            "delta_Ha",
            "osc",
            "alpha0_contribution",
            "alpha0_fraction",
            "single_channel_c6",
            "single_channel_c6_fraction",
            "cross_inclusive_c6",
            "cross_inclusive_c6_fraction",
            "is_residual",
        ]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_channels)
    return summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ar basis convergence for radial EFT spectral channels.")
    parser.add_argument("--basis", action="append", choices=BASIS_LIST, help="Basis to run. Repeatable.")
    parser.add_argument("--output-root", default="results/ar")
    parser.add_argument("--n-grid", type=int, default=2400)
    args = parser.parse_args(argv)

    summaries = run_basis_convergence(args.basis or BASIS_LIST, Path(args.output_root), n_grid=args.n_grid)
    print(
        "basis,alpha0,C6,sum_osc_discrete/N_core,sum_osc_total/N_core,"
        "3p_to_3d_delta,3p_to_3d_osc,3p_to_3d_alpha_fraction,residual_fraction"
    )
    for row in summaries:
        print(
            f"{row['basis']},"
            f"{row['alpha0']:.8f},"
            f"{row['C6']:.8f},"
            f"{row['sum_osc_discrete_over_N_core']:.8f},"
            f"{row['sum_osc_total_over_N_core']:.8f},"
            f"{row['3p_to_3d_delta']:.12f},"
            f"{row['3p_to_3d_osc']:.12f},"
            f"{row['3p_to_3d_alpha_fraction']:.8f},"
            f"{row['residual_fraction']:.8f}"
        )


if __name__ == "__main__":
    main()
