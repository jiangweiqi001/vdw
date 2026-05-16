import argparse
import csv
from pathlib import Path

from pyscf_export_ar_mo_oscillators import export_ar_mo_oscillators, write_mo_channels
from pyscf_export_ar_tdhf_decomposed import export_to_dir
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from split_mo_channels_by_shell import decompose_mo_channels


REFERENCE_C6 = 64.3


def percent_error(c6):
    return 100.0 * (float(c6) - REFERENCE_C6) / REFERENCE_C6


def analyze_channel_file(path, prefix):
    alpha_rows = build_alpha_rows(path)
    write_alpha_rows(f"{prefix}_alpha_c6_table.csv", alpha_rows)
    write_c6_rows(f"{prefix}_c6_table.csv", build_c6_rows(path))
    return next(row for row in alpha_rows if row["atom"] == "Ar")


def result_row(method, component, alpha_row, note=""):
    c6 = float(alpha_row["C6_self_au"])
    return {
        "method": method,
        "component": component,
        "alpha0": float(alpha_row["alpha0_au"]),
        "C6": c6,
        "C6_error_pct": percent_error(c6),
        "n_channels": int(alpha_row["n_channels"]),
        "note": note,
    }


def write_rows(path, rows):
    fieldnames = ["method", "component", "alpha0", "C6", "C6_error_pct", "n_channels", "note"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_decomposition(basis="aug-cc-pVQZ", nstates=200, output_root="results/ar"):
    output_root = Path(output_root)
    mo_dir = output_root / "core_valence_mo"
    tdhf_dir = output_root / "core_valence_tdhf"
    mo_dir.mkdir(parents=True, exist_ok=True)
    tdhf_dir.mkdir(parents=True, exist_ok=True)

    mo_channels, _mo_summary = export_ar_mo_oscillators(basis)
    mo_all_path = mo_dir / "ar_mo_all_channels.csv"
    write_mo_channels(mo_all_path, mo_channels)
    mo_summaries = decompose_mo_channels(mo_all_path, mo_dir, basis)

    export_to_dir(tdhf_dir, basis, nstates)
    tdhf_files = {
        "all": tdhf_dir / "ar_tdhf_all_channels.csv",
        "valence": tdhf_dir / "ar_tdhf_valence_projected_channels.csv",
        "core": tdhf_dir / "ar_tdhf_core_projected_channels.csv",
        "cross": tdhf_dir / "ar_tdhf_cross_projected_channels.csv",
    }
    tdhf_summaries = {
        component: analyze_channel_file(path, tdhf_dir / component)
        for component, path in tdhf_files.items()
    }

    rows = []
    for component, summary in mo_summaries.items():
        rows.append(result_row("MO", component, summary, "occupied-shell split"))
    mo_all_c6 = float(mo_summaries["all"]["C6_self_au"])
    mo_val_c6 = float(mo_summaries["valence"]["C6_self_au"])
    rows.append(
        {
            "method": "MO",
            "component": "Delta_C6_core",
            "alpha0": "",
            "C6": mo_all_c6 - mo_val_c6,
            "C6_error_pct": "",
            "n_channels": "",
            "note": f"relative_core_contribution={(mo_all_c6 - mo_val_c6) / mo_all_c6:.8f}",
        }
    )

    for component, summary in tdhf_summaries.items():
        rows.append(result_row("TDHF", component, summary, "transition-dipole projection"))
    tdhf_all_c6 = float(tdhf_summaries["all"]["C6_self_au"])
    tdhf_val_c6 = float(tdhf_summaries["valence"]["C6_self_au"])
    tdhf_core_c6 = float(tdhf_summaries["core"]["C6_self_au"])
    tdhf_val_core_rows = []
    with open(tdhf_dir / "ar_tdhf_valence_plus_core_channels.csv", "w", newline="", encoding="utf-8") as fp:
        fieldnames = ["atom", "channel", "delta_Ha", "osc", "is_core", "source", "component"]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for component in ("valence", "core"):
            with open(tdhf_files[component], newline="", encoding="utf-8") as in_fp:
                for row in csv.DictReader(in_fp):
                    writer.writerow(row)
                    tdhf_val_core_rows.append(row)
    val_core_summary = analyze_channel_file(tdhf_dir / "ar_tdhf_valence_plus_core_channels.csv", tdhf_dir / "valence_plus_core")
    rows.append(
        {
            "method": "TDHF",
            "component": "cross_alpha0",
            "alpha0": float(tdhf_summaries["all"]["alpha0_au"]) - float(tdhf_summaries["valence"]["alpha0_au"]) - float(tdhf_summaries["core"]["alpha0_au"]),
            "C6": "",
            "C6_error_pct": "",
            "n_channels": "",
            "note": "alpha0_all - alpha0_valence - alpha0_core",
        }
    )
    rows.append(
        {
            "method": "TDHF",
            "component": "cross_C6_effect",
            "alpha0": "",
            "C6": tdhf_all_c6 - float(val_core_summary["C6_self_au"]),
            "C6_error_pct": "",
            "n_channels": "",
            "note": "C6_all - C6[valence+core_without_cross]",
        }
    )
    rows.append(
        {
            "method": "TDHF",
            "component": "Delta_C6_core",
            "alpha0": "",
            "C6": tdhf_all_c6 - tdhf_val_c6,
            "C6_error_pct": "",
            "n_channels": "",
            "note": f"relative_core_contribution={(tdhf_all_c6 - tdhf_val_c6) / tdhf_all_c6:.8f}; core_only_C6={tdhf_core_c6:.8f}",
        }
    )

    write_rows(output_root / "ar_core_valence_decomposition.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Ar MO and TDHF core/valence decomposition.")
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--nstates", type=int, default=200)
    parser.add_argument("--output-root", default="results/ar")
    args = parser.parse_args(argv)

    rows = run_decomposition(args.basis, args.nstates, args.output_root)
    print("method,component,alpha0,C6,C6_error_pct,n_channels,note")
    for row in rows:
        print(
            f"{row['method']},"
            f"{row['component']},"
            f"{row['alpha0']},"
            f"{row['C6']},"
            f"{row['C6_error_pct']},"
            f"{row['n_channels']},"
            f"{row['note']}"
        )


if __name__ == "__main__":
    main()
