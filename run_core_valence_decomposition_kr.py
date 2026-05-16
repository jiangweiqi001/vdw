import argparse
import csv
from pathlib import Path

from pyscf_export_ar_tdhf_decomposed import export_decomposed_tdhf, write_component_channels
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows


KR_VALENCE_SHELLS = {"4s", "4p"}
KR_CORE_SHELLS = {"1s", "2s", "2p", "3s", "3p", "3d"}
REFERENCE_C6 = 129.6


def percent_error(c6):
    return 100.0 * (float(c6) - REFERENCE_C6) / REFERENCE_C6


def analyze_channel_file(path, prefix):
    alpha_rows = build_alpha_rows(path)
    write_alpha_rows(f"{prefix}_alpha_c6_table.csv", alpha_rows)
    write_c6_rows(f"{prefix}_c6_table.csv", build_c6_rows(path))
    return next(row for row in alpha_rows if row["atom"] == "Kr")


def result_row(component, alpha_row, note="transition-dipole projection"):
    c6 = float(alpha_row["C6_self_au"])
    return {
        "method": "TDHF",
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


def run_decomposition(basis="aug-cc-pVQZ", nstates=200, output_root="results/kr"):
    output_root = Path(output_root)
    tdhf_dir = output_root / "core_valence_tdhf"
    tdhf_dir.mkdir(parents=True, exist_ok=True)

    component_rows, _summary = export_decomposed_tdhf(
        atom="Kr",
        basis=basis,
        nstates=nstates,
        valence_shells=KR_VALENCE_SHELLS,
        core_shells=KR_CORE_SHELLS,
    )
    files = {
        "all": tdhf_dir / "kr_tdhf_all_channels.csv",
        "valence": tdhf_dir / "kr_tdhf_valence_projected_channels.csv",
        "core": tdhf_dir / "kr_tdhf_core_projected_channels.csv",
        "cross": tdhf_dir / "kr_tdhf_cross_projected_channels.csv",
    }
    for component, rows in component_rows.items():
        write_component_channels(files[component], rows)

    summaries = {
        component: analyze_channel_file(path, tdhf_dir / component)
        for component, path in files.items()
    }

    with open(tdhf_dir / "kr_tdhf_valence_plus_core_channels.csv", "w", newline="", encoding="utf-8") as fp:
        fieldnames = ["atom", "channel", "delta_Ha", "osc", "is_core", "source", "component"]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for component in ("valence", "core"):
            with open(files[component], newline="", encoding="utf-8") as in_fp:
                for row in csv.DictReader(in_fp):
                    writer.writerow(row)
    val_core_summary = analyze_channel_file(
        tdhf_dir / "kr_tdhf_valence_plus_core_channels.csv",
        tdhf_dir / "valence_plus_core",
    )

    rows = [result_row(component, summary) for component, summary in summaries.items()]
    all_alpha = float(summaries["all"]["alpha0_au"])
    val_alpha = float(summaries["valence"]["alpha0_au"])
    core_alpha = float(summaries["core"]["alpha0_au"])
    all_c6 = float(summaries["all"]["C6_self_au"])
    val_c6 = float(summaries["valence"]["C6_self_au"])
    core_c6 = float(summaries["core"]["C6_self_au"])
    rows.extend(
        [
            {
                "method": "TDHF",
                "component": "cross_alpha0",
                "alpha0": all_alpha - val_alpha - core_alpha,
                "C6": "",
                "C6_error_pct": "",
                "n_channels": "",
                "note": "alpha0_all - alpha0_valence - alpha0_core",
            },
            {
                "method": "TDHF",
                "component": "cross_C6_effect",
                "alpha0": "",
                "C6": all_c6 - float(val_core_summary["C6_self_au"]),
                "C6_error_pct": "",
                "n_channels": "",
                "note": "C6_all - C6[valence+core_without_cross]",
            },
            {
                "method": "TDHF",
                "component": "Delta_C6_core",
                "alpha0": "",
                "C6": all_c6 - val_c6,
                "C6_error_pct": "",
                "n_channels": "",
                "note": f"relative_core_contribution={(all_c6 - val_c6) / all_c6:.8f}; core_only_C6={core_c6:.8f}",
            },
        ]
    )

    write_rows(output_root / "kr_core_valence_decomposition.csv", rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run Kr TDHF core/valence decomposition.")
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--nstates", type=int, default=200)
    parser.add_argument("--output-root", default="results/kr")
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
