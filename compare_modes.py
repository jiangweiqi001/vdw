import argparse
import csv
import shutil
import subprocess
import sys

from compare_alpha_c6 import compare_tables


def format_mode_comparison_rows(mode, comparison_rows):
    rows = []
    for row in comparison_rows:
        rows.append(
            {
                "mode": mode,
                "atom": row["atom"],
                "alpha0_eft": row["alpha0_eft"],
                "alpha0_ref": row["alpha0_ref"],
                "err_alpha_pct": row["err_alpha_pct"],
                "C6_eft": row["C6_eft"],
                "C6_ref": row["C6_ref"],
                "err_C6_pct": row["err_C6_pct"],
            }
        )
    return rows


def _run(command):
    subprocess.run(command, check=True)


def _write_mode_comparison(path, rows):
    fieldnames = [
        "mode",
        "atom",
        "alpha0_eft",
        "alpha0_ref",
        "err_alpha_pct",
        "C6_eft",
        "C6_ref",
        "err_C6_pct",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_mode(mode, python_executable=sys.executable):
    _run([python_executable, "build_eft_channels.py", "--mode", mode, "--output", "atomic_channels.csv"])
    _run([python_executable, "run_alpha_table.py"])
    shutil.copyfile("alpha_c6_table.csv", f"alpha_c6_table_{mode}.csv")
    _run([python_executable, "run_c6_table.py"])
    shutil.copyfile("c6_table.csv", f"c6_table_{mode}.csv")
    return format_mode_comparison_rows(mode, compare_tables("alpha_c6_table.csv", "reference_alpha_c6.csv"))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run calibrated and spectral pipelines and compare both to references.")
    parser.add_argument("--output", default="mode_alpha_c6_errors.csv")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args(argv)

    all_rows = []
    for mode in ("calibrated", "spectral"):
        all_rows.extend(run_mode(mode, args.python))

    _write_mode_comparison(args.output, all_rows)
    print("mode,atom,alpha0_eft,alpha0_ref,err_alpha_pct,C6_eft,C6_ref,err_C6_pct")
    for row in all_rows:
        print(
            f"{row['mode']},"
            f"{row['atom']},"
            f"{row['alpha0_eft']:.8f},"
            f"{row['alpha0_ref']:.8f},"
            f"{row['err_alpha_pct']:.6f},"
            f"{row['C6_eft']:.8f},"
            f"{row['C6_ref']:.8f},"
            f"{row['err_C6_pct']:.6f}"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
