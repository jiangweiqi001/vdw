import argparse
import csv
from pathlib import Path


HARTREE_TO_MEV = 27211.386245988
R_VALUES = [8, 10, 12, 15, 20, 30, 40]


def energy_hartree(c6, r_bohr):
    return -float(c6) / float(r_bohr) ** 6


def percent_error(value, reference):
    return 100.0 * (value - reference) / reference


def tail_row(r_bohr, c6_ref, c6_tdhf, c6_mo, c6_calibrated):
    e_ref = energy_hartree(c6_ref, r_bohr)
    e_tdhf = energy_hartree(c6_tdhf, r_bohr)
    e_mo = energy_hartree(c6_mo, r_bohr)
    e_calibrated = energy_hartree(c6_calibrated, r_bohr)
    return {
        "R_bohr": float(r_bohr),
        "E_ref_Ha": e_ref,
        "E_ref_meV": e_ref * HARTREE_TO_MEV,
        "E_tdhf_Ha": e_tdhf,
        "E_tdhf_meV": e_tdhf * HARTREE_TO_MEV,
        "E_mo_Ha": e_mo,
        "E_mo_meV": e_mo * HARTREE_TO_MEV,
        "E_calibrated_Ha": e_calibrated,
        "E_calibrated_meV": e_calibrated * HARTREE_TO_MEV,
        "err_tdhf_pct": percent_error(c6_tdhf, c6_ref),
        "err_mo_pct": percent_error(c6_mo, c6_ref),
    }


def build_tail_rows(c6_ref=64.3, c6_tdhf=60.73027908, c6_mo=76.13748508, c6_calibrated=64.3):
    return [tail_row(r, c6_ref, c6_tdhf, c6_mo, c6_calibrated) for r in R_VALUES]


def write_tail_rows(path, rows):
    fieldnames = [
        "R_bohr",
        "E_ref_Ha",
        "E_ref_meV",
        "E_tdhf_Ha",
        "E_tdhf_meV",
        "E_mo_Ha",
        "E_mo_meV",
        "E_calibrated_Ha",
        "E_calibrated_meV",
        "err_tdhf_pct",
        "err_mo_pct",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare Ar2 long-range -C6/R^6 tails.")
    parser.add_argument("--output", default="results/ar/ar2_tail_comparison.csv")
    parser.add_argument("--c6-ref", type=float, default=64.3)
    parser.add_argument("--c6-tdhf", type=float, default=60.73027908)
    parser.add_argument("--c6-mo", type=float, default=76.13748508)
    parser.add_argument("--c6-calibrated", type=float, default=64.3)
    args = parser.parse_args(argv)

    rows = build_tail_rows(args.c6_ref, args.c6_tdhf, args.c6_mo, args.c6_calibrated)
    write_tail_rows(args.output, rows)
    print("R_bohr,E_ref_meV,E_tdhf_meV,E_mo_meV,E_calibrated_meV,err_tdhf_pct,err_mo_pct")
    for row in rows:
        print(
            f"{row['R_bohr']:.0f},"
            f"{row['E_ref_meV']:.12e},"
            f"{row['E_tdhf_meV']:.12e},"
            f"{row['E_mo_meV']:.12e},"
            f"{row['E_calibrated_meV']:.12e},"
            f"{row['err_tdhf_pct']:.6f},"
            f"{row['err_mo_pct']:.6f}"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
