import argparse
import csv
from pathlib import Path

from compute_core_sternheimer import compute_core_sternheimer_channels, write_channels
from run_alpha_table import build_alpha_rows
from semicore_c6_targets import SEMICORE_C6_TARGETS


DEFAULT_BASES = ["def2-TZVP", "def2-TZVPP", "ano-rcc"]


def core_electron_count(atom):
    return sum(
        2 if shell.endswith("s") else 6 if shell.endswith("p") else 10 if shell.endswith("d") else 0
        for shell in SEMICORE_C6_TARGETS[atom].correction_shells
    )


def diagnostic_row(atom, basis, output_root):
    target = SEMICORE_C6_TARGETS[atom]
    output_root = Path(output_root)
    channels_path = output_root / "channels" / f"{atom.lower()}_{basis.replace('-', '_')}_core_channels.csv"
    rows = compute_core_sternheimer_channels(atom, basis, target.correction_shells)
    write_channels(channels_path, rows)
    alpha_row = next(row for row in build_alpha_rows(channels_path) if row["atom"] == atom)
    sum_osc = sum(float(row["osc"]) for row in rows)
    expected = core_electron_count(atom)
    return {
        "atom": atom,
        "basis": basis,
        "correction_shells": ";".join(sorted(target.correction_shells)),
        "n_channels": len(rows),
        "sum_osc": sum_osc,
        "expected_core_electrons": expected,
        "sum_osc_over_expected": sum_osc / expected if expected else "",
        "alpha0_core": float(alpha_row["alpha0_au"]),
        "C6_self_core": float(alpha_row["C6_self_au"]),
        "channels_path": str(channels_path),
    }


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "basis",
        "correction_shells",
        "n_channels",
        "sum_osc",
        "expected_core_electrons",
        "sum_osc_over_expected",
        "alpha0_core",
        "C6_self_core",
        "channels_path",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Diagnose finite-basis core Sternheimer oscillator closure.")
    parser.add_argument("--atom", action="append", choices=sorted(SEMICORE_C6_TARGETS), required=True)
    parser.add_argument("--basis", action="append")
    parser.add_argument("--output-root", default="results/core_sternheimer_basis_diagnostic")
    args = parser.parse_args(argv)

    bases = args.basis or DEFAULT_BASES
    rows = []
    for atom in args.atom:
        for basis in bases:
            try:
                rows.append(diagnostic_row(atom, basis, args.output_root))
            except Exception as exc:
                rows.append(
                    {
                        "atom": atom,
                        "basis": basis,
                        "correction_shells": ";".join(sorted(SEMICORE_C6_TARGETS[atom].correction_shells)),
                        "n_channels": "",
                        "sum_osc": "",
                        "expected_core_electrons": core_electron_count(atom),
                        "sum_osc_over_expected": "",
                        "alpha0_core": "",
                        "C6_self_core": "",
                        "channels_path": f"{type(exc).__name__}: {exc}",
                    }
                )
    output = Path(args.output_root) / "summary.csv"
    write_rows(output, rows)
    print("atom,basis,correction_shells,n_channels,sum_osc,expected_core_electrons,sum_osc_over_expected,alpha0_core,C6_self_core")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['basis']},"
            f"{row['correction_shells']},"
            f"{row['n_channels']},"
            f"{row['sum_osc']},"
            f"{row['expected_core_electrons']},"
            f"{row['sum_osc_over_expected']},"
            f"{row['alpha0_core']},"
            f"{row['C6_self_core']}"
        )
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()
