import argparse
import csv

import numpy as np


def _load_radial_orbitals(path):
    orbitals = {}
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            key = (row["atom"].strip(), row["orbital"].strip())
            orbitals.setdefault(key, {"r": [], "u": []})
            orbitals[key]["r"].append(float(row["r_bohr"]))
            orbitals[key]["u"].append(float(row["u"]))
    return orbitals


def _orbital_report(atom, orbital, values, norm_tol):
    r = np.asarray(values["r"], dtype=float)
    u = np.asarray(values["u"], dtype=float)
    monotonic = bool(np.all(np.diff(r) > 0.0))
    if monotonic:
        norm = float(np.trapezoid(u * u, r))
        r_min = float(r[0])
        r_max = float(r[-1])
    else:
        order = np.argsort(r)
        r_sorted = r[order]
        u_sorted = u[order]
        norm = float(np.trapezoid(u_sorted * u_sorted, r_sorted))
        r_min = float(r_sorted[0])
        r_max = float(r_sorted[-1])

    norm_error = abs(norm - 1.0)
    return {
        "atom": atom,
        "orbital": orbital,
        "n_grid": int(len(r)),
        "r_min": r_min,
        "r_max": r_max,
        "norm": norm,
        "norm_error": norm_error,
        "norm_ok": bool(norm_error <= norm_tol),
        "monotonic_grid": monotonic,
    }


def _coverage_warnings(reports, coverage_tol):
    warnings = []
    by_atom = {}
    for row in reports:
        by_atom.setdefault(row["atom"], []).append(row)

    for atom, rows in sorted(by_atom.items()):
        min_values = [row["r_min"] for row in rows]
        max_values = [row["r_max"] for row in rows]
        if max(min_values) - min(min_values) > coverage_tol or max(max_values) - min(max_values) > coverage_tol:
            warnings.append(
                {
                    "atom": atom,
                    "message": "radial ranges differ across orbitals",
                    "r_min_min": min(min_values),
                    "r_min_max": max(min_values),
                    "r_max_min": min(max_values),
                    "r_max_max": max(max_values),
                }
            )
    return warnings


def check_radial_orbitals(path="radial_orbitals.csv", norm_tol=1e-4, coverage_tol=1e-8):
    radial = _load_radial_orbitals(path)
    reports = [
        _orbital_report(atom, orbital, values, norm_tol)
        for (atom, orbital), values in sorted(radial.items())
    ]
    return {
        "orbitals": reports,
        "coverage_warnings": _coverage_warnings(reports, coverage_tol),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate reduced radial orbital grids u(r).")
    parser.add_argument("--input", default="radial_orbitals.csv")
    parser.add_argument("--norm-tol", type=float, default=1e-4)
    parser.add_argument("--coverage-tol", type=float, default=1e-8)
    args = parser.parse_args(argv)

    report = check_radial_orbitals(args.input, args.norm_tol, args.coverage_tol)
    print("atom,orbital,n_grid,r_min,r_max,norm,norm_error,norm_ok,monotonic_grid")
    for row in report["orbitals"]:
        print(
            f"{row['atom']},"
            f"{row['orbital']},"
            f"{row['n_grid']},"
            f"{row['r_min']:.12f},"
            f"{row['r_max']:.12f},"
            f"{row['norm']:.12f},"
            f"{row['norm_error']:.12e},"
            f"{str(row['norm_ok']).lower()},"
            f"{str(row['monotonic_grid']).lower()}"
        )

    if report["coverage_warnings"]:
        print("\ncoverage_warning,atom,message,r_min_min,r_min_max,r_max_min,r_max_max")
        for warning in report["coverage_warnings"]:
            print(
                f"coverage_warning,"
                f"{warning['atom']},"
                f"{warning['message']},"
                f"{warning['r_min_min']:.12f},"
                f"{warning['r_min_max']:.12f},"
                f"{warning['r_max_min']:.12f},"
                f"{warning['r_max_max']:.12f}"
            )


if __name__ == "__main__":
    main()
