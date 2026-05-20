import argparse
import csv
from pathlib import Path

import numpy as np

from convert_atomic_solver_output import convert_rows


DEFAULT_TARGETS = {
    "Mg": ["2s"],
    "Ca": ["3s"],
}


def cumulative_trapezoid(y, x):
    out = np.zeros_like(y, dtype=float)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1]))
    return out


def orbital_hartree_potential(r, u):
    r = np.asarray(r, dtype=float)
    u2 = np.asarray(u, dtype=float) ** 2
    inner = cumulative_trapezoid(u2, r) / r
    tail_integrand = np.divide(u2, r, out=np.zeros_like(u2), where=r != 0.0)
    outer_total = np.trapezoid(tail_integrand, r)
    outer = outer_total - cumulative_trapezoid(tail_integrand, r)
    return inner + outer


def self_coulomb(r, u):
    vh = orbital_hartree_potential(r, u)
    return float(np.trapezoid(np.asarray(u) ** 2 * vh, r))


def fK_values(r, u, vh, jc, k_grid):
    r = np.asarray(r, dtype=float)
    u = np.asarray(u, dtype=float)
    kernel = u * (vh - jc)
    values = []
    f0 = float(np.sqrt(4.0 * np.pi) * np.trapezoid(kernel * r, r))
    for k in k_grid:
        k = float(k)
        if abs(k) < 1e-12:
            values.append(f0)
        else:
            values.append(float(np.sqrt(4.0 * np.pi) / k * np.trapezoid(kernel * np.sin(k * r), r)))
    return f0, values


def _radial_from_shell(atom, basis, shell_label, r_min=1e-5, r_max=40.0, n_grid=2400):
    from pyscf.scf import atom_hf
    from pyscf_export_ar_radials import grouped_shells, make_atom_molecule

    mol = make_atom_molecule(atom, basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} atom HF did not converge for {basis}.")

    shell = next(shell for shell in grouped_shells(mol, mf.mo_energy, mf.mo_occ, mf.mo_coeff) if shell["label"] == shell_label)
    r_grid = np.geomspace(r_min, r_max, n_grid)
    coords = np.column_stack([np.zeros_like(r_grid), np.zeros_like(r_grid), r_grid])
    ao_values = mol.eval_gto("GTOval_sph", coords)
    mo_values = ao_values @ mf.mo_coeff
    l_value = shell["l"]
    y_l0_at_z = np.sqrt((2 * l_value + 1) / (4.0 * np.pi))
    column = max(shell["columns"], key=lambda idx: np.max(np.abs(mo_values[:, idx])))
    radial_R = mo_values[:, column] / y_l0_at_z
    converted = convert_rows(
        [
            {"atom": atom, "orbital": shell_label, "r_bohr": f"{r:.12g}", "R": f"{R:.12g}"}
            for r, R in zip(r_grid, radial_R)
        ],
        input_kind="R",
        normalize=True,
    )
    r = np.asarray([row["r_bohr"] for row in converted], dtype=float)
    u = np.asarray([row["u"] for row in converted], dtype=float)
    return r, u, shell


def compute_wilson_row(atom, shell_label, basis, k_grid):
    r, u, shell = _radial_from_shell(atom, basis, shell_label)
    vh = orbital_hartree_potential(r, u)
    jc = self_coulomb(r, u)
    f0, fk = fK_values(r, u, vh, jc, k_grid)
    delta = -float(shell["energy_Ha"])
    return {
        "atom": atom,
        "shell": shell_label,
        "basis": basis,
        "Delta_E_Ha": delta,
        "occupation": float(shell["occupation"]),
        "Jc_Ha": jc,
        "f0": f0,
        "f0_over_delta": f0 / delta if delta else 0.0,
        "K_grid": ";".join(f"{k:.8g}" for k in k_grid),
        "fK_values": ";".join(f"{v:.12g}" for v in fk),
        "source": "PRL_SCALAR_WILSON_KOOPMANS_HF",
    }


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "shell",
        "basis",
        "Delta_E_Ha",
        "occupation",
        "Jc_Ha",
        "f0",
        "f0_over_delta",
        "K_grid",
        "fK_values",
        "source",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute PRL scalar core Wilson coefficients for s shells.")
    parser.add_argument("--atom", action="append")
    parser.add_argument("--shell", action="append")
    parser.add_argument("--basis", default="cc-pVQZ")
    parser.add_argument("--k-max", type=float, default=8.0)
    parser.add_argument("--n-k", type=int, default=81)
    parser.add_argument("--output", default="results/core_wilson_coefficients.csv")
    args = parser.parse_args(argv)

    k_grid = np.linspace(0.0, args.k_max, args.n_k)
    rows = []
    if args.atom and args.shell:
        for atom, shell in zip(args.atom, args.shell):
            rows.append(compute_wilson_row(atom, shell, args.basis, k_grid))
    else:
        for atom, shells in DEFAULT_TARGETS.items():
            for shell in shells:
                rows.append(compute_wilson_row(atom, shell, args.basis, k_grid))
    write_rows(args.output, rows)
    print("atom,shell,basis,Delta_E_Ha,Jc_Ha,f0,f0_over_delta,source")
    for row in rows:
        print(
            f"{row['atom']},{row['shell']},{row['basis']},"
            f"{row['Delta_E_Ha']:.8f},{row['Jc_Ha']:.8f},"
            f"{row['f0']:.8f},{row['f0_over_delta']:.8f},{row['source']}"
        )


if __name__ == "__main__":
    main()
