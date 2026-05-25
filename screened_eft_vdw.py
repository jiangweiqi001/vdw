import argparse
import csv
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from screened_pairwise_vdw import radial_derivatives


def load_channel_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def group_channels_by_atom(rows):
    grouped = {}
    for row in rows:
        atom = row["atom"].strip()
        grouped.setdefault(atom, {"delta": [], "osc": []})
        grouped[atom]["delta"].append(float(row["delta_Ha"]))
        grouped[atom]["osc"].append(float(row["osc"]))
    return grouped


def alpha_iw(delta, osc, xi):
    delta = np.asarray(delta, dtype=float)
    osc = np.asarray(osc, dtype=float)
    return float(np.sum(osc / (delta**2 + xi**2)))


def quadrature_grid(n_quad=120):
    x, w = leggauss(n_quad)
    t = 0.25 * np.pi * (x + 1.0)
    xi = np.tan(t)
    jac = 0.25 * np.pi / np.cos(t) ** 2
    return xi, w * jac


def screened_dipole_tensor(r_vec, screening=None):
    screening = screening or {"model": "bare"}
    r_vec = np.asarray(r_vec, dtype=float)
    r = float(np.linalg.norm(r_vec))
    if r <= 0.0:
        raise ValueError("Dipole tensor is undefined for zero separation.")
    n = r_vec / r
    model = screening.get("model", "bare")
    epsilon = float(screening.get("epsilon", 1.0))
    kappa = float(screening.get("kappa", 0.0))
    phi_prime, phi_second = radial_derivatives(r, model=model, epsilon=epsilon, kappa=kappa)
    a = phi_prime / r
    b = phi_second - phi_prime / r
    return a * np.eye(3) + b * np.outer(n, n)


def build_coupling_matrix(alphas, positions, screening=None):
    n_atoms = len(alphas)
    mat = np.zeros((3 * n_atoms, 3 * n_atoms), dtype=float)
    positions = np.asarray(positions, dtype=float)
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i == j:
                continue
            tensor = screened_dipole_tensor(positions[i] - positions[j], screening=screening)
            block = np.sqrt(alphas[i] * alphas[j]) * tensor
            mat[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = block
    return mat


def logdet_integrand(coupling, expansion_order=None):
    if expansion_order == 2:
        return -0.5 * float(np.trace(coupling @ coupling))
    sign, logdet = np.linalg.slogdet(np.eye(coupling.shape[0]) - coupling)
    if sign <= 0:
        # Large model couplings can leave the perturbative domain. Fall back to
        # real part of slogdet-compatible eigenvalue expression for diagnostics.
        eig = np.linalg.eigvals(coupling)
        logdet = float(np.sum(np.log(np.abs(1.0 - eig))))
    return logdet + float(np.trace(coupling))


def logdet_vdw_energy(channel_rows, positions_bohr, atom_order=None, screening=None, n_quad=120, expansion_order=None):
    grouped = group_channels_by_atom(channel_rows)
    if atom_order is None:
        atom_order = [row["atom"].strip() for row in channel_rows]
        # Preserve first occurrence of each atom label.
        atom_order = list(dict.fromkeys(atom_order))
        if len(atom_order) != len(positions_bohr):
            if len(grouped) == 1:
                atom_order = [next(iter(grouped))] * len(positions_bohr)
            else:
                raise ValueError("--atoms is required when multiple atom labels are present.")
    if len(atom_order) != len(positions_bohr):
        raise ValueError("atom_order and positions must have the same length.")

    xi_grid, weights = quadrature_grid(n_quad)
    total = 0.0
    for xi, weight in zip(xi_grid, weights):
        alphas = [
            alpha_iw(grouped[atom]["delta"], grouped[atom]["osc"], xi)
            for atom in atom_order
        ]
        coupling = build_coupling_matrix(alphas, positions_bohr, screening=screening)
        total += weight * logdet_integrand(coupling, expansion_order=expansion_order)
    return total / (2.0 * np.pi)


def read_positions(path):
    positions = []
    atoms = []
    with open(path, newline="", encoding="utf-8") as fp:
        for row in csv.DictReader(fp):
            atoms.append(row["atom"].strip())
            positions.append([float(row["x_bohr"]), float(row["y_bohr"]), float(row["z_bohr"])])
    return atoms, positions


def main(argv=None):
    parser = argparse.ArgumentParser(description="Finite-system screened EFT-vdW logdet prototype.")
    parser.add_argument("--channels", required=True)
    parser.add_argument("--positions", required=True, help="CSV: atom,x_bohr,y_bohr,z_bohr")
    parser.add_argument("--model", choices=["bare", "dielectric", "yukawa"], default="bare")
    parser.add_argument("--epsilon", type=float, default=1.0)
    parser.add_argument("--kappa", type=float, default=0.0)
    parser.add_argument("--n-quad", type=int, default=120)
    parser.add_argument("--second-order", action="store_true")
    parser.add_argument("--output", default="results/screened_eft_vdw_energy.csv")
    args = parser.parse_args(argv)

    atoms, positions = read_positions(args.positions)
    screening = {"model": args.model, "epsilon": args.epsilon, "kappa": args.kappa}
    energy = logdet_vdw_energy(
        load_channel_rows(args.channels),
        positions,
        atom_order=atoms,
        screening=screening,
        n_quad=args.n_quad,
        expansion_order=2 if args.second_order else None,
    )
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=["model", "epsilon", "kappa", "energy_Ha"])
        writer.writeheader()
        writer.writerow({"model": args.model, "epsilon": args.epsilon, "kappa": args.kappa, "energy_Ha": energy})
    print("model,epsilon,kappa,energy_Ha")
    print(f"{args.model},{args.epsilon},{args.kappa},{energy:.12e}")


if __name__ == "__main__":
    main()
