import argparse
import csv
from pathlib import Path

import numpy as np


HARTREE_TO_MEV = 27211.386245988


def radial_derivatives(r, model="bare", epsilon=1.0, kappa=0.0):
    r = float(r)
    epsilon = float(epsilon)
    kappa = float(kappa)
    if r <= 0.0:
        raise ValueError("r must be positive.")
    if model == "bare":
        phi_prime = -1.0 / r**2
        phi_second = 2.0 / r**3
    elif model == "dielectric":
        phi_prime = -1.0 / (epsilon * r**2)
        phi_second = 2.0 / (epsilon * r**3)
    elif model == "yukawa":
        expkr = np.exp(-kappa * r)
        phi_prime = -expkr * (kappa / r + 1.0 / r**2) / epsilon
        phi_second = expkr * (kappa**2 / r + 2.0 * kappa / r**2 + 2.0 / r**3) / epsilon
    else:
        raise ValueError(f"Unknown screening model: {model}")
    return phi_prime, phi_second


def dipole_tensor_trace_square(r, model="bare", epsilon=1.0, kappa=0.0):
    phi_prime, phi_second = radial_derivatives(r, model, epsilon, kappa)
    a = phi_prime / float(r)
    b = phi_second - phi_prime / float(r)
    # Hessian eigenvalues for a radial potential are a,a,a+b.
    return float(2.0 * a**2 + (a + b) ** 2)


def screened_pair_energy(c6, r_bohr, model="bare", epsilon=1.0, kappa=0.0):
    """Pairwise isotropic energy using a model screened dipole tensor.

    For the bare Coulomb tensor, Tr[T T] = 6/R^6 and this reduces to -C6/R^6.
    """
    trace_square = dipole_tensor_trace_square(r_bohr, model, epsilon, kappa)
    return -float(c6) * trace_square / 6.0


def build_tail_rows(c6, r_values, model="bare", epsilon=1.0, kappa=0.0):
    rows = []
    for r in r_values:
        e = screened_pair_energy(c6, r, model, epsilon, kappa)
        rows.append(
            {
                "R_bohr": float(r),
                "model": model,
                "epsilon": float(epsilon),
                "kappa": float(kappa),
                "C6": float(c6),
                "E_Ha": e,
                "E_meV": e * HARTREE_TO_MEV,
            }
        )
    return rows


def write_rows(path, rows):
    fieldnames = ["R_bohr", "model", "epsilon", "kappa", "C6", "E_Ha", "E_meV"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Model-screened pairwise -C6 tail using a screened W_v dipole tensor.")
    parser.add_argument("--c6", type=float, required=True)
    parser.add_argument("--r", type=float, action="append", default=[8, 10, 12, 15, 20, 30, 40])
    parser.add_argument("--model", choices=["bare", "dielectric", "yukawa"], default="bare")
    parser.add_argument("--epsilon", type=float, default=1.0)
    parser.add_argument("--kappa", type=float, default=0.0)
    parser.add_argument("--output", default="results/screened_pairwise_tail.csv")
    args = parser.parse_args(argv)

    rows = build_tail_rows(args.c6, args.r, args.model, args.epsilon, args.kappa)
    write_rows(args.output, rows)
    print("R_bohr,model,epsilon,kappa,C6,E_Ha,E_meV")
    for row in rows:
        print(
            f"{row['R_bohr']},"
            f"{row['model']},"
            f"{row['epsilon']},"
            f"{row['kappa']},"
            f"{row['C6']},"
            f"{row['E_Ha']:.12e},"
            f"{row['E_meV']:.12e}"
        )


if __name__ == "__main__":
    main()
