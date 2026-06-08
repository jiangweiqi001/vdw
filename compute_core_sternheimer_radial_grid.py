import argparse
import csv
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


L_LETTER = {0: "s", 1: "p", 2: "d", 3: "f"}
SOURCE = "PAULI_PROJECTED_TRK_CONSTRAINED_CPHF_RADIAL_GRID"


def angular_factor(l_i, l_a):
    if abs(l_a - l_i) != 1:
        return 0.0
    return max(l_i, l_a) / (2.0 * l_i + 1.0)


def normalize_u(r, u):
    norm = np.trapezoid(u * u, r)
    if norm <= 0.0:
        raise ValueError("Radial orbital norm must be positive.")
    return u / np.sqrt(norm)


def finite_difference_hamiltonian(r, l_value, potential):
    r = np.asarray(r, dtype=float)
    potential = np.asarray(potential, dtype=float)
    if len(r) < 5:
        raise ValueError("Radial grid must contain at least five points.")
    dr = float(r[1] - r[0])
    if not np.allclose(np.diff(r), dr, rtol=1e-6, atol=1e-12):
        raise ValueError("Radial Sternheimer grid must be uniform.")

    diag = np.full(len(r), 1.0 / dr**2)
    diag += l_value * (l_value + 1.0) / (2.0 * r**2)
    diag += potential
    off = np.full(len(r) - 1, -0.5 / dr**2)
    return np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)


def sternheimer_alpha_from_hamiltonian(r, hamiltonian, epsilon_occ, source_u, occupation, l_occ, l_resp, xi):
    """Solve radial imaginary-frequency Sternheimer response for one l channel.

    The solve uses the real form

        [(H - eps)^2 + xi^2] x = (H - eps) q

    with q(r) = r u_c(r).  This evaluates the same linear response as a
    complete spectral expansion of the supplied radial-grid Hamiltonian, but it
    does not require selecting bound virtual states.
    """
    factor = angular_factor(l_occ, l_resp)
    if factor == 0.0:
        return np.zeros_like(np.atleast_1d(xi), dtype=float)

    r = np.asarray(r, dtype=float)
    q = r * np.asarray(source_u, dtype=float)
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    dr = float(r[1] - r[0])
    a_matrix = np.asarray(hamiltonian, dtype=float) - float(epsilon_occ) * np.eye(len(r))
    aq = a_matrix @ q
    alpha = []
    for value in xi:
        lhs = a_matrix @ a_matrix + float(value) ** 2 * np.eye(len(r))
        response = np.linalg.solve(lhs, aq)
        contraction = float(np.dot(q, response) * dr)
        alpha.append((2.0 / 3.0) * float(occupation) * factor * contraction)
    return np.asarray(alpha, dtype=float)


def spectral_alpha_from_hamiltonian(r, hamiltonian, epsilon_occ, source_u, occupation, l_occ, l_resp, xi):
    """Diagnostic spectral expansion of the same radial-grid Hamiltonian."""
    factor = angular_factor(l_occ, l_resp)
    if factor == 0.0:
        return np.zeros_like(np.atleast_1d(xi), dtype=float), 0.0

    r = np.asarray(r, dtype=float)
    q = r * np.asarray(source_u, dtype=float)
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    dr = float(r[1] - r[0])
    eigvals, eigvecs = np.linalg.eigh(np.asarray(hamiltonian, dtype=float))
    deltas = eigvals - float(epsilon_occ)
    mask = deltas > 1e-10
    deltas = deltas[mask]
    vectors = eigvecs[:, mask] / np.sqrt(dr)
    dipoles = np.asarray([np.dot(q, vectors[:, idx]) * dr for idx in range(vectors.shape[1])])
    osc = (2.0 / 3.0) * float(occupation) * factor * deltas * dipoles**2
    alpha = np.sum(osc[None, :] / (deltas[None, :] ** 2 + xi[:, None] ** 2), axis=1)
    return alpha, float(np.sum(osc))


def spectral_components_from_hamiltonian(r, hamiltonian, epsilon_occ, source_u, occupation, l_occ, l_resp, xi):
    """Return signed and Pauli-positive spectral response components.

    This is a radial-box continuum discretization audit, not a Gaussian virtual
    MO sum.  The positive component represents Pauli-allowed absorption channels
    in the one-shell central-field picture.  Lower-energy occupied-state
    contributions are retained separately as the negative component so the
    signed commutator/TRK audit is not confused with the C6 absorption channel.
    """
    factor = angular_factor(l_occ, l_resp)
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    zeros = np.zeros_like(xi, dtype=float)
    if factor == 0.0:
        return {
            "signed_alpha": zeros,
            "positive_alpha": zeros,
            "negative_alpha": zeros,
            "signed_sum_osc": 0.0,
            "positive_sum_osc": 0.0,
            "negative_sum_osc": 0.0,
            "n_negative_delta": 0,
        }

    r = np.asarray(r, dtype=float)
    q = r * np.asarray(source_u, dtype=float)
    dr = float(r[1] - r[0])
    eigvals, eigvecs = np.linalg.eigh(np.asarray(hamiltonian, dtype=float))
    deltas = eigvals - float(epsilon_occ)
    vectors = eigvecs / np.sqrt(dr)
    dipoles = np.asarray([np.dot(q, vectors[:, idx]) * dr for idx in range(vectors.shape[1])])
    osc = (2.0 / 3.0) * float(occupation) * factor * deltas * dipoles**2
    positive = deltas > 1e-10
    negative = deltas < -1e-10

    def alpha_from_mask(mask):
        if not np.any(mask):
            return zeros
        return np.sum(
            osc[None, mask] / (deltas[None, mask] ** 2 + xi[:, None] ** 2),
            axis=1,
        )

    return {
        "signed_alpha": alpha_from_mask(np.ones_like(deltas, dtype=bool)),
        "positive_alpha": alpha_from_mask(positive),
        "negative_alpha": alpha_from_mask(negative),
        "signed_sum_osc": float(np.sum(osc)),
        "positive_sum_osc": float(np.sum(osc[positive])),
        "negative_sum_osc": float(np.sum(osc[negative])),
        "n_negative_delta": int(np.count_nonzero(negative)),
    }


def positive_spectral_channels_from_hamiltonian(
    atom,
    shell_label,
    r,
    hamiltonian,
    epsilon_occ,
    source_u,
    occupation,
    l_occ,
    l_resp,
    min_osc=1e-14,
):
    factor = angular_factor(l_occ, l_resp)
    if factor == 0.0:
        return []

    r = np.asarray(r, dtype=float)
    q = r * np.asarray(source_u, dtype=float)
    dr = float(r[1] - r[0])
    eigvals, eigvecs = np.linalg.eigh(np.asarray(hamiltonian, dtype=float))
    deltas = eigvals - float(epsilon_occ)
    vectors = eigvecs / np.sqrt(dr)
    dipoles = np.asarray([np.dot(q, vectors[:, idx]) * dr for idx in range(vectors.shape[1])])
    osc = (2.0 / 3.0) * float(occupation) * factor * deltas * dipoles**2
    rows = []
    for idx, (delta, strength, dipole) in enumerate(zip(deltas, osc, dipoles)):
        if delta <= 1e-10 or abs(strength) <= min_osc:
            continue
        rows.append(
            {
                "atom": atom,
                "channel": f"radial_cphf_{shell_label}_to_{L_LETTER.get(l_resp, l_resp)}_{idx:04d}",
                "delta_Ha": float(delta),
                "osc": float(strength),
                "is_core": "true",
                "source": SOURCE,
                "from_shell": shell_label,
                "response_l": int(l_resp),
                "response_basis_index": int(idx),
                "occ_energy_Ha": float(epsilon_occ),
                "response_energy_Ha": float(eigvals[idx]),
                "radial_dipole": float(dipole),
                "cphf_trk_scale": 1.0,
            }
        )
    return rows


def spectral_moments_from_hamiltonian(r, hamiltonian, epsilon_occ, source_u, occupation, l_occ, l_resp):
    """Return signed and absorption-only oscillator sums for a radial operator.

    The signed sum is the TRK/commutator audit of the central-field operator.
    The positive-only sum is the absorption strength that would enter a naive
    independent-particle channel model. For inner d shells these can differ
    because lower occupied p states carry negative oscillator strength in the
    one-shell audit and are not valid absorption channels.
    """
    factor = angular_factor(l_occ, l_resp)
    if factor == 0.0:
        return {
            "signed_sum_osc": 0.0,
            "positive_sum_osc": 0.0,
            "negative_sum_osc": 0.0,
            "n_negative_delta": 0,
        }

    components = spectral_components_from_hamiltonian(
        r,
        hamiltonian,
        epsilon_occ,
        source_u,
        occupation,
        l_occ,
        l_resp,
        [0.0],
    )
    return {key: components[key] for key in ["signed_sum_osc", "positive_sum_osc", "negative_sum_osc", "n_negative_delta"]}


def invert_local_potential(r, u, epsilon, l_value, node_tol=1e-5):
    r = np.asarray(r, dtype=float)
    u = normalize_u(r, np.asarray(u, dtype=float))
    dr = float(r[1] - r[0])
    second = (u[:-2] - 2.0 * u[1:-1] + u[2:]) / dr**2
    ri = r[1:-1]
    ui = u[1:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        potential = float(epsilon) + 0.5 * second / ui - l_value * (l_value + 1.0) / (2.0 * ri**2)
    valid = np.isfinite(potential) & (np.abs(ui) > node_tol * np.max(np.abs(u)))
    if np.count_nonzero(valid) < 5:
        raise ValueError("Too few stable points to invert local radial potential.")
    potential = np.interp(ri, ri[valid], potential[valid])
    return ri, u[1:-1], potential


def shell_radial_from_pyscf(atom, basis, shell_label, r_min=1e-4, r_max=60.0, n_grid=700):
    from pyscf.scf import atom_hf
    from pyscf_export_ar_radials import grouped_shells, make_atom_molecule

    mol = make_atom_molecule(atom, basis)
    mf = atom_hf.AtomSphAverageRHF(mol)
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"{atom} atom HF did not converge for {basis}.")

    shell = next(shell for shell in grouped_shells(mol, mf.mo_energy, mf.mo_occ, mf.mo_coeff) if shell["label"] == shell_label)
    r = np.linspace(r_min, r_max, int(n_grid))
    coords = np.column_stack([np.zeros_like(r), np.zeros_like(r), r])
    ao_values = mol.eval_gto("GTOval_sph", coords)
    mo_values = ao_values @ mf.mo_coeff
    l_value = int(shell["l"])
    y_l0_at_z = np.sqrt((2 * l_value + 1) / (4.0 * np.pi))
    column = max(shell["columns"], key=lambda idx: np.max(np.abs(mo_values[:, idx])))
    radial_r = mo_values[:, column] / y_l0_at_z
    u = normalize_u(r, r * radial_r)
    return {
        "atom": atom,
        "basis": basis,
        "shell": shell_label,
        "l": l_value,
        "occupation": float(shell["occupation"]),
        "epsilon": float(shell["energy_Ha"]),
        "r": r,
        "u": u,
    }


def quadrature_grid(n_quad=120):
    x, w = leggauss(int(n_quad))
    t = 0.25 * np.pi * (x + 1.0)
    xi = np.tan(t)
    weights = w * (0.25 * np.pi / np.cos(t) ** 2)
    return xi, weights


def self_c6_from_alpha_grid(alpha, weights):
    return float((3.0 / np.pi) * np.sum(weights * alpha * alpha))


def radial_grid_sternheimer_summary(atom, basis, shell_label, n_grid=700, r_max=60.0, n_quad=120):
    shell = shell_radial_from_pyscf(atom, basis, shell_label, r_max=r_max, n_grid=n_grid)
    r, u, potential = invert_local_potential(shell["r"], shell["u"], shell["epsilon"], shell["l"])
    xi, weights = quadrature_grid(n_quad)

    alpha_total = np.zeros_like(xi)
    positive_alpha_total = np.zeros_like(xi)
    signed_sum_osc = 0.0
    positive_sum_osc = 0.0
    negative_sum_osc = 0.0
    n_negative_delta = 0
    response_channels = []
    for l_resp in [shell["l"] - 1, shell["l"] + 1]:
        if l_resp < 0:
            continue
        hamiltonian = finite_difference_hamiltonian(r, l_resp, potential)
        alpha = sternheimer_alpha_from_hamiltonian(
            r,
            hamiltonian,
            shell["epsilon"],
            u,
            shell["occupation"],
            shell["l"],
            l_resp,
            xi,
        )
        components = spectral_components_from_hamiltonian(
            r,
            hamiltonian,
            shell["epsilon"],
            u,
            shell["occupation"],
            shell["l"],
            l_resp,
            xi,
        )
        alpha_total += alpha
        positive_alpha_total += components["positive_alpha"]
        signed_sum_osc += components["signed_sum_osc"]
        positive_sum_osc += components["positive_sum_osc"]
        negative_sum_osc += components["negative_sum_osc"]
        n_negative_delta += components["n_negative_delta"]
        response_channels.append(L_LETTER.get(l_resp, str(l_resp)))

    signed_ratio = signed_sum_osc / shell["occupation"] if shell["occupation"] else np.nan
    positive_ratio = positive_sum_osc / shell["occupation"] if shell["occupation"] else np.nan
    cphf_trk_scale = shell["occupation"] / positive_sum_osc if positive_sum_osc > 0.0 else np.nan
    cphf_alpha_total = positive_alpha_total * cphf_trk_scale
    cphf_positive_sum_osc = positive_sum_osc * cphf_trk_scale
    cphf_positive_ratio = cphf_positive_sum_osc / shell["occupation"] if shell["occupation"] else np.nan
    operator_trk_pass = abs(signed_ratio - 1.0) <= 0.05
    raw_absorption_trk_pass = abs(positive_ratio - 1.0) <= 0.05
    cphf_absorption_trk_pass = abs(cphf_positive_ratio - 1.0) <= 0.05
    return {
        "atom": atom,
        "basis": basis,
        "shell": shell_label,
        "response_l": ";".join(response_channels),
        "occupation": shell["occupation"],
        "epsilon_Ha": shell["epsilon"],
        "n_grid": int(n_grid),
        "r_max": float(r_max),
        "signed_alpha0_core": float(alpha_total[0]),
        "signed_C6_self_core": self_c6_from_alpha_grid(alpha_total, weights),
        "raw_positive_alpha0_core": float(positive_alpha_total[0]),
        "raw_positive_C6_self_core": self_c6_from_alpha_grid(positive_alpha_total, weights),
        "cphf_trk_scale": cphf_trk_scale,
        "cphf_alpha0_core": float(cphf_alpha_total[0]),
        "cphf_C6_self_core": self_c6_from_alpha_grid(cphf_alpha_total, weights),
        "signed_sum_osc": signed_sum_osc,
        "signed_sum_osc_over_occupation": signed_ratio,
        "positive_sum_osc": positive_sum_osc,
        "positive_sum_osc_over_occupation": positive_ratio,
        "cphf_positive_sum_osc": cphf_positive_sum_osc,
        "cphf_positive_sum_osc_over_occupation": cphf_positive_ratio,
        "negative_sum_osc": negative_sum_osc,
        "n_negative_delta_states": n_negative_delta,
        "operator_trk_pass": str(operator_trk_pass).lower(),
        "raw_absorption_trk_pass": str(raw_absorption_trk_pass).lower(),
        "cphf_absorption_trk_pass": str(cphf_absorption_trk_pass).lower(),
        "psp_go_no_go_ready": str(operator_trk_pass and cphf_absorption_trk_pass).lower(),
        "method": "PAULI_PROJECTED_TRK_CONSTRAINED_CPHF_EXPERIMENTAL",
        "note": "Radial-box central-field response with Pauli-positive absorption separated and constrained to the shell TRK sum before PSP C6 use.",
    }


def radial_grid_cphf_channels(atom, basis, shell_label, n_grid=700, r_max=60.0, min_osc=1e-14):
    shell = shell_radial_from_pyscf(atom, basis, shell_label, r_max=r_max, n_grid=n_grid)
    r, u, potential = invert_local_potential(shell["r"], shell["u"], shell["epsilon"], shell["l"])
    rows = []
    for l_resp in [shell["l"] - 1, shell["l"] + 1]:
        if l_resp < 0:
            continue
        hamiltonian = finite_difference_hamiltonian(r, l_resp, potential)
        rows.extend(
            positive_spectral_channels_from_hamiltonian(
                atom,
                shell_label,
                r,
                hamiltonian,
                shell["epsilon"],
                u,
                shell["occupation"],
                shell["l"],
                l_resp,
                min_osc=min_osc,
            )
        )

    positive_sum = sum(row["osc"] for row in rows)
    scale = shell["occupation"] / positive_sum if positive_sum > 0.0 else 1.0
    for row in rows:
        row["osc"] *= scale
        row["cphf_trk_scale"] = scale
    return rows


def write_channels(path, rows):
    fieldnames = [
        "atom",
        "channel",
        "delta_Ha",
        "osc",
        "is_core",
        "source",
        "from_shell",
        "response_l",
        "response_basis_index",
        "occ_energy_Ha",
        "response_energy_Ha",
        "radial_dipole",
        "cphf_trk_scale",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "atom": row["atom"],
                    "channel": row["channel"],
                    "delta_Ha": f"{row['delta_Ha']:.12f}",
                    "osc": f"{row['osc']:.12f}",
                    "is_core": row["is_core"],
                    "source": row["source"],
                    "from_shell": row["from_shell"],
                    "response_l": row["response_l"],
                    "response_basis_index": row["response_basis_index"],
                    "occ_energy_Ha": f"{row['occ_energy_Ha']:.12f}",
                    "response_energy_Ha": f"{row['response_energy_Ha']:.12f}",
                    "radial_dipole": f"{row['radial_dipole']:.12f}",
                    "cphf_trk_scale": f"{row['cphf_trk_scale']:.12f}",
                }
            )


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "basis",
        "shell",
        "response_l",
        "occupation",
        "epsilon_Ha",
        "n_grid",
        "r_max",
        "signed_alpha0_core",
        "signed_C6_self_core",
        "raw_positive_alpha0_core",
        "raw_positive_C6_self_core",
        "cphf_trk_scale",
        "cphf_alpha0_core",
        "cphf_C6_self_core",
        "signed_sum_osc",
        "signed_sum_osc_over_occupation",
        "positive_sum_osc",
        "positive_sum_osc_over_occupation",
        "cphf_positive_sum_osc",
        "cphf_positive_sum_osc_over_occupation",
        "negative_sum_osc",
        "n_negative_delta_states",
        "operator_trk_pass",
        "raw_absorption_trk_pass",
        "cphf_absorption_trk_pass",
        "psp_go_no_go_ready",
        "method",
        "note",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experimental radial-grid Sternheimer core response diagnostic.")
    parser.add_argument("--atom", required=True)
    parser.add_argument("--basis", default="ano-rcc")
    parser.add_argument("--shell", action="append", required=True)
    parser.add_argument("--n-grid", type=int, default=700)
    parser.add_argument("--r-max", type=float, default=60.0)
    parser.add_argument("--n-quad", type=int, default=120)
    parser.add_argument("--output", default="results/radial_grid_sternheimer/summary.csv")
    parser.add_argument("--channels-output")
    args = parser.parse_args(argv)

    rows = [
        radial_grid_sternheimer_summary(args.atom, args.basis, shell, args.n_grid, args.r_max, args.n_quad)
        for shell in args.shell
    ]
    write_rows(args.output, rows)
    if args.channels_output:
        channel_rows = []
        for shell in args.shell:
            channel_rows.extend(radial_grid_cphf_channels(args.atom, args.basis, shell, args.n_grid, args.r_max))
        write_channels(args.channels_output, channel_rows)
    print(
        "atom,basis,shell,response_l,cphf_alpha0_core,cphf_C6_self_core,"
        "signed_sum_osc_over_occupation,positive_sum_osc_over_occupation,"
        "cphf_positive_sum_osc_over_occupation,operator_trk_pass,psp_go_no_go_ready"
    )
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['basis']},"
            f"{row['shell']},"
            f"{row['response_l']},"
            f"{row['cphf_alpha0_core']},"
            f"{row['cphf_C6_self_core']},"
            f"{row['signed_sum_osc_over_occupation']},"
            f"{row['positive_sum_osc_over_occupation']},"
            f"{row['cphf_positive_sum_osc_over_occupation']},"
            f"{row['operator_trk_pass']},"
            f"{row['psp_go_no_go_ready']}"
        )
    print(f"\nWrote {args.output}")
    if args.channels_output:
        print(f"Wrote {args.channels_output}")


if __name__ == "__main__":
    main()
