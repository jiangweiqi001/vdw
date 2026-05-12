import csv
import numpy as np
from numpy.polynomial.legendre import leggauss

# ---------- polarizability and C6 ----------

def alpha_iw_from_osc(xi, delta, osc):
    """Compute alpha(i xi) in atomic units from oscillator strengths.

    Parameters
    ----------
    xi : array_like
        Imaginary frequencies in Hartree.
    delta : array_like
        Excitation energies in Hartree.
    osc : array_like
        Oscillator strengths.

    Returns
    -------
    alpha : ndarray
        Frequency-dependent polarizability values at each xi.
    """
    xi = np.atleast_1d(xi)[:, None]
    delta = np.asarray(delta)[None, :]
    osc = np.asarray(osc)[None, :]
    return np.sum(osc / (delta**2 + xi**2), axis=1)


def c6_from_alpha(delta_a, osc_a, delta_b, osc_b, n_quad=200):
    """Compute C6_AB from imaginary-frequency polarizabilities."""
    x, w = leggauss(n_quad)
    t = 0.25 * np.pi * (x + 1.0)
    xi = np.tan(t)
    jac = 0.25 * np.pi / np.cos(t) ** 2

    aa = alpha_iw_from_osc(xi, delta_a, osc_a)
    ab = alpha_iw_from_osc(xi, delta_b, osc_b)
    return (3.0 / np.pi) * np.sum(w * jac * aa * ab)


def london_single_oscillator(alpha0_a, omega_a, alpha0_b, omega_b):
    """Analytic C6 for a single effective oscillator model."""
    return 1.5 * alpha0_a * alpha0_b * omega_a * omega_b / (omega_a + omega_b)


def alpha0_from_osc(delta, osc):
    """Static polarizability alpha(0) from oscillator strengths."""
    return alpha_iw_from_osc([0.0], delta, osc)[0]


def self_c6_from_osc(delta, osc, n_quad=200):
    """Self C6 for a given atomic oscillator distribution."""
    return c6_from_alpha(delta, osc, delta, osc, n_quad=n_quad)


# ---------- EFT dipole form factor helper ----------

def dipole_vector_from_form_factor(form_factor, eps=1e-5):
    """Extract dipole vector d_lambda from Coulomb-dressed EFT form factor f(q).

    Uses:
        tau(q) = q^2/(4*pi) * f(q)
        d_i = i * partial_{q_i} tau(q)|_{q=0}

    form_factor(q) may be singular like 1/q near q=0; that is OK,
    because tau(q) should be regular.
    """
    def tau(q):
        q = np.asarray(q, dtype=float)
        q2 = np.dot(q, q)
        if q2 == 0.0:
            return 0.0j
        return q2 / (4.0 * np.pi) * form_factor(q)

    grad = np.zeros(3, dtype=complex)
    for i in range(3):
        dq_vec = np.zeros(3)
        dq_vec[i] = eps
        grad[i] = (tau(dq_vec) - tau(-dq_vec)) / (2.0 * eps)

    return 1j * grad


def oscillator_strength_from_dipole(delta, d_vector):
    """Convert dipole vector to oscillator strength in the dipole channel."""
    return (2.0 / 3.0) * delta * np.dot(d_vector.conjugate(), d_vector).real


def tau_from_form_factor(form_factor, q):
    """EFT dipole form factor tau(q) = q^2/(4pi) * f(q)."""
    q = np.asarray(q)
    q2 = np.dot(q, q)
    return q2 / (4.0 * np.pi) * form_factor(q)


def load_channels_csv(path):
    """Load EFT atomic channel data from a CSV file."""
    data = {}
    with open(path, newline="", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            atom = row["atom"].strip()
            delta_val = float(row["delta_Ha"])
            osc_val = row.get("osc", "")
            d2_val = row.get("d2", "")
            if atom not in data:
                data[atom] = {"delta": [], "osc": []}
            has_osc = osc_val.strip() != ""
            has_d2 = d2_val.strip() != ""
            if has_osc and has_d2:
                # Spectral builders may include d2 for auditability and osc with
                # shell occupation already folded in. The backend consumes osc.
                osc = float(osc_val)
            elif has_osc:
                osc = float(osc_val)
            elif has_d2:
                d2 = float(d2_val)
                osc = (2.0 / 3.0) * delta_val * d2
            else:
                raise ValueError(
                    f"Row for atom {atom} must contain osc or d2."
                )
            data[atom]["delta"].append(delta_val)
            data[atom]["osc"].append(osc)
    for atom in data:
        data[atom]["delta"] = np.asarray(data[atom]["delta"], dtype=float)
        data[atom]["osc"] = np.asarray(data[atom]["osc"], dtype=float)
    return data


# ---------- pairwise dispersion energy ----------

def pairwise_energy(positions_bohr, c6_matrix, damping=None):
    """Compute the pairwise -C6/R^6 dispersion energy."""
    e = 0.0
    n = len(positions_bohr)
    positions = np.asarray(positions_bohr, dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            r = np.linalg.norm(positions[i] - positions[j])
            fdamp = 1.0 if damping is None else damping(r, i, j)
            e -= fdamp * c6_matrix[i, j] / r**6
    return e


def fermi_damping(a, r0):
    """Return a simple Fermi-like damping function for R-dependent damping."""
    def damping(r, i=None, j=None):
        return 1.0 / (1.0 + np.exp(-a * (r / r0 - 1.0)))
    return damping


# ---------- convenience / sample utilities ----------

def effective_single_oscillator(delta, alpha0):
    """Return oscillator strength for a single effective oscillator."""
    return alpha0 * delta ** 2


def sample_atomic_parameters():
    """Return a sample dictionary of atomic effective oscillator parameters."""
    # Values are rough effective parameters in atomic units.
    return {
        "Na": {"alpha0": 162.7, "omega": 0.078},
        "K": {"alpha0": 290.0, "omega": 0.056},
        "Mg": {"alpha0": 71.3, "omega": 0.20},
        "Ca": {"alpha0": 157.0, "omega": 0.12},
        "Ar": {"alpha0": 11.1, "omega": 0.26},
    }


def build_oscillator_data(atom_params):
    """Build oscillator arrays for one-oscillator atoms."""
    delta = np.asarray([atom_params["omega"]])
    osc = np.asarray([effective_single_oscillator(atom_params["omega"], atom_params["alpha0"])])
    return delta, osc


def benchmark_c6():
    atoms = sample_atomic_parameters()
    names = list(atoms.keys())

    print("Sample effective single-oscillator atomic parameters:")
    for name in names:
        p = atoms[name]
        print(f"  {name}: alpha0={p['alpha0']:.2f} a.u., omega={p['omega']:.3f} Ha")

    print("\nC6 estimates using numerical integration and London analytic formula:")
    for i, a in enumerate(names):
        for b in names[i:]:
            delta_a, osc_a = build_oscillator_data(atoms[a])
            delta_b, osc_b = build_oscillator_data(atoms[b])
            c6_num = c6_from_alpha(delta_a, osc_a, delta_b, osc_b)
            c6_analytic = london_single_oscillator(
                atoms[a]["alpha0"], atoms[a]["omega"], atoms[b]["alpha0"], atoms[b]["omega"]
            )
            print(f"  {a}-{b}: C6_num={c6_num:.1f}, C6_analytic={c6_analytic:.1f}")


def sanity_check_single_oscillator():
    alpha0 = 10.0
    omega = 0.2
    delta = np.array([omega])
    osc = np.array([effective_single_oscillator(omega, alpha0)])
    c6_num = c6_from_alpha(delta, osc, delta, osc)
    c6_analytic = london_single_oscillator(alpha0, omega, alpha0, omega)
    print("Sanity check single oscillator:")
    print(f"  alpha0={alpha0}, omega={omega}")
    print(f"  C6 numerical = {c6_num:.8f}")
    print(f"  C6 analytic  = {c6_analytic:.8f}")
    print(f"  relative error = {abs(c6_num - c6_analytic) / c6_analytic:.3e}")


def sanity_check_dipole_extraction():
    d_ref = np.array([1.2, -0.4, 0.7])

    def form_factor(q):
        q = np.asarray(q, dtype=float)
        q2 = np.dot(q, q)
        if q2 == 0.0:
            return 0.0j

        tau = -1j * np.dot(q, d_ref) * np.exp(-0.5 * q2)
        return 4.0 * np.pi * tau / q2

    d_num = dipole_vector_from_form_factor(form_factor, eps=1e-5)
    print("Sanity check dipole extraction:")
    print(f"  d_ref = {d_ref}")
    print(f"  d_num = {d_num.real}")
    print(f"  error = {np.linalg.norm(d_num.real - d_ref):.3e}")


def main():
    sanity_check_single_oscillator()
    sanity_check_dipole_extraction()
    print("\nBenchmark sample atoms:\n")
    benchmark_c6()

    # Example pairwise energy for Ar2 at 8 bohr using bare C6
    ar = sample_atomic_parameters()["Ar"]
    delta_ar, osc_ar = build_oscillator_data(ar)
    c6_ar = c6_from_alpha(delta_ar, osc_ar, delta_ar, osc_ar)
    positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 8.0]])
    energy = pairwise_energy(positions, np.array([[0.0, c6_ar], [c6_ar, 0.0]]))
    print(f"\nAr2 bare -C6/R^6 at R=8 bohr: E={energy:.6f} Ha")


if __name__ == "__main__":
    main()
