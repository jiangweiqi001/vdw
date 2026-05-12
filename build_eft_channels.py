import csv
import argparse
import numpy as np
from eft_alpha import dipole_vector_from_form_factor
from build_eft_channels_spectral import build_spectral_channels_from_csv, build_spectral_channels_from_json


def write_channels(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["atom", "channel", "delta_Ha", "d2", "osc", "is_core"],
        )
        writer.writeheader()
        writer.writerows(rows)


def build_channel_row(atom, channel, delta, d_vector, is_core=True):
    d2 = np.dot(d_vector.conjugate(), d_vector).real
    return {
        "atom": atom,
        "channel": channel,
        "delta_Ha": f"{delta:.8f}",
        "d2": f"{d2:.12f}",
        "osc": "",
        "is_core": str(is_core).lower(),
    }


def sample_dipole_form_factor(d_ref):
    """Return a form factor that encodes a dipole vector d_ref in the l=1 channel."""
    def form_factor(q):
        q = np.asarray(q, dtype=float)
        q2 = np.dot(q, q)
        if q2 == 0.0:
            return 0.0j
        tau = -1j * np.dot(q, d_ref) * np.exp(-0.5 * q2)
        return 4.0 * np.pi * tau / q2

    return form_factor


def build_sample_channels():
    """Generate sample EFT dipole channel rows for atoms.

    This is a builder for dipole-channel data. In the true EFT implementation,
    `form_factor` should be replaced by the core dipole form factor derived from
    the paper's Eq. (6)/(7) and then converted into `d2`.
    """
    rows = []
    atoms = [
        ("Na", "core_dipole_1", 0.078, np.array([4.36, 0.0, 0.0])),
        ("K", "core_dipole_1", 0.056, np.array([5.39, 0.0, 0.0])),
        ("Mg", "core_dipole_1", 0.200, np.array([2.96, 0.0, 0.0])),
        ("Ca", "core_dipole_1", 0.120, np.array([4.70, 0.0, 0.0])),
        ("Ar", "core_dipole_1", 0.260, np.array([1.16, 0.0, 0.0])),
    ]
    for atom, channel, delta, d_ref in atoms:
        form_factor = sample_dipole_form_factor(d_ref)
        d_num = dipole_vector_from_form_factor(form_factor)
        rows.append(build_channel_row(atom, channel, delta, d_num, is_core=True))
    return rows


def _read_reference_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def build_calibrated_channels(reference_path="reference_alpha_c6.csv"):
    """Build one effective oscillator per atom from reference alpha0 and self C6."""
    rows = []
    for ref in _read_reference_rows(reference_path):
        atom = ref["atom"].strip()
        alpha0 = ref.get("alpha0_ref", "").strip()
        c6 = ref.get("C6_ref", "").strip()
        if not atom or not alpha0 or not c6:
            continue

        alpha0 = float(alpha0)
        c6 = float(c6)
        if alpha0 <= 0.0 or c6 <= 0.0:
            raise ValueError(f"Reference alpha0 and C6 must be positive for {atom}.")

        omega_eff = 4.0 * c6 / (3.0 * alpha0**2)
        osc = alpha0 * omega_eff**2
        rows.append(
            {
                "atom": atom,
                "channel": "calibrated_single_oscillator",
                "delta_Ha": f"{omega_eff:.12f}",
                "d2": "",
                "osc": f"{osc:.12f}",
                "is_core": "true",
            }
        )
    return rows


def build_channels_for_mode(
    mode,
    reference_path="reference_alpha_c6.csv",
    spectral_input=None,
    spectral_orbitals="atomic_spectral_input.csv",
    spectral_dipoles="radial_dipoles.csv",
    spectral_residuals="residual_oscillators.csv",
    add_residual_oscillator=True,
):
    if mode == "toy":
        return build_sample_channels()
    if mode == "calibrated":
        return build_calibrated_channels(reference_path)
    if mode == "spectral":
        if spectral_input is not None:
            return build_spectral_channels_from_json(spectral_input)
        return build_spectral_channels_from_csv(
            spectral_orbitals,
            spectral_dipoles,
            add_residual_oscillator=add_residual_oscillator,
            residual_path=spectral_residuals,
        )
    raise ValueError(f"Unknown channel build mode: {mode}")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build EFT atomic dipole channel rows.")
    parser.add_argument(
        "--mode",
        choices=["toy", "calibrated", "spectral"],
        default="toy",
        help="Channel source: toy placeholder, calibrated single oscillator, or spectral orbital input.",
    )
    parser.add_argument("--reference", default="reference_alpha_c6.csv")
    parser.add_argument("--spectral-input", default=None)
    parser.add_argument("--spectral-orbitals", default="atomic_spectral_input.csv")
    parser.add_argument("--spectral-dipoles", default="radial_dipoles.csv")
    parser.add_argument("--spectral-residuals", default="residual_oscillators.csv")
    parser.add_argument(
        "--no-residual-oscillator",
        action="store_true",
        help="Disable TRK residual oscillator in spectral CSV mode.",
    )
    parser.add_argument("--output", default="atomic_channels.csv")
    args = parser.parse_args(argv)

    rows = build_channels_for_mode(
        args.mode,
        args.reference,
        args.spectral_input,
        args.spectral_orbitals,
        args.spectral_dipoles,
        args.spectral_residuals,
        not args.no_residual_oscillator,
    )
    write_channels(args.output, rows)
    print(f"Wrote {args.output} with {args.mode} dipole-channel rows.")
    for row in rows:
        value = f"osc={row['osc']}" if row["osc"] else f"d2={row['d2']}"
        print(f"{row['atom']} {row['channel']} delta={row['delta_Ha']} {value}")


if __name__ == "__main__":
    main()
