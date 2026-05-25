import argparse
import csv
from pathlib import Path

import numpy as np

from compute_multipole_core_wilson import compute_multipole_core_channels


def z_axis_tau_and_form_factor(atom, state_index, q, delta, dipole_z, oscillator_strength):
    """Small-q z-axis l=1 transition form factor.

    Along q = q zhat and density convention exp(i q.r), the leading transition
    density is tau(q) = i q d_z + O(q^3). The Coulomb-dressed coefficient is
    F(q) = 4 pi tau(q) / q^2.
    """
    q = float(q)
    dipole_z = float(dipole_z)
    tau = 1j * q * dipole_z
    form_factor = 4.0 * np.pi * tau / q**2
    return {
        "atom": atom,
        "state": int(state_index),
        "q": q,
        "delta_Ha": float(delta),
        "tau_q_real": float(tau.real),
        "tau_q_imag": float(tau.imag),
        "F_q_real": float(form_factor.real),
        "F_q_imag": float(form_factor.imag),
        "dipole_from_tau_z": float((tau / (1j * q)).real),
        "osc": float(oscillator_strength),
        "source": "EFT_CORE_MULTIPOLE_FORM_FACTOR_SMALL_Q",
    }


def build_form_factor_rows(channel_rows, q_values):
    rows = []
    for channel in channel_rows:
        for q in q_values:
            rows.append(
                z_axis_tau_and_form_factor(
                    atom=channel["atom"],
                    state_index=int(str(channel["channel"]).split("_")[-1]),
                    q=q,
                    delta=channel["delta_Ha"],
                    dipole_z=channel["d_z"],
                    oscillator_strength=channel["osc"],
                )
            )
    return rows


def write_rows(path, rows):
    fieldnames = [
        "atom",
        "state",
        "q",
        "delta_Ha",
        "tau_q_real",
        "tau_q_imag",
        "F_q_real",
        "F_q_imag",
        "dipole_from_tau_z",
        "osc",
        "source",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export small-q l=1 multipole form factors from core-ion TDHF dipoles.")
    parser.add_argument("--atom", default="Mg")
    parser.add_argument("--charge", type=int, default=2)
    parser.add_argument("--basis", default="aug-cc-pVQZ")
    parser.add_argument("--nstates", type=int, default=100)
    parser.add_argument("--q", action="append", type=float, default=[1e-4, 2e-4, 5e-4, 1e-3])
    parser.add_argument("--output", default="results/multipole_form_factors_mg_core.csv")
    args = parser.parse_args(argv)

    channels, _summary = compute_multipole_core_channels(args.atom, args.charge, args.basis, args.nstates)
    rows = build_form_factor_rows(channels, args.q)
    write_rows(args.output, rows)
    print("atom,state,q,delta_Ha,tau_q_imag,F_q_imag,dipole_from_tau_z,osc")
    for row in rows[:20]:
        print(
            f"{row['atom']},{row['state']},{row['q']:.8g},{row['delta_Ha']:.12f},"
            f"{row['tau_q_imag']:.12e},{row['F_q_imag']:.12e},"
            f"{row['dipole_from_tau_z']:.12f},{row['osc']:.12f}"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
