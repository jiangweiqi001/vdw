import numpy as np
import csv
from eft_alpha import pairwise_energy


def load_c6_table(path="c6_table.csv"):
    c6 = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            a = row["A"]
            b = row["B"]
            val = float(row["C6_au"])
            c6[(a, b)] = val
            c6[(b, a)] = val
    return c6


def dimer_tail(atom_a, atom_b, r_min=10.0, r_max=30.0, n=11):
    c6 = load_c6_table()
    c6_ab = c6[(atom_a, atom_b)]

    print("A,B,R_bohr,E_Ha,E_meV,E_times_R6")
    for r in np.linspace(r_min, r_max, n):
        pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, r]])
        c6mat = np.array([[0.0, c6_ab], [c6_ab, 0.0]])
        e = pairwise_energy(pos, c6mat)
        print(f"{atom_a},{atom_b},{r:.6f},{e:.12e},{e*27211.386:.8f},{e*r**6:.8f}")


if __name__ == "__main__":
    dimer_tail("Na", "Na")
