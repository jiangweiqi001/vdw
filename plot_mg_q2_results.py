import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def plot_results(input_dir="results/mg_q2", output="results/mg_q2/benchmark.png"):
    input_dir = Path(input_dir)
    summary = read_rows(input_dir / "summary.csv")[0]
    tail = read_rows(input_dir / "tail.csv")

    labels = ["PSP", "PSP+EFT", "all-e PBE"]
    c6_values = [
        float(summary["C6_PSP"]),
        float(summary["C6_PSP_EFT"]),
        float(summary["C6_all_e"]),
    ]

    r = [float(row["R_bohr"]) for row in tail]
    e_all = [float(row["E_all_e_pbe_meV"]) for row in tail]
    e_psp = [float(row["E_psp_q2_meV"]) for row in tail]
    e_eft = [float(row["E_psp_plus_eft_meV"]) for row in tail]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(labels, c6_values, color=["#7aa6c2", "#83b692", "#333333"])
    axes[0].set_ylabel("C6 (a.u.)")
    axes[0].set_title("Mg q2 C6 closure")

    axes[1].plot(r, e_all, "k-", label="all-e PBE")
    axes[1].plot(r, e_psp, "o--", label="PSP q2")
    axes[1].plot(r, e_eft, "s--", label="PSP+EFT")
    axes[1].set_xlabel("R (Bohr)")
    axes[1].set_ylabel("E(R) (meV)")
    axes[1].set_title("Mg2 long-range tail")
    axes[1].legend()

    fig.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot Mg q2 clean benchmark results.")
    parser.add_argument("--input-dir", default="results/mg_q2")
    parser.add_argument("--output", default="results/mg_q2/benchmark.png")
    args = parser.parse_args(argv)

    print(plot_results(args.input_dir, args.output))


if __name__ == "__main__":
    main()
