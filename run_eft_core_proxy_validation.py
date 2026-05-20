import argparse
import csv
from pathlib import Path

from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows


BEST_PSP = {
    "Mg": "results/psp_rpa/mg/gth-pbe_gth-qzv2p_pbe_tddft/mg_psp_channels.csv",
    "Ca": "results/psp_rpa/ca/gth-lda_gth-dzvp-molopt-sr_lda_tddft/ca_psp_channels.csv",
}


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_combined_channels(path, rows):
    fieldnames = ["atom", "channel", "delta_Ha", "osc", "is_core", "source"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})


def c6_for_channels(path, prefix):
    alpha_rows = build_alpha_rows(path)
    write_alpha_rows(f"{prefix}_alpha_c6_table.csv", alpha_rows)
    write_c6_rows(f"{prefix}_c6_table.csv", build_c6_rows(path))
    return {row["atom"]: float(row["C6_self_au"]) for row in alpha_rows}


def run_validation(proxy_path="results/eft_core_scalar_proxy_channels.csv", output="results/eft_core_validation_summary.csv"):
    proxy_rows = read_csv(proxy_path)
    rows = []
    for proxy in proxy_rows:
        atom = proxy["atom"]
        psp_path = BEST_PSP.get(atom)
        if not psp_path or not Path(psp_path).exists():
            continue
        combined_path = Path("results/eft_core_proxy") / atom.lower() / f"{atom.lower()}_psp_plus_proxy_channels.csv"
        combined = read_csv(psp_path) + [proxy]
        write_combined_channels(combined_path, combined)
        c6_psp = c6_for_channels(psp_path, Path("results/eft_core_proxy") / atom.lower() / "psp")[atom]
        c6_eft = c6_for_channels(combined_path, Path("results/eft_core_proxy") / atom.lower() / "psp_plus_proxy")[atom]
        rows.append(
            {
                "atom": atom,
                "C6_psp": c6_psp,
                "C6_psp_plus_proxy": c6_eft,
                "Delta_C6_proxy": c6_eft - c6_psp,
                "proxy_channel": proxy["channel"],
                "note": "scalar proxy from PRL f0/Delta; diagnostic not final EFT-core",
            }
        )
    with open(output, "w", newline="", encoding="utf-8") as fp:
        fieldnames = ["atom", "C6_psp", "C6_psp_plus_proxy", "Delta_C6_proxy", "proxy_channel", "note"]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate scalar-proxy EFT-core channels against PSP C6.")
    parser.add_argument("--proxy", default="results/eft_core_scalar_proxy_channels.csv")
    parser.add_argument("--output", default="results/eft_core_validation_summary.csv")
    args = parser.parse_args(argv)

    rows = run_validation(args.proxy, args.output)
    print("atom,C6_psp,C6_psp_plus_proxy,Delta_C6_proxy,proxy_channel,note")
    for row in rows:
        print(
            f"{row['atom']},{row['C6_psp']},{row['C6_psp_plus_proxy']},"
            f"{row['Delta_C6_proxy']},{row['proxy_channel']},{row['note']}"
        )


if __name__ == "__main__":
    main()
