import argparse
import csv
from pathlib import Path


def scalar_proxy_channel(row):
    delta = float(row["Delta_E_Ha"])
    occupation = float(row.get("occupation", 1.0))
    f0 = float(row["f0"])
    osc = occupation * (f0 / delta) ** 2
    return {
        "atom": row["atom"],
        "channel": f"eft_core_scalar_proxy_{row['shell']}",
        "delta_Ha": f"{delta:.12f}",
        "d2": "",
        "osc": f"{osc:.12f}",
        "is_core": "true",
        "source": "EFT_CORE_SCALAR_PROXY",
    }


def read_wilson_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_channels(path, rows):
    fieldnames = ["atom", "channel", "delta_Ha", "d2", "osc", "is_core", "source"]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_channels(input_path="results/core_wilson_coefficients.csv", output_path="results/eft_core_scalar_proxy_channels.csv"):
    rows = [scalar_proxy_channel(row) for row in read_wilson_rows(input_path)]
    write_channels(output_path, rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build scalar-proxy EFT-core oscillator channels from Wilson coefficients.")
    parser.add_argument("--input", default="results/core_wilson_coefficients.csv")
    parser.add_argument("--output", default="results/eft_core_scalar_proxy_channels.csv")
    args = parser.parse_args(argv)

    rows = build_channels(args.input, args.output)
    print("atom,channel,delta_Ha,osc,source")
    for row in rows:
        print(f"{row['atom']},{row['channel']},{row['delta_Ha']},{row['osc']},{row['source']}")


if __name__ == "__main__":
    main()
