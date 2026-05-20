import argparse
import csv
from pathlib import Path


def _read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def comparison_row(atom, all_e_c6, psp_c6, psp_status, partition):
    if psp_status != "ok" or psp_c6 in {"", None}:
        return {
            "atom": atom,
            "partition": partition,
            "C6_all_e": float(all_e_c6),
            "C6_psp": "",
            "delta_C6_missing": "",
            "relative_missing_pct": "",
            "status": psp_status,
        }
    all_e_c6 = float(all_e_c6)
    psp_c6 = float(psp_c6)
    delta = all_e_c6 - psp_c6
    return {
        "atom": atom,
        "partition": partition,
        "C6_all_e": all_e_c6,
        "C6_psp": psp_c6,
        "delta_C6_missing": delta,
        "relative_missing_pct": 100.0 * delta / all_e_c6 if all_e_c6 else 0.0,
        "status": psp_status,
    }


def _all_e_rows(path):
    rows = {}
    for row in _read_rows(path):
        if row["xc"] == "lda" and row["method"] == "TDDFT" and int(row["nstates"]) == 200:
            rows[row["atom"]] = row
        elif row["method"] == "TDDFT" and int(row["nstates"]) == 200:
            rows[(row["atom"], row["xc"])] = row
    return rows


def build_summary(all_e_path="results/all_e_rpa_summary.csv", psp_path="results/psp_rpa_summary.csv"):
    all_e = _all_e_rows(all_e_path)
    rows = []
    for psp in _read_rows(psp_path):
        atom = psp["atom"]
        all_e_row = all_e.get((atom, psp["xc"])) or all_e.get(atom)
        if all_e_row is None:
            continue
        rows.append(
            comparison_row(
                atom=atom,
                all_e_c6=all_e_row["C6"],
                psp_c6=psp["C6_psp"],
                psp_status=psp["status"],
                partition=f"{psp['psp']}/{psp['basis']}/{psp['xc']}",
            )
        )
    return rows


def write_summary(path, rows):
    fieldnames = [
        "atom",
        "partition",
        "C6_all_e",
        "C6_psp",
        "delta_C6_missing",
        "relative_missing_pct",
        "status",
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare all-electron LDA TDDFT C6 to PSP-RPA C6.")
    parser.add_argument("--all-e", default="results/all_e_rpa_summary.csv")
    parser.add_argument("--psp", default="results/psp_rpa_summary.csv")
    parser.add_argument("--output", default="results/all_e_vs_psp_rpa_summary.csv")
    args = parser.parse_args(argv)

    rows = build_summary(args.all_e, args.psp)
    write_summary(args.output, rows)
    print("atom,partition,C6_all_e,C6_psp,delta_C6_missing,relative_missing_pct,status")
    for row in rows:
        print(
            f"{row['atom']},"
            f"{row['partition']},"
            f"{row['C6_all_e']},"
            f"{row['C6_psp']},"
            f"{row['delta_C6_missing']},"
            f"{row['relative_missing_pct']},"
            f"{row['status']}"
        )


if __name__ == "__main__":
    main()
