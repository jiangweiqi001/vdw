import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from compute_dipole_wilson import compute_dipole_wilson_channels, write_channels as write_dipole_channels
from run_all_e_rpa_atom import export_all_e_response
from run_all_e_rpa_atom import run_atom as run_all_e_atom
from run_all_e_rpa_atom import write_channels as write_all_e_channels
from run_all_e_rpa_atom import write_summary as write_all_e_summary
from run_alpha_table import build_alpha_rows, write_alpha_rows
from run_c6_table import build_c6_rows, write_c6_rows
from run_eft_core_dipole_validation import c6_for, correction_shells, double_counting_status, write_combined
from run_psp_rpa_atom import run_atom as run_psp_atom
from run_psp_rpa_atom import write_summary as write_psp_summary


@dataclass(frozen=True)
class BenchmarkSpec:
    case_id: str
    atom: str
    psp: str
    psp_basis: str
    basis_file: str
    pseudo_file: str
    xc: str
    method: str
    psp_nstates: int
    all_e_basis: str
    all_e_nstates: int
    active_electrons: int
    explicit_shells: set
    correction_basis: str
    correction_shells: set
    role: str
    note: str


BENCHMARK_SPECS = {
    "Be_q2": BenchmarkSpec(
        case_id="Be_q2",
        atom="Be",
        psp="GTH-LDA-q2",
        psp_basis="TZV2P-MOLOPT-PBE-GTH-q2",
        basis_file="external_data/cp2k/BASIS_MOLOPT_UZH",
        pseudo_file="external_data/cp2k/GTH_POTENTIALS",
        xc="pbe",
        method="TDDFT",
        psp_nstates=100,
        all_e_basis="aug-cc-pVQZ",
        all_e_nstates=100,
        active_electrons=2,
        explicit_shells={"2s"},
        correction_basis="aug-cc-pVQZ",
        correction_shells={"1s"},
        role="light_alkaline_earth_q2_clean_diagnostic",
        note="Clean Be q2 diagnostic; local CP2K data has LDA/PADE q2 pseudo, not PBE q2 pseudo.",
    ),
    "Be_q2_LDA": BenchmarkSpec(
        case_id="Be_q2_LDA",
        atom="Be",
        psp="GTH-LDA-q2",
        psp_basis="TZV2P-MOLOPT-PBE-GTH-q2",
        basis_file="external_data/cp2k/BASIS_MOLOPT_UZH",
        pseudo_file="external_data/cp2k/GTH_POTENTIALS",
        xc="lda",
        method="TDDFT",
        psp_nstates=100,
        all_e_basis="aug-cc-pVQZ",
        all_e_nstates=100,
        active_electrons=2,
        explicit_shells={"2s"},
        correction_basis="aug-cc-pVQZ",
        correction_shells={"1s"},
        role="be_q2_lda_consistency_check",
        note="LDA consistency check for Be q2; uses GTH-LDA-q2 pseudo and LDA all-electron/PSP TDDFT.",
    ),
    "Kr_q8": BenchmarkSpec(
        case_id="Kr_q8",
        atom="Kr",
        psp="GTH-PBE-q8",
        psp_basis="TZV2P-MOLOPT-PBE-GTH-q8",
        basis_file="external_data/cp2k/BASIS_MOLOPT_UZH",
        pseudo_file="external_data/cp2k/GTH_POTENTIALS",
        xc="pbe",
        method="TDDFT",
        psp_nstates=100,
        all_e_basis="aug-cc-pVQZ",
        all_e_nstates=200,
        active_electrons=8,
        explicit_shells={"4s", "4p"},
        correction_basis="aug-cc-pVQZ",
        correction_shells={"1s", "2s", "2p", "3s", "3p", "3d"},
        role="noble_gas_core_clean_benchmark",
        note="Clean q8 noble-gas benchmark; correction excludes explicit 4s/4p PSP valence.",
    ),
    "Ca_q10": BenchmarkSpec(
        case_id="Ca_q10",
        atom="Ca",
        psp="GTH-PBE-q10",
        psp_basis="TZV2P-MOLOPT-PBE-GTH-q10",
        basis_file="external_data/cp2k/BASIS_MOLOPT_UZH",
        pseudo_file="external_data/cp2k/GTH_POTENTIALS",
        xc="pbe",
        method="TDDFT",
        psp_nstates=100,
        all_e_basis="cc-pVQZ",
        all_e_nstates=200,
        active_electrons=10,
        explicit_shells={"3s", "3p", "4s"},
        correction_basis="cc-pVQZ",
        correction_shells={"1s", "2s", "2p"},
        role="deep_core_only_clean_diagnostic",
        note="Clean deep-core diagnostic; Ca q10 keeps 3s/3p/4s explicit, so semicore EFT is excluded.",
    ),
    "Ca_q2_PBE_adapted": BenchmarkSpec(
        case_id="Ca_q2_PBE_adapted",
        atom="Ca",
        psp="GTH-PBE-q2",
        psp_basis="TZV2P-MOLOPT-PBE-GTH-q10",
        basis_file="external_data/cp2k/BASIS_MOLOPT_UZH",
        pseudo_file="external_data/cp2k/POTENTIAL_UZH_CASR_Q2",
        xc="pbe",
        method="TDDFT",
        psp_nstates=100,
        all_e_basis="cc-pVQZ",
        all_e_nstates=200,
        active_electrons=2,
        explicit_shells={"4s"},
        correction_basis="cc-pVQZ",
        correction_shells={"3s", "3p"},
        role="large_core_q2_adapted_basis_diagnostic",
        note="Ca q2 PBE diagnostic using official UZH GTH-PBE-q2 pseudo with q10 UZH basis as an explicitly adapted q2 candidate.",
    ),
}


def summarize_all_e_without_reference(atom, xc, basis, nstates, method, alpha_row):
    return {
        "atom": atom,
        "xc": xc,
        "basis": basis,
        "nstates": int(nstates),
        "method": method,
        "alpha0": float(alpha_row["alpha0_au"]),
        "C6": float(alpha_row["C6_self_au"]),
        "alpha0_ref": "",
        "alpha0_err": "",
        "C6_ref": "",
        "C6_err": "",
        "n_channels": int(alpha_row["n_channels"]),
    }


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def write_rows(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def upsert_row(path, row, keys, writer):
    existing = read_rows(path) if Path(path).exists() else []
    filtered = [old for old in existing if not all(str(old.get(key)) == str(row.get(key)) for key in keys)]
    filtered.append(row)
    writer(path, filtered)


def sorted_shells(value):
    if isinstance(value, set):
        return sorted(value)
    return sorted(shell for shell in str(value).split(";") if shell)


def select_one(rows, **criteria):
    matches = [
        row for row in rows
        if all(str(row.get(key)) == str(value) for key, value in criteria.items())
    ]
    if not matches:
        raise ValueError(f"No row matching {criteria}")
    if len(matches) > 1:
        raise ValueError(f"Multiple rows matching {criteria}")
    return matches[0]


def psp_channels_path(spec):
    leaf = f"{spec.psp}_{spec.psp_basis}_{spec.xc}_{spec.method.lower()}".replace("/", "_")
    return Path("results/psp_rpa") / spec.atom.lower() / leaf / f"{spec.atom.lower()}_psp_channels.csv"


def benchmark_row(spec, psp_row, all_e_row, eft_row):
    c6_psp = float(psp_row["C6_psp"])
    c6_all_e = float(all_e_row["C6"])
    c6_eft = float(eft_row["C6_psp_plus_dipole"])
    delta_psp = c6_all_e - c6_psp
    delta_eft = c6_all_e - c6_eft
    closure = (c6_eft - c6_psp) / delta_psp if delta_psp else 0.0
    explicit = set(sorted_shells(eft_row.get("psp_explicit_valence_shells", "")))
    corrected = set(sorted_shells(eft_row.get("correction_shells", "")))
    clean = eft_row.get("double_counting_status") == "clean" and corrected.isdisjoint(explicit)
    return {
        "case_id": spec.case_id,
        "atom": spec.atom,
        "psp": psp_row["psp"],
        "psp_basis": psp_row["basis"],
        "xc": psp_row["xc"],
        "method": psp_row["method"],
        "active_electrons": int(psp_row["active_electrons"]),
        "active_shells": psp_row["active_shells"],
        "eft_shells": ";".join(sorted_shells(eft_row["correction_shells"])),
        "C6_psp": c6_psp,
        "C6_psp_plus_eft": c6_eft,
        "C6_all_e": c6_all_e,
        "delta_C6_missing": delta_psp,
        "residual_C6": delta_eft,
        "closure_fraction": closure,
        "closure_pct": 100.0 * closure,
        "double_counting_status": eft_row["double_counting_status"],
        "benchmark_status": "clean_candidate" if clean else "needs_audit",
        "candidate_role": spec.role,
        "note": spec.note,
    }


def audit_row(spec, psp_path, psp_row, all_e_row, eft_row):
    correction = set(sorted_shells(eft_row.get("correction_shells", "")))
    explicit = set(sorted_shells(eft_row.get("psp_explicit_valence_shells", "")))
    shell_overlap = bool(correction & explicit)
    checks = {
        "active_electron_count": str(psp_row.get("active_electrons")) == str(spec.active_electrons),
        "active_shells_match": set(sorted_shells(psp_row.get("active_shells", ""))) == set(spec.explicit_shells),
        "correction_shells_match": correction == set(spec.correction_shells),
        "no_shell_overlap": not shell_overlap,
        "psp_pbe_tddft": psp_row.get("xc") == spec.xc and psp_row.get("method") == spec.method,
        "all_e_pbe_tddft": all_e_row.get("xc") == spec.xc and all_e_row.get("method") == spec.method,
        "placeholder_path_not_used": "placeholder" not in str(psp_path),
        "double_counting_clean": eft_row.get("double_counting_status") == "clean",
    }
    return {
        "case_id": spec.case_id,
        "atom": spec.atom,
        "psp_path": str(psp_path),
        "placeholder_path_used": str(not checks["placeholder_path_not_used"]).lower(),
        "shell_overlap": str(shell_overlap).lower(),
        "audit_status": "pass" if all(checks.values()) else "fail",
        **{key: str(value).lower() for key, value in checks.items()},
    }


def validate_case(spec, psp_path, correction_rows):
    output_root = Path("results/non_q2_clean_benchmarks") / spec.case_id
    combined_path = output_root / f"{spec.atom.lower()}_psp_plus_dipole_channels.csv"
    write_combined(combined_path, read_rows(psp_path) + correction_rows)
    c6_psp = c6_for(psp_path, output_root / "psp")[spec.atom]
    c6_total = c6_for(combined_path, output_root / "psp_plus_dipole")[spec.atom]
    corrected = correction_shells(correction_rows)
    status = double_counting_status(corrected, spec.explicit_shells)
    return {
        "case_id": spec.case_id,
        "atom": spec.atom,
        "C6_psp": c6_psp,
        "C6_psp_plus_dipole": c6_total,
        "Delta_C6_dipole": c6_total - c6_psp,
        "n_dipole_channels": len(correction_rows),
        "correction_shells": ";".join(sorted(corrected)),
        "psp_explicit_valence_shells": ";".join(sorted(spec.explicit_shells)),
        "double_counting_status": status,
        "note": "l=1 MO dipole Wilson approximation; unscreened additive test",
    }


def select_or_run_all_e(spec, rerun=False):
    path = "results/all_e_rpa_summary.csv"
    if not rerun and Path(path).exists():
        rows = read_rows(path)
        try:
            return select_one(
                rows,
                atom=spec.atom,
                xc=spec.xc,
                basis=spec.all_e_basis,
                method=spec.method,
                nstates=str(spec.all_e_nstates),
            )
        except ValueError:
            pass
    try:
        row = run_all_e_atom(spec.atom, spec.xc, spec.all_e_basis, spec.all_e_nstates, spec.method)
    except KeyError:
        row = run_all_e_without_reference(spec)
    upsert_row(path, row, keys=["atom", "xc", "basis", "nstates", "method"], writer=write_all_e_summary)
    return row


def run_all_e_without_reference(spec):
    output_dir = Path("results/all_e_rpa") / spec.atom.lower() / f"{spec.xc}_{spec.method.lower()}_{spec.all_e_basis.lower()}"
    output_dir.mkdir(parents=True, exist_ok=True)
    channels_path = output_dir / f"{spec.atom.lower()}_all_e_channels.csv"
    rows, meta = export_all_e_response(spec.atom, spec.xc, spec.all_e_basis, spec.all_e_nstates, spec.method)
    write_all_e_channels(channels_path, rows)
    alpha_rows = build_alpha_rows(channels_path)
    write_alpha_rows(output_dir / "alpha_c6_table.csv", alpha_rows)
    write_c6_rows(output_dir / "c6_table.csv", build_c6_rows(channels_path))
    alpha_row = next(row for row in alpha_rows if row["atom"] == spec.atom)
    return summarize_all_e_without_reference(
        atom=spec.atom,
        xc=spec.xc,
        basis=spec.all_e_basis,
        nstates=spec.all_e_nstates,
        method=meta["method"],
        alpha_row=alpha_row,
    )


def run_psp(spec):
    row = run_psp_atom(
        atom=spec.atom,
        psp=spec.psp,
        basis=spec.psp_basis,
        xc=spec.xc,
        nstates=spec.psp_nstates,
        method=spec.method,
        basis_file=spec.basis_file,
        basis_name=spec.psp_basis,
        pseudo_file=spec.pseudo_file,
        pseudo_name=spec.psp,
    )
    upsert_row(
        "results/psp_rpa_summary.csv",
        row,
        keys=["atom", "psp", "basis", "xc", "nstates", "method"],
        writer=write_psp_summary,
    )
    return row


def run_case(spec, rerun_all_e=False):
    psp_row = run_psp(spec)
    all_e_row = select_or_run_all_e(spec, rerun=rerun_all_e)
    correction_rows = compute_dipole_wilson_channels(spec.atom, spec.correction_basis, spec.correction_shells)
    dipole_path = Path("results/non_q2_clean_benchmarks") / spec.case_id / f"{spec.atom.lower()}_dipole_wilson_channels.csv"
    write_dipole_channels(dipole_path, correction_rows)
    eft_row = validate_case(spec, psp_channels_path(spec), correction_rows)
    return benchmark_row(spec, psp_row, all_e_row, eft_row), audit_row(spec, psp_channels_path(spec), psp_row, all_e_row, eft_row), eft_row


def write_markdown(path, benchmark_rows, audit_rows):
    lines = [
        "# Non-q2 Clean Benchmark Trial",
        "",
        "Be q2, Be q2 LDA, Kr q8, and Ca q10 test whether the additive dipole correction remains clean",
        "when the corrected shells are absent from the explicit PSP valence space.",
        "",
        "```text",
        "case,C6_PSP,C6_PSP+EFT,C6_all_e,closure_pct,residual_C6,audit,status",
    ]
    audits = {row["case_id"]: row for row in audit_rows}
    for row in benchmark_rows:
        lines.append(
            f"{row['case_id']},{row['C6_psp']:.8f},{row['C6_psp_plus_eft']:.8f},"
            f"{row['C6_all_e']:.8f},{row['closure_pct']:.6f},{row['residual_C6']:.8f},"
            f"{audits[row['case_id']]['audit_status']},{row['benchmark_status']}"
        )
    lines.extend(["```", ""])
    lines.append("Ca q10 is a deep-core-only diagnostic because 3s/3p/4s are explicit.")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def run_benchmarks(case_ids, rerun_all_e=False):
    benchmark_rows = []
    audit_rows = []
    eft_rows = []
    for case_id in case_ids:
        bench, audit, eft = run_case(BENCHMARK_SPECS[case_id], rerun_all_e=rerun_all_e)
        benchmark_rows.append(bench)
        audit_rows.append(audit)
        eft_rows.append(eft)
    write_rows(
        "results/non_q2_clean_benchmark_summary.csv",
        benchmark_rows,
        [
            "case_id", "atom", "psp", "psp_basis", "xc", "method", "active_electrons",
            "active_shells", "eft_shells", "C6_psp", "C6_psp_plus_eft", "C6_all_e",
            "delta_C6_missing", "residual_C6", "closure_fraction", "closure_pct",
            "double_counting_status", "benchmark_status", "candidate_role", "note",
        ],
    )
    write_rows(
        "results/non_q2_clean_benchmark_audit.csv",
        audit_rows,
        [
            "case_id", "atom", "psp_path", "placeholder_path_used", "shell_overlap",
            "audit_status", "active_electron_count", "active_shells_match",
            "correction_shells_match", "no_shell_overlap", "psp_pbe_tddft",
            "all_e_pbe_tddft", "placeholder_path_not_used", "double_counting_clean",
        ],
    )
    write_rows(
        "results/non_q2_clean_benchmark_eft_summary.csv",
        eft_rows,
        [
            "case_id", "atom", "C6_psp", "C6_psp_plus_dipole", "Delta_C6_dipole",
            "n_dipole_channels", "correction_shells", "psp_explicit_valence_shells",
            "double_counting_status", "note",
        ],
    )
    write_markdown("docs/non_q2_clean_benchmark_trial.md", benchmark_rows, audit_rows)
    return benchmark_rows, audit_rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="Try clean non-q2 PSP + EFT-core benchmark candidates.")
    parser.add_argument("--case", action="append", choices=sorted(BENCHMARK_SPECS), help="Case id to run; default runs all.")
    parser.add_argument("--rerun-all-e", action="store_true", help="Rerun all-electron PBE TDDFT references instead of reusing matching summary rows.")
    args = parser.parse_args(argv)
    case_ids = args.case or sorted(BENCHMARK_SPECS)
    bench_rows, audit_rows = run_benchmarks(case_ids, rerun_all_e=args.rerun_all_e)
    audits = {row["case_id"]: row for row in audit_rows}
    print("case_id,atom,C6_psp,C6_psp_plus_eft,C6_all_e,closure_fraction,residual_C6,audit_status,benchmark_status")
    for row in bench_rows:
        print(
            f"{row['case_id']},{row['atom']},{row['C6_psp']:.8f},{row['C6_psp_plus_eft']:.8f},"
            f"{row['C6_all_e']:.8f},{row['closure_fraction']:.12f},{row['residual_C6']:.8f},"
            f"{audits[row['case_id']]['audit_status']},{row['benchmark_status']}"
        )


if __name__ == "__main__":
    main()
