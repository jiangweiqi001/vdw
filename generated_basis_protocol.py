import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


FORBIDDEN_TARGET_PATTERNS = (
    "alpha0",
    "alpha_0",
    "polarizability closure",
    "c6",
    "closure",
    "vdw",
    "dispersion",
    "psp+eft",
    "semicore correction",
)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fp:
        for chunk in iter(lambda: fp.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_allowed_targets(targets):
    bad = []
    for target in targets:
        normalized = str(target).strip().lower()
        if any(pattern in normalized for pattern in FORBIDDEN_TARGET_PATTERNS):
            bad.append(target)
    if bad:
        raise ValueError(f"Forbidden vdW target in generated-basis protocol: {', '.join(bad)}")


def build_freeze_record(
    element,
    pseudo_name,
    basis_name,
    basis_path,
    generation_method,
    allowed_targets,
    notes="",
    frozen_at_utc=None,
):
    validate_allowed_targets(allowed_targets)
    basis_path = Path(basis_path)
    if not basis_path.exists():
        raise FileNotFoundError(basis_path)

    frozen_at_utc = frozen_at_utc or datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    return {
        "element": element,
        "pseudo_name": pseudo_name,
        "basis_name": basis_name,
        "basis_path": str(basis_path),
        "basis_sha256": sha256_file(basis_path),
        "generation_method": generation_method,
        "allowed_targets": ";".join(allowed_targets),
        "forbidden_targets_used": "false",
        "frozen_at_utc": frozen_at_utc,
        "frozen_before_vdw_validation": "true",
        "benchmark_label": "generated_protocol_frozen",
        "notes": notes,
    }


def write_freeze_record(path, record):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fp:
        json.dump(record, fp, indent=2, sort_keys=True)
        fp.write("\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Freeze a generated PSP basis before vdW validation.")
    parser.add_argument("--element", required=True)
    parser.add_argument("--pseudo", required=True)
    parser.add_argument("--basis-name", required=True)
    parser.add_argument("--basis-file", required=True)
    parser.add_argument("--generation-method", required=True)
    parser.add_argument("--allowed-target", action="append", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--note", default="")
    args = parser.parse_args(argv)

    record = build_freeze_record(
        element=args.element,
        pseudo_name=args.pseudo,
        basis_name=args.basis_name,
        basis_path=args.basis_file,
        generation_method=args.generation_method,
        allowed_targets=args.allowed_target,
        notes=args.note,
    )
    write_freeze_record(args.output, record)
    print(f"wrote freeze record: {args.output}")
    return record


if __name__ == "__main__":
    main()
