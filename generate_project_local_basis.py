import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from generated_basis_protocol import sha256_file
from run_psp_rpa_atom import extract_cp2k_named_block


def generate_basis_from_cp2k_seed(
    atom,
    source_basis_file,
    source_basis_name,
    generated_basis_name,
    output_basis_file,
    provenance_file,
):
    block = extract_cp2k_named_block(source_basis_file, atom, source_basis_name)
    lines = block.splitlines()
    if not lines:
        raise ValueError(f"Empty CP2K basis block for {atom} {source_basis_name}.")

    lines[0] = f"{atom} {generated_basis_name}"
    output_basis_file = Path(output_basis_file)
    output_basis_file.parent.mkdir(parents=True, exist_ok=True)
    output_basis_file.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")

    record = {
        "atom": atom,
        "basis_provenance": "project_local_generated_from_cp2k_seed",
        "source_basis_file": str(source_basis_file),
        "source_basis_name": source_basis_name,
        "generated_basis_file": str(output_basis_file),
        "generated_basis_name": generated_basis_name,
        "generated_basis_sha256": sha256_file(output_basis_file),
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "generation_rule": "copy numeric CP2K seed block unchanged and replace header with project-local generated basis name",
        "vdw_targets_used": "false",
    }

    provenance_file = Path(provenance_file)
    provenance_file.parent.mkdir(parents=True, exist_ok=True)
    with open(provenance_file, "w", encoding="utf-8") as fp:
        json.dump(record, fp, indent=2, sort_keys=True)
        fp.write("\n")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate a project-local CP2K basis block from a fixed seed block.")
    parser.add_argument("--atom", required=True)
    parser.add_argument("--source-basis-file", required=True)
    parser.add_argument("--source-basis-name", required=True)
    parser.add_argument("--generated-basis-name", required=True)
    parser.add_argument("--output-basis-file", required=True)
    parser.add_argument("--provenance-file", required=True)
    args = parser.parse_args(argv)

    record = generate_basis_from_cp2k_seed(
        atom=args.atom,
        source_basis_file=args.source_basis_file,
        source_basis_name=args.source_basis_name,
        generated_basis_name=args.generated_basis_name,
        output_basis_file=args.output_basis_file,
        provenance_file=args.provenance_file,
    )
    print(f"wrote generated basis: {record['generated_basis_file']}")
    print(f"wrote provenance: {args.provenance_file}")
    return record


if __name__ == "__main__":
    main()
