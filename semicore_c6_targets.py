from dataclasses import dataclass


@dataclass(frozen=True)
class SemicoreC6Target:
    atom: str
    dimer: str
    explicit_shells: set
    correction_shells: set
    active_electrons: int
    note: str


SEMICORE_C6_TARGETS = {
    "Sr": SemicoreC6Target(
        atom="Sr",
        dimer="Sr2",
        explicit_shells={"5s"},
        correction_shells={"4s", "4p"},
        active_electrons=2,
        note="large-core alkaline-earth target; PSP valence should be outer 5s only",
    ),
    "Zn": SemicoreC6Target(
        atom="Zn",
        dimer="Zn2",
        explicit_shells={"4s"},
        correction_shells={"3d"},
        active_electrons=2,
        note="large-core d10 target; q12/small-core PSPs already include 3d10 and are not valid for this correction",
    ),
    "Cd": SemicoreC6Target(
        atom="Cd",
        dimer="Cd2",
        explicit_shells={"5s"},
        correction_shells={"4d"},
        active_electrons=2,
        note="large-core d10 target; q12/small-core PSPs already include 4d10 and are not valid for this correction",
    ),
}


def parse_shells(value):
    if isinstance(value, set):
        return set(value)
    return {shell.strip() for shell in str(value).replace(",", ";").split(";") if shell.strip()}


def audit_semicore_target(atom, active_electrons, active_shells):
    target = SEMICORE_C6_TARGETS[atom]
    active_shell_set = parse_shells(active_shells)
    overlap = active_shell_set & target.correction_shells
    active_count_ok = int(active_electrons) == int(target.active_electrons)
    active_shells_ok = active_shell_set == target.explicit_shells
    no_double_count = not overlap
    status = "pass" if active_count_ok and active_shells_ok and no_double_count else "fail"
    return {
        "atom": atom,
        "dimer": target.dimer,
        "active_electrons": int(active_electrons),
        "expected_active_electrons": target.active_electrons,
        "active_shells": ";".join(sorted(active_shell_set)),
        "expected_active_shells": ";".join(sorted(target.explicit_shells)),
        "correction_shells": ";".join(sorted(target.correction_shells)),
        "shell_overlap": ";".join(sorted(overlap)),
        "active_count_ok": str(active_count_ok).lower(),
        "active_shells_ok": str(active_shells_ok).lower(),
        "no_double_count": str(no_double_count).lower(),
        "audit_status": status,
        "note": target.note,
    }
