# Generated Basis Protocol

## 1. Purpose

This document defines the protocol for generating project-local pseudopotential
basis sets, especially for candidate large-core `q2` cases such as Ca and Sr.

The goal is to allow generated basis sets without turning the van der Waals
benchmark into a fitted result. Basis generation and vdW validation must be
strictly separated.

The central rule is:

```text
The generated basis must not be optimized against alpha0, C6, PSP+EFT closure,
or any vdW energy result.
```

If a generated basis is later used in a PSP+EFT-core benchmark, the benchmark is
valid only if the basis was fixed before any vdW validation result was used.

## 2. Scope

This protocol applies to any generated or modified basis used for:

- large-core PSP-valence TDDFT/RPA calculations
- PSP+EFT-core or PSP+EFT-semicore C6 tests
- future screened or periodic vdW validation

It is especially relevant for cases where official matched CP2K/PySCF basis
sets are unavailable, for example:

```text
Ca q2: explicit valence = 4s; frozen semicore = 3s, 3p
Sr q2: explicit valence = 5s; frozen semicore = 4s, 4p
```

Generated basis sets must not be described as official MOLOPT, CP2K, or PySCF
matched basis sets.

## 3. Separation Of Stages

### Stage 1: Basis Generation

Allowed basis-generation targets are ordinary PSP electronic-structure targets:

- atomic total-energy convergence
- SCF convergence and numerical stability
- occupied orbital eigenvalue stability
- low-lying virtual orbital stability, if used only as a generic response-space
  quality check
- density or radial-tail stability
- transferability across simple charge states, such as `Ca`, `Ca+`, and `Ca2+`
- comparison to standard PSP atomic or small-molecule convergence behavior

Forbidden basis-generation targets are any quantities this project wants to
validate:

- `alpha0`
- `C6`
- `C6_PSP+EFT`
- PSP+EFT closure percentage
- long-range vdW tail energy
- agreement with all-electron TDDFT/RPA vdW response
- agreement with a target semicore correction

The basis-generation script or notes must record which targets were used.

### Stage 2: Freeze The Basis

Before any vdW validation is run, the generated basis must be frozen.

The frozen record should include:

- element and pseudopotential name
- generated basis name
- generation method and parameters
- source script or command
- date of generation
- file path
- file hash
- short note describing the allowed optimization targets

After this freeze point, exponents, contractions, or diffuse functions must not
be edited in response to vdW validation results.

If a change is needed later, it must create a new basis version with a new freeze
record. The previous vdW results must not be silently overwritten.

### Stage 3: vdW Validation

Only after the basis is frozen may the project compute:

- PSP-valence TDDFT/RPA response
- EFT-core or EFT-semicore oscillator corrections
- additive PSP+EFT C6
- all-electron comparison
- closure percentage
- pairwise or screened vdW tail diagnostics

The validation report must state that the generated basis was frozen before vdW
quantities were evaluated.

## 4. Required Audits

A generated basis should pass these checks before it is promoted from
diagnostic to benchmark status.

### 4.1 Basis Convergence

The basis should be part of a systematic sequence, for example:

```text
DZ-like -> TZ-like -> TZ2P-like -> QZ-like or extra-diffuse diagnostic
```

The benchmark should not depend on a single hand-selected exponent set.

### 4.2 Response Sanity Checks

For response calculations, check:

- oscillator-strength sum is not grossly inconsistent with the active electron
  count
- no unphysical near-zero excitation dominates `alpha0` or `C6`
- `alpha0` and `C6` are stable under reasonable basis enlargement
- adding diffuse functions improves or stabilizes the response rather than
  causing runaway growth
- PSP-valence response is reproducible with the frozen basis file

### 4.3 Double-Counting Audit

The basis must be paired with a clear shell partition.

For a large-core benchmark:

```text
explicit PSP valence shells  ∩  EFT correction shells = empty
```

For example:

```text
Ca q2:
explicit PSP valence shells = 4s
EFT correction shells       = 3s, 3p
```

If shell overlap exists, the result must be labeled diagnostic double counting,
not a clean benchmark.

### 4.4 Transferability

At minimum, check one or more independent states that were not used to tune vdW
closure, such as:

- neutral atom
- singly ionized atom
- doubly ionized atom
- simple dimer or small molecule, if available

The exact transferability tests may depend on the element, but they must be
declared before looking at vdW closure.

## 5. Naming Rules

Generated basis names must make provenance explicit.

Allowed examples:

```text
Ca-q2-generated-v1
Ca-q2-psp-energy-optimized-v1
Sr-q2-generated-v1
```

Disallowed examples:

```text
Ca TZV2P-MOLOPT-PBE-GTH-q2
Ca official q2 basis
Ca matched MOLOPT q2
```

Unless the basis comes from an official external source, it must be called
`generated`, `project-local`, or equivalent.

## 6. Benchmark Labels

Use these labels consistently:

```text
official_matched
    Official pseudo and official matching basis are both available.

generated_protocol_frozen
    Basis was generated by this protocol, frozen before vdW validation, and
    audited.

adapted_diagnostic
    Basis was borrowed from another q value or otherwise not generated by a
    fixed protocol.

invalid_for_prediction
    Response sanity checks fail, or the result violates the protocol.
```

Only `official_matched` and `generated_protocol_frozen` should be considered
clean benchmark categories.

## 7. Reporting Requirements

Any generated-basis benchmark page must include:

- basis provenance and freeze record
- pseudopotential name and source
- active PSP valence shell definition
- EFT correction shell definition
- statement that `alpha0`, `C6`, and closure were not used in basis generation
- basis convergence table
- response sanity checks
- double-counting status
- limitations

Recommended wording:

```text
This is a project-local generated q2 basis produced by a fixed protocol and
frozen before vdW validation. It is not an official CP2K/MOLOPT q2 basis and was
not fitted to alpha0, C6, or PSP+EFT closure.
```

## 8. Freeze Record Command

Use `generated_basis_protocol.py` to freeze a generated basis before any vdW
validation is run.

Example:

```bash
python3 generated_basis_protocol.py \
  --element Ca \
  --pseudo GTH-PBE-q2 \
  --basis-name Ca-q2-generated-v1 \
  --basis-file generated_basis/Ca-q2-generated-v1.bas \
  --generation-method "fixed non-vdW protocol" \
  --allowed-target "atomic total energy convergence" \
  --allowed-target "SCF robustness" \
  --allowed-target "orbital eigenvalue stability" \
  --output generated_basis/Ca-q2-generated-v1.freeze.json \
  --note "Frozen before PSP-RPA or C6 validation."
```

The command writes a JSON record containing:

- basis file hash
- pseudopotential name
- generated basis name
- generation method
- allowed non-vdW targets
- `benchmark_label = generated_protocol_frozen`
- `frozen_before_vdw_validation = true`

The command rejects forbidden targets such as `C6`, `alpha0`, `closure`, `vdW`,
or `PSP+EFT`.

## 9. Practical Rule For This Project

For the next Ca/Sr q2 attempt:

1. Write or select the basis-generation procedure.
2. Run only non-vdW generation and transferability checks.
3. Freeze the generated basis and record its hash.
4. Run PSP-RPA, EFT-core, and all-electron comparison.
5. Report the result whether closure is good or bad.

If the basis is changed after seeing vdW closure, the old validation result must
be treated as exploratory and a new frozen basis version must be created before
claiming a clean benchmark.
