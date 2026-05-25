# Generated Basis Protocol Execution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the first executable pieces of the generated-basis protocol so Ca/Sr q2 basis work can be frozen and audited before vdW validation.

**Architecture:** Add a small metadata module plus CLI that creates a freeze record from a basis file and refuses forbidden vdW optimization targets. Keep basis generation itself separate from vdW validation; this first slice does not compute C6 or tune exponents.

**Tech Stack:** Python standard library, existing `unittest` tests in `tests/test_eft_pipeline.py`, CSV/JSON files for records.

---

### Task 1: Generated Basis Freeze Records

**Files:**
- Create: `generated_basis_protocol.py`
- Modify: `tests/test_eft_pipeline.py`

- [ ] **Step 1: Write failing tests**

Add tests that:
- create a temporary basis file
- build a freeze record with allowed non-vdW targets
- verify the record includes a SHA256 hash and `generated_protocol_frozen`
- verify forbidden targets such as `C6` and `closure` are rejected

- [ ] **Step 2: Run tests and verify RED**

Run:

```bash
python3 -m unittest tests.test_eft_pipeline.EftPipelineTests.test_generated_basis_freeze_record_includes_hash_and_protocol_label tests.test_eft_pipeline.EftPipelineTests.test_generated_basis_freeze_record_rejects_vdw_targets -v
```

Expected: fail because `generated_basis_protocol` does not exist.

- [ ] **Step 3: Implement minimal module**

Create `generated_basis_protocol.py` with:
- `FORBIDDEN_TARGET_PATTERNS`
- `validate_allowed_targets(targets)`
- `sha256_file(path)`
- `build_freeze_record(...)`
- `write_freeze_record(path, record)`

- [ ] **Step 4: Run tests and verify GREEN**

Run the same two tests and confirm they pass.

### Task 2: CLI Entry Point

**Files:**
- Modify: `generated_basis_protocol.py`
- Modify: `tests/test_eft_pipeline.py`

- [ ] **Step 1: Write failing CLI test**

Add a test that calls `main([...])` with a temporary basis file and freeze output path, then verifies JSON output exists.

- [ ] **Step 2: Run test and verify RED**

Expected: fail because `main` does not yet implement CLI writing.

- [ ] **Step 3: Implement CLI**

Add argparse options:
- `--element`
- `--pseudo`
- `--basis-name`
- `--basis-file`
- `--generation-method`
- `--allowed-target`
- `--output`
- `--note`

- [ ] **Step 4: Run targeted tests**

Run generated-basis tests and confirm they pass.

### Task 3: Documentation Link

**Files:**
- Modify: `docs/generated_basis_protocol.md`

- [ ] **Step 1: Add usage section**

Document the freeze-record command and expected output file role.

- [ ] **Step 2: Run placeholder scan**

Run a placeholder scan on `docs/generated_basis_protocol.md`.

Expected: no unresolved marker matches.
