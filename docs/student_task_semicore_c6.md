# Student task (≈1 week): non-empirical semicore dispersion for large-core PSPs

## The direction (why this is the promising one)

The plain vdW/C6 chain — monomer `α(iξ)` → Casimir–Polder `C6 = (3/π)∫α_A α_B dξ` — is
textbook for valence-only small molecules, and CCSD-response benchmarks already nail it to
~1%. There is **no contribution there**.

The open, defensible niche is **heavy elements**. A large-core pseudopotential / ECP freezes a
*polarizable semicore* shell and silently drops its contribution to the dynamic polarizability,
hence to dispersion. The standard fix is the **core-polarization potential (CPP)**, which carries
an **empirical cutoff radius `r_c`**. Our angle:

> compute the frozen-core dipole response `Δα_core(iξ)` **non-empirically** from the PSP
> partition — no cutoff — and validate it where the effect is large.

Scope it correctly (from our own controls): the effect is **negligible for deep cores**
(Be 1s → 2%, Kr → 4% closure) and for organics, and **large only for diffuse semicore shells**
(Mg 2s2p, Ca 3s3p closed the all-electron−PSP gap). The biggest and most novel case is the
**`d¹⁰` semicore of Zn / Cd**, which dominates the dispersion of the real, dispersion-bound
**Zn₂ / Cd₂** van der Waals dimers.

## What already exists — do not redo

- `vdw/eft_core_alpha_derivation.md`: the partition `α_total = α_valence^PSP + Δα_core`, the
  double-counting rules (only correct shells absent from the PSP valence space), and clean
  large-core results — **Mg q2: 638.6 → 647.6 ≈ all-e 647.59**; Ca q2 large.
- `tddft/c6_calc/`: the pyscf `α(iξ)` → Casimir–Polder pipeline.

Two known weaknesses to fix:
1. The current `Δα_core` reuses **all-electron MO dipole channels** — a stand-in, not a
   prediction.
2. It uses **sum-over-bound-states**, which undershot atomic Na `α_core` by ~70% because the
   oscillator strength sits in the **continuum** (`tddft/c6_calc/c6_hf_result_v2.md`).

## This week — make `Δα_core(iξ)` genuinely non-empirical, then validate

**Theory (with agent help, ~2–3 d).**
1. Write `Δα_core(iξ)` as the **frozen-core `l = 1` dipole polarizability on the imaginary
   axis**, summed over the occupied semicore shells that are in the PSP **core** but **not** in
   the valence.
2. Compute it by **core linear response (Sternheimer / coupled-perturbed)** — solve
   `(H_core − ε_c ± iξ)|δψ_c^±⟩ = −r̂|ψ_c⟩` and form
   `α_core(iξ) = −Σ_c ⟨ψ_c|r̂|δψ_c^+ + δψ_c^-⟩/3`.
   **Not** sum-over-bound-states: the Sternheimer solve implicitly includes the full spectrum,
   so it captures the continuum oscillator strength that the bound-state sum missed for Na.
   This is the single most important methodological fix.
3. State honestly the relation to the PRL scalar form factor `f_K^c`: the vdW object is a
   **two-point core response (a polarization)**, *not* the single-propagator band self-energy.
   Derive the response directly; do **not** assume it equals the bandwidth coefficient.

**Numerics (~2–3 d).**
4. Implement `Δα_core(iξ)` via the core Sternheimer solve in the existing pyscf pipeline.
5. Build `α(iξ) = α_valence^{ECP}(iξ) + Δα_core(iξ)` and `C6` by Casimir–Polder. Run **Sr, Zn,
   Cd**, using a **large-core ECP that freezes the target shell** (valence = outer `ns`; core
   correction = Sr `4s4p`, Zn `3d¹⁰`, Cd `4d¹⁰`). If only a small-core ECP is available the
   shell is already in valence → nothing to correct (the Ca q10 situation); pick the PSP so the
   shell is frozen.

**Validation (~1 d).**
6. `C6(Zn₂)`, `C6(Cd₂)`, `C6(Sr₂)` with and without the correction, vs published values
   (Zn₂/Cd₂ are well-characterized dispersion dimers — look up CCSD(T) / spectroscopic `C6`;
   do not guess the numbers).
7. Checks: right **sign**; **closes** the all-electron−ECP gap; **no double-count** (target shell
   absent from valence); Sternheimer vs sum-over-states (confirm the continuum is now captured).

## Go / no-go criterion

Promising if, for **Zn/Cd**, the non-empirical core correction (a) has the right sign, (b) closes
most of the all-electron−ECP `C6` gap, and (c) brings `C6(Zn₂/Cd₂)` within ~10% of published —
**with no fitted cutoff**. If yes, this is a non-empirical alternative to CPP for heavy-element
dispersion, headlined by the `d¹⁰` metals, and worth writing up.

## Out of scope this week (follow-ons, once the atomic proof-of-concept holds)

Environmental screening of `α_core`; the core–valence cross term `α_cv`; the full screened
functional; periodic / solid-state systems.
