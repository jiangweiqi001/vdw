# Derivation Note: EFT Core Contribution to `alpha(i xi)`

This note records the working derivation target for Milestone C. It is not yet a
final theorem. Its purpose is to define the minimal correction that should be
implemented and tested against the all-electron and PSP-RPA baselines.

## 1. Starting Point

The PRL EFT derivation integrates out frozen core states and gives the dynamic
part of the pseudopotential as a separable pole:

```text
delta V_pp^dyn(K, K'; i omega)
    = sum_c f_K^c f_K'^c / (i omega + Delta E_c)
```

The same coefficients enter the bandwidth renormalization:

```text
z_core^{-1} - 1 ~ sum_c |F_c|^2 / Delta E_c^2
```

For vdW, the target is instead a density response or dipole polarizability on
the imaginary axis:

```text
Delta alpha_core^EFT(i xi)
```

## 2. What Is Already Counted By A Static PSP

The integrated-out core contributes both static and dynamic terms.

Static pieces include:

- core Hartree screening of the nucleus
- core exchange/nonlocal orthogonality effects
- the zero-frequency/static part of the fitted pseudopotential

These are already part of a conventional norm-conserving pseudopotential and
must not be added again.

Therefore the EFT-vdW correction should not be a correction to the static
Hamiltonian. It should be a correction to the missing frequency-dependent core
response:

```text
alpha_total(i xi)
    = alpha_valence^PSP-RPA(i xi)
    + Delta alpha_core^EFT(i xi)
    + possible cross/mixing terms
```

The precise status of the cross/mixing term is an open derivation question. In a
strict PSP calculation, static core effects are already present in the valence
orbitals. The dynamic `alpha_cv` piece may be partly encoded in the PSP-RPA leg
and partly missing. This must be checked, not assumed.

## 3. From The Dynamic Pole To An Even-Frequency Response

The one-particle dynamic term has a pole:

```text
1 / (i omega + Delta E_c)
```

A density response or polarizability is even on the imaginary axis. The
corresponding oscillator kernel is obtained from the symmetric combination:

```text
1 / (i xi + Delta)
  + 1 / (-i xi + Delta)
  = 2 Delta / (xi^2 + Delta^2)
```

Thus any EFT core response channel should have the generic imaginary-frequency
shape:

```text
weight_c * 2 Delta_c / (xi^2 + Delta_c^2)
```

Equivalently, if represented in the existing oscillator-strength backend:

```text
Delta alpha_c(i xi) = f_c^osc / (Delta_c^2 + xi^2)
```

with an effective oscillator strength:

```text
f_c^osc = 2 Delta_c |d_c^eff|^2 / 3
```

for an isotropic dipole channel.

## 4. Relation Between `f_K^c` And Dipole Response

The PRL coefficient is a Coulomb-dressed form factor:

```text
f_K^c = sqrt(4 pi) / K * integral u_c(r) [V_H,c(r) - J_c] sin(K r) dr
```

For a density response derivation, introduce the transition density form factor:

```text
tau_lambda(q) = <0|rho_c(q)|lambda>
```

and its Coulomb-dressed version:

```text
f_lambda(q) = v(q) tau_lambda(q)
            = 4 pi tau_lambda(q) / q^2
```

If this identification is valid for the EFT Wilson coefficient, then:

```text
tau_lambda(q) = q^2 f_lambda(q) / (4 pi)
```

and the transition dipole is:

```text
d_lambda,i = i * partial_{q_i} tau_lambda(q)|_{q=0}
           = i * partial_{q_i} [q^2 f_lambda(q) / (4 pi)]_{q=0}
```

Important caveat:

The PRL `f_K^c` used for bandwidth is mainly an `s`-like scalar channel relevant
to quasiparticle frequency renormalization near the Fermi level. A leading
dipole vdW response requires the `l = 1` multipole-resolved version. The scalar
`f_K^c` is not automatically a dipole polarizability.

Therefore the implementation should distinguish two cases:

1. **Scalar bandwidth Wilson coefficient**: already used for `z_core`.
2. **Dipole EFT response coefficient**: still to be derived or computed as the
   `l = 1` density-response Wilson coefficient.

## 5. Minimal Implementable Model

The first code implementation should use a conservative oscillator-channel
model:

```text
Delta alpha_core^EFT(i xi)
    = sum_c f_c^EFT / (Delta_c^2 + xi^2)
```

where:

- `Delta_c` is the core or semicore excitation energy from an atomic calculation.
- `f_c^EFT` is a Wilson-derived oscillator strength or fitted diagnostic weight.
- the backend consumes `delta_Ha, osc` exactly like current TDHF/MO channels.

For validation, keep three labels separate:

```text
source = PySCF_TDHF       # all-electron response reference
source = PySCF_PSP_TDDFT  # PSP valence response
source = EFT_CORE_MODEL   # Wilson-derived core correction
```

No code should silently reuse all-electron TDHF core transition strengths and
call them EFT-core. That would be a diagnostic, not the proposed method.

## 6. Double-Counting Rules

When adding `Delta alpha_core^EFT` to `alpha_valence^PSP-RPA`:

1. Do not add static pseudopotential terms to the Hamiltonian.
2. Add only nonlocal/dynamic response pieces to the dispersion calculation.
3. If a shell is explicit in the PSP valence space, exclude that shell from the
   EFT-core correction.
4. Track the partition metadata:

```text
atom, psp, zval, explicit_valence_shells, eft_core_shells
```

For example:

```text
Ca large-core: explicit valence = 4s; EFT semicore may include 3s,3p
Ca small-core: explicit valence = 3s,3p,4s; EFT semicore must exclude 3s,3p
```

## 7. Validation Targets

The correction should be judged against:

```text
C6_all_e_RPA
C6_PSP_RPA
C6_PSP_RPA + Delta C6_EFT
C6_reference
```

Success criteria:

1. `C6_EFT - C6_PSP` has the right sign.
2. It closes a significant part of the all-electron minus PSP gap.
3. It does not overshoot because of double counting.
4. The same atomic Wilson coefficients are consistent with the PRL bandwidth
   `z_core` trends.

## 8. Immediate Coding Consequence

Before writing a screened functional, implement:

```text
compute_core_wilson.py
build_eft_core_correction.py
results/eft_core_validation_summary.csv
```

The first target should be Mg/Ca semicore shells because the PSP-RPA gap is
measurable and the PRL already identifies the outer `s` semicore channel as
important for bandwidth renormalization.

## 8.1 Implemented Diagnostic: Scalar Proxy

The first code pass implements only a diagnostic scalar proxy:

```text
compute_core_wilson.py
build_eft_core_correction.py
run_eft_core_proxy_validation.py
```

It computes the PRL scalar Wilson coefficient for selected s shells:

```text
Mg 2s
Ca 3s
```

using:

```text
f_K^c = sqrt(4 pi) / K * integral u_c(r) [V_H,c(r) - J_c] sin(K r) dr
```

with `Delta_E` approximated by a Koopmans-like all-electron HF orbital energy.

The diagnostic oscillator proxy is:

```text
osc_proxy = occupation * (f0 / Delta_E)^2
```

This is **not** the final EFT dipole correction. It uses the scalar bandwidth
coefficient as a proxy weight and is only meant to test direction and order of
magnitude.

Earlier small-core diagnostic result:

```text
atom  C6_PSP      C6_PSP+proxy  Delta_C6_proxy
Mg    637.3518    637.8663      +0.5146
Ca    1110.2731   1114.7195     +4.4464
```

The correction has the right sign for closing the PSP missing-C6 gap, but it is
far too small for Ca relative to the current all-electron-minus-PSP gap. This is
expected because the scalar `f0` proxy is not yet the dipole Wilson coefficient.
The next implementation must derive or compute the `l=1` dipole Wilson channel
rather than reusing the scalar bandwidth channel.

## 8.2 Implemented Diagnostic: l=1 Dipole Wilson Approximation

The next code pass implements a true dipole-channel approximation:

```text
compute_dipole_wilson.py
run_eft_core_dipole_validation.py
```

For selected core/semicore occupied shells, it computes ordinary 3D dipole
matrix elements:

```text
d_ia = <phi_i | r | phi_a>
```

and creates oscillator channels:

```text
f_ia = (2/3) (epsilon_a - epsilon_i) n_i |d_ia|^2
```

This is a real `l=1` dipole channel, unlike the scalar `f0` proxy. However, it is
still an **MO approximation** to the EFT dipole Wilson coefficient:

```text
source = EFT_CORE_DIPOLE_WILSON_MO_APPROX
```

It is not yet screened, and it uses all-electron atomic MO transition energies
and dipoles rather than a closed-form finite-xi Wilson coefficient derived from
the PRL scalar kernel.

Current diagnostic result:

```text
atom  C6_PSP      C6_PSP+dipole  Delta_C6_dipole  double_counting_status
Mg    637.3518    646.5804       +9.2286          diagnostic_double_counting
Ca    1424.6731   1596.6169      +171.9438        diagnostic_double_counting
```

Compared to the scalar proxy, the dipole channel has the expected much larger
effect. For Ca it recovers a noticeable fraction of the PSP missing-C6 gap while
remaining below the all-electron PBE-TDDFT value.

However, the current Mg and Ca PSPs are 10-electron small-core GTH
pseudopotentials. Their explicit valence spaces already include:

```text
Mg q10: 2s,2p,3s
Ca q10: 3s,3p,4s
```

Adding Mg `2s,2p` or Ca `3s,3p` dipole channels therefore double-counts
semicore response. The validation script fails by default in this situation and
only writes these rows when `--diagnostic-double-counting` is passed.

Interpretation:

- The sign is correct for the PSP+core correction test.
- The magnitude is now physically relevant for Ca.
- The result is still an unscreened additive test. It should not yet be called a
  final EFT-vdW correction until the relation between this MO dipole channel and
  the PRL Wilson coefficient is made explicit.
- Because the current Mg/Ca PSPs already contain the corrected semicore shells,
  these rows are diagnostics only. A clean test needs a true large-core PSP where
  the corrected shell is absent from the explicit valence space.

Current clean large-core validation result:

```text
atom  C6_PSP      C6_PSP+dipole  Delta_C6_dipole  explicit PSP shells  EFT shells  double_counting_status
Mg    638.6202    647.6079       +8.9878          3s                   2s,2p       clean
Ca    1972.5599   2149.5021      +176.9422        4s                   3s,3p       clean
```

Interpretation of the clean rows:

- Mg q2 is the first clean benchmark candidate: the PSP valence space is `3s`
  and the EFT correction shells are `2s,2p`, so the additive correction does not
  overlap the explicit PSP response. The current chain is
  `C6_PSP = 638.6202 -> C6_PSP+dipole = 647.6079 -> C6_all-e_PBE = 647.5881`.
- Ca q2 is clean by shell overlap, but the current row uses `GTH-LDA-q2` with
  `cc-pVQZ`; it should remain a clean diagnostic until a matched large-core Ca
  PSP basis route is established.
- These clean rows move the project from a pure double-counting diagnostic into
  a benchmark loop for at least Mg q2. The dipole channels are still an
  unscreened MO approximation to the EFT dipole Wilson coefficient, not the
  final screened EFT-vdW functional.

Additional clean negative controls now bracket this result:

```text
case        explicit PSP shells  EFT shells       closure
Be q2 PBE   2s                   1s               2.21%
Be q2 LDA   2s                   1s               2.26%
Kr q8 PBE   4s,4p                core through 3d  4.28%
Ca q10 PBE  3s,3p,4s             1s,2s,2p         0.081%
```

These controls pass the double-counting audit but do not close the all-electron
minus PSP gap. The Be LDA row is especially useful: because it uses
`GTH-LDA-q2` and LDA-TDDFT on both the PSP and all-electron legs, it shows that
Be is not a weak case because of a PBE/LDA mismatch. It is weak because the
frozen `1s` dipole contribution is small. Kr q8 and Ca q10 give the same lesson
for deeper frozen cores.

The second strong benchmark should therefore not be another deep-core clean
case. It should be a large-core semicore case:

```text
Ca q2: explicit valence = 4s; EFT semicore = 3s,3p
Sr q2: explicit valence = 5s; EFT semicore = 4s,4p
```

The immediate implementation task is to import or construct a matched q2 basis
for Ca or Sr and keep the basis provenance explicit in the benchmark metadata.

## 9. Open Questions

Open derivation questions to resolve before claiming the full EFT-vdW method:

1. Is the finite-`xi` coefficient identical to the bandwidth Wilson coefficient,
   or does the matching change away from `xi -> 0`?
2. What is the correct `l = 1` dipole version of the PRL scalar form factor?
3. Is `alpha_cv` captured by PSP-RPA, or does it need an explicit EFT cross term?
4. What is the controlled error parameter for `C6`, where relevant `xi` may be
   larger than the valence Fermi scale?

These questions should be answered in the derivation and checked numerically
before calling the correction a final screened EFT-vdW functional.
