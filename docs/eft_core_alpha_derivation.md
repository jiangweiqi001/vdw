# Derivation Note: EFT Core Contribution to `alpha(i xi)`

This note records the working derivation target for Milestone C. It is not yet a
final theorem. Its purpose is to define the minimal correction that should be
implemented and tested against the all-electron and PSP-RPA baselines.

## 0. Milestone C Closure

The derivation is closed at the level needed for the prototype if the EFT-core
correction is formulated in terms of **response Wilson coefficients**, not by
reusing the scalar bandwidth coefficient directly.

The final working statement is:

```text
alpha_total(i xi)
  = alpha_valence^PSP-RPA(i xi)
  + alpha_core^EFT(i xi)
```

with

```text
alpha_core^EFT(i xi)
  = sum_lambda f_lambda^EFT / (Delta_lambda^2 + xi^2)
```

where `lambda` labels core/semicore dipole-response channels that are absent
from the explicit PSP valence space.

Equivalently,

```text
f_lambda^EFT = (2/3) Delta_lambda |d_lambda^EFT|^2
```

for an isotropic atom.

The dipole Wilson coefficient is defined by the long-wavelength density-response
form factor:

```text
d_lambda,i^EFT = i partial_{q_i} tau_lambda(q)|_{q=0}
```

and not by the scalar bandwidth `f_K^c` alone.

This closes the main ambiguity as follows:

1. **finite imaginary xi**: once the response channel `lambda` and its Wilson
   coefficient are matched, the finite-`xi` dependence is fixed by the spectral
   pole `1 / (Delta_lambda^2 + xi^2)`. No extra fitting of the `xi` dependence is
   introduced.
2. **scalar to dipole**: the scalar PRL coefficient and the dipole vdW
   coefficient are different Wilson coefficients in the same integrated-out-core
   EFT. The scalar coefficient controls the quasiparticle frequency derivative;
   the dipole coefficient controls density response to an external electric
   field. They share the same core scales and matching logic but are not
   algebraically identical.
3. **alpha_cv**: the first additive implementation does not introduce a separate
   `alpha_cv` oscillator. PSP-RPA is defined as the valence response in the
   static pseudopotential background. The Casimir-Polder integral of
   `alpha_valence + alpha_core` automatically includes the energy cross term
   `2 alpha_valence alpha_core`. Additional collective valence-core mixing beyond
   this additive response belongs to the future screened `W_v`/vertex-correction
   stage.
4. **double counting**: any shell explicitly present in the PSP valence space is
   projected out of `alpha_core^EFT`. The correction is added only for absent
   frozen shells.

The current code implements two approximations to this closed form:

```text
EFT_CORE_SCALAR_PROXY              # diagnostic only
EFT_CORE_DIPOLE_WILSON_MO_APPROX   # current l=1 dipole prototype
```

The second one is the physically relevant prototype because it uses a true
dipole channel. It is still an MO approximation to `d_lambda^EFT`, not the final
analytic multipole Wilson coefficient.

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

## 4.0 Multipole Derivation Of The `l = 1` Wilson Coefficient

This section closes the formal relation between a multipole-resolved
core-density form factor and the dipole oscillator channel used by the C6
backend.

Use the density convention

```text
rho(q) = integral d^3r exp(i q.r) rho(r).
```

For a core transition `lambda`, expand the transition density in spherical
harmonics:

```text
rho_lambda(r)
  = sum_lm rho_lambda,lm(r) Y_lm(rhat).
```

The Fourier-space transition density is

```text
tau_lambda(q)
  = <0|rho(q)|lambda>
  = 4 pi sum_lm i^l Y_lm(qhat)
      integral_0^infty dr r^2 j_l(q r) rho_lambda,lm(r).
```

The Coulomb-dressed multipole Wilson coefficient is

```text
F_lambda,lm(q)
  = v(q) tau_lambda,lm(q)
  = (4 pi / q^2) tau_lambda,lm(q).
```

For `l = 1`, use the small-argument limit

```text
j_1(q r) = q r / 3 + O(q^3).
```

Then

```text
tau_lambda,1m(q)
  = 4 pi i Y_1m(qhat) (q/3)
      integral_0^infty dr r^3 rho_lambda,1m(r)
    + O(q^3).
```

Therefore the Coulomb-dressed `l = 1` Wilson coefficient has the long-wavelength
singularity

```text
F_lambda,1m(q)
  = (16 pi^2 i / 3)
    Y_1m(qhat)
    [ integral dr r^3 rho_lambda,1m(r) ]
    / q
    + O(q).
```

The singularity is physical: a dipole transition density has no monopole
moment, so `tau(q) ~ q`, while the Coulomb dressing contributes `1/q^2`.

The undressed quantity that enters polarizability is recovered by removing the
Coulomb factor:

```text
tau_lambda(q) = q^2 F_lambda(q) / (4 pi).
```

The Cartesian dipole Wilson vector is therefore

```text
d_lambda,i^EFT
  = (1/i) partial_{q_i} tau_lambda(q)|_{q=0}
```

up to the sign convention for `exp(+-i q.r)`. The sign is irrelevant for C6
because the oscillator strength depends on `|d_lambda|^2`.

Equivalently, in spherical tensor notation,

```text
d_lambda,m^EFT
  = sqrt(4 pi / 3)
    integral_0^infty dr r^3 rho_lambda,1m(r).
```

This is the strict `l = 1` dipole Wilson coefficient. The scalar PRL coefficient
is the `l = 0` object:

```text
F_c,00(q)  -> bandwidth / z_core matching
F_lambda,1m(q) -> vdW dipole polarizability
```

They are different multipole projections of the integrated-out-core response.
They share core excitation energies and the same EFT matching logic, but the
dipole Wilson coefficient must be obtained from the `l = 1` transition density,
not from the scalar `l = 0` form factor.

The oscillator strength used by the backend is then

```text
f_lambda^osc
  = (2/3) Delta_lambda sum_i |d_lambda,i^EFT|^2
```

or, for a closed-shell isotropic atom, the equivalent spherical-tensor sum over
`m = -1,0,1`.

The current `compute_dipole_wilson.py` implements this structure with atomic MO
transition densities:

```text
rho_lambda(r) ~= phi_i*(r) phi_a(r)
```

for selected frozen shells. This is why its source label remains
`EFT_CORE_DIPOLE_WILSON_MO_APPROX`: it is a genuine `l=1` dipole channel, but it
is not yet an analytic core-only multipole Wilson coefficient computed directly
from a standalone core transition-density solver.

An intermediate improvement was implemented first:

```text
compute_core_tdhf_wilson.py
```

It computes the response of the isolated closed-shell core ion, e.g.

```text
Mg2+  -> Ne-like core
Ca2+  -> Ar-like core
```

and exports the resulting TDHF transition dipoles as

```text
source = EFT_CORE_DIPOLE_WILSON_CORE_TDHF
```

This is closer to the desired EFT object than the neutral-atom MO transition
approximation because the response is computed inside the core sector itself.
It is still a TDHF proxy for the core transition density, not an analytic
closed-form multipole Wilson coefficient.

Current clean q2 validation:

```text
atom  C6_PSP      C6_PSP+core-ion-TDHF  Delta_C6  status
Mg    638.6202    643.7337              +5.1136   clean
Ca    1972.5599   2135.8960             +163.3361 clean
```

The core-ion TDHF proxy gives the same sign as the MO dipole approximation and a
similar Ca magnitude, but a smaller Mg correction. This bracket is useful for
estimating uncertainty in the current dipole Wilson approximation.

The same transition-density contraction has now been made explicit:

```text
compute_multipole_core_wilson.py
```

For each TDHF excitation, the code forms the transition-density amplitude

```text
Gamma_lambda,ia = X_lambda,ia + Y_lambda,ia
```

and contracts it with the occupied-virtual dipole integrals:

```text
d_lambda = 2 sum_ia Gamma_lambda,ia <i|r|a>.
```

This is the discrete MO-basis version of

```text
d_lambda = integral r rho_lambda(r) d^3r
```

or equivalently the `q -> 0` derivative of `tau_lambda(q)`. The resulting
oscillator strengths exactly reproduce the PySCF TDHF oscillator-strength route
for the same core ion:

```text
Mg2+ sum_osc = 0.888436598687
Ca2+ sum_osc = 5.333290310908
```

The output is:

```text
results/eft_core_multipole_wilson_channels.csv
source = EFT_CORE_MULTIPOLE_TDENSITY_TDHF
```

This is the strictest implemented version of the `l=1` transition-density
Wilson channel in the current repo. It still obtains the transition density from
core-ion TDHF rather than from an analytic closed-form core solver, but it now
implements the multipole/small-q definition explicitly rather than relying on a
black-box oscillator-strength shortcut.

## 4.1 Detailed Additive Derivation

This section gives the closed additive derivation used by the current prototype.
It is deliberately limited to the unscreened atomic C6 level; screened `W_v`
appears only as the next extension.

### 4.1.1 Core integration with an external source

Introduce a scalar source `phi` coupled to the electronic density before
integrating out the frozen core. For one ionic center `A`, the core sees

```text
Phi_A(r, tau)
  = phi_ext(r, tau)
  + integral dr' v(r-r') rho_v(r', tau)
```

where `rho_v` is the valence density and `v(q) = 4 pi / q^2`.

After integrating out the core, the source-dependent effective action has the
cumulant expansion

```text
S_core^eff[Phi]
  = S_core^0
  + integral n_c^0 Phi
  - (1/2) integral Phi chi_c Phi
  + O(Phi^3)
```

The first two terms are static: they contribute to core charge screening,
exchange/orthogonality, and the fitted static pseudopotential. The quadratic
term is the dynamic core density response. Only the nonlocal dynamic part of
this quadratic response is eligible to generate a vdW correction.

### 4.1.2 Lehmann representation of the core response

Let `|lambda>` be an excited state of the isolated frozen-core sector with
excitation energy

```text
Delta_lambda = E_lambda - E_0 > 0.
```

Define the transition density form factor

```text
tau_lambda(q) = <0 | rho_c(q) | lambda>.
```

Then the imaginary-frequency density response of an isolated center is

```text
chi_c^A(q, q'; i xi)
  = sum_lambda
      2 Delta_lambda
      tau_lambda(q) tau_lambda*(q')
      / (Delta_lambda^2 + xi^2)
      exp[-i(q-q') . R_A].
```

This step resolves the finite-`xi` question for a matched response channel: once
`Delta_lambda` and `tau_lambda` are specified, the entire imaginary-axis
dependence is fixed by the spectral denominator. No additional finite-`xi` fit
is introduced.

### 4.1.3 Coulomb-dressed Wilson coefficient

The core response enters the valence problem through Coulomb coupling. Define
the Coulomb-dressed transition form factor

```text
F_lambda(q) = v(q) tau_lambda(q)
            = 4 pi tau_lambda(q) / q^2.
```

Equivalently,

```text
tau_lambda(q) = q^2 F_lambda(q) / (4 pi).
```

The PRL bandwidth coefficient `f_K^c` has the same Coulomb-dressed structure,
but it is the scalar channel relevant to the valence one-particle self-energy.
For vdW, the required object is the `l = 1` component of the response form
factor. We therefore denote the vdW coefficient by

```text
F_lambda,1m(q)
```

to avoid identifying it with the scalar bandwidth coefficient.

The dipole Wilson coefficient is the long-wavelength derivative of the
undressed transition density:

```text
d_lambda,i^EFT
  = i partial_{q_i} tau_lambda(q)|_{q=0}
  = i partial_{q_i} [q^2 F_lambda(q) / (4 pi)]_{q=0}.
```

This is the precise scalar-to-dipole resolution:

```text
scalar f_K^c        -> bandwidth z_core
dipole F_lambda,1m  -> alpha_core(i xi)
```

They are related by the same integrated-out-core EFT and the same core energy
scales, but they are not the same Wilson coefficient.

### 4.1.4 Dipole polarizability

The core contribution to the dipole polarizability tensor is

```text
alpha_c,ij(i xi)
  = sum_lambda
      2 Delta_lambda
      d_lambda,i d_lambda,j*
      / (Delta_lambda^2 + xi^2).
```

For an isotropic closed-shell atom,

```text
alpha_c(i xi)
  = (1/3) Tr alpha_c,ij(i xi)
  = sum_lambda f_lambda^EFT / (Delta_lambda^2 + xi^2),
```

with

```text
f_lambda^EFT
  = (2/3) Delta_lambda |d_lambda^EFT|^2.
```

This is the backend-compatible oscillator-channel form.

### 4.1.5 Additive PSP+EFT response and double counting

For a pseudopotential calculation, the baseline response is

```text
alpha_valence^PSP-RPA(i xi),
```

computed from the explicit valence electrons in the static pseudopotential
background.

The first EFT-vdW implementation uses the additive response

```text
alpha_EFT(i xi)
  = alpha_valence^PSP-RPA(i xi)
  + P_frozen alpha_core^EFT(i xi) P_frozen,
```

where `P_frozen` is a shell-space projector that keeps only shells absent from
the explicit PSP valence space.

In practice the projector is implemented by metadata:

```text
explicit_valence_shells ∩ eft_core_shells = empty.
```

If the intersection is non-empty, the correction is a double-counting
diagnostic, not a valid EFT-core benchmark.

### 4.1.6 What happens to alpha_cv?

The all-electron decomposition often contains

```text
alpha_all = alpha_v + alpha_c + alpha_cv.
```

In the additive PSP+EFT formulation, the explicit PSP-RPA response is defined as
the valence response in the static pseudopotential background. Thus the first
implementation assumes:

```text
alpha_cv,static-like -> encoded in PSP-valence orbitals and PSP-RPA baseline
alpha_cv,dynamic collective screening -> deferred to screened W_v / vertex stage
```

The Casimir-Polder energy formed from

```text
alpha_valence + alpha_core
```

already contains the energy-level cross term

```text
2 alpha_valence alpha_core
```

in the product `alpha_A alpha_B`. What it does not contain is a new independent
collective excitation where valence and core amplitudes mix dynamically inside
one atom. That missing physics is outside the additive prototype and belongs to
the screened/vertex-corrected theory.

### 4.1.7 Current implementable formula

The implemented benchmark formula is therefore:

```text
C6_PSP+EFT
  = (3/pi) integral_0^infty dxi
      [alpha_PSP(i xi) + alpha_core^EFT(i xi)]^2
```

with

```text
alpha_core^EFT(i xi)
  = sum_{lambda in frozen shells}
      f_lambda^EFT / (Delta_lambda^2 + xi^2).
```

This formula is the current Milestone C closure for atomic long-range C6. It is
not yet the screened EFT-vdW functional.

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
- Ca q2 is clean by shell overlap, but the available rows are basis diagnostics.
  The strongest adapted row uses `GTH-PBE-q2` with
  `TZV2P-MOLOPT-PBE-GTH-q10`, i.e. a q2 pseudopotential with a q10 UZH basis.
  This pseudo-basis mismatch must be reported explicitly; it is diagnostic only
  until a matched large-core Ca q2 basis is established.
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

## 9. Remaining Limitations

The basic additive core-response derivation above is sufficient for the current
prototype, but several limitations remain before claiming a full screened
EFT-vdW functional:

1. The current code computes `d_lambda^EFT` with an all-electron MO dipole
   approximation. A fully analytic multipole Wilson coefficient derived directly
   from core form factors remains to be implemented.
2. The scalar bandwidth `f_K^c` and the dipole response Wilson coefficient are
   distinct. Cross-observable validation should compare their shell trends, not
   assume equality.
3. PSP-RPA plus additive `alpha_core^EFT` includes the Casimir-Polder energy
   cross term but not full collective valence-core screening. That belongs to the
   future `W_v` stage.
4. The controlled error parameter for `C6` may differ from the bandwidth
   parameter because the relevant imaginary frequencies can be larger than the
   valence Fermi scale.

These limitations do not block the current PSP+EFT-core benchmark loop. They do
block any claim that the present code is a final screened EFT-vdW functional.
