# Screened Pairwise `W_v` Prototype Status

This is the first minimal screened-pairwise layer. It does not implement an
ab initio valence-screened interaction yet.

Implemented:

```text
screened_pairwise_vdw.py
```

Available model interactions:

```text
bare        W(r) = 1/r
dielectric  W(r) = 1/(epsilon r)
yukawa      W(r) = exp(-kappa r)/r
```

The isotropic pair energy is evaluated from the screened dipole tensor:

```text
E_AB = - C6 * Tr[T_scr T_scr] / 6
```

For `bare`, this reduces exactly to:

```text
E_AB = -C6/R^6
```

Generated example outputs for the clean Mg q2 benchmark:

```text
results/mg_q2_screened_tail_bare.csv
results/mg_q2_screened_tail_dielectric2.csv
results/mg_q2_screened_tail_yukawa02.csv
```

The same Mg q2 PSP+core channels were also passed through the finite-system
second-order logdet interface:

```text
results/mg_q2_eft_logdet_bare_R20.csv
results/mg_q2_eft_logdet_dielectric2_R20.csv
```

For the bare model at `R=20 Bohr`, the second-order logdet energy matches the
`-C6/R^6` tail. The dielectric `epsilon=2` result is reduced by `1/epsilon^2`,
as expected for two screened dipole tensors.

This stage is only a model-screening interface. It does not yet compute:

```text
W_v(i xi) = [v^-1 - chi_v^irr(i xi)]^-1
```

from an ab initio valence response, and it does not yet implement the full
finite-system or periodic log determinant.
