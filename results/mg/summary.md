# Mg Semicore Partition Benchmark

## Partition

```text
valence = 3s
semicore = 2s,2p
deep_core = 1s
```

## TDHF Result

The Mg partition benchmark used `aug-cc-pVQZ` with `nstates=200`.

```text
C6_all                   = 757.83551633
C6_valence_3s            = 772.64572894
C6_valence_plus_semicore = 757.83785149
Delta_C6_semicore        = -14.80787745
relative contribution    = -1.95%
```

Interpretation: Mg has a smaller semicore/mixing correction than Ca. The
semicore effect is still visible, but it is closer to the noble-gas core/mixing
scale than to the stronger Ca semicore partition effect.
