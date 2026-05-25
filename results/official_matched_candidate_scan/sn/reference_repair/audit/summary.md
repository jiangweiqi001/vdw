# Sn ANO TDDFT100 vs TDA120 Oscillator Audit

## Files

```text
TDDFT100 channels = ../sn_all_e_ano_nstates100_channels.csv
TDA120 channels   = ../sn_all_e_ano_tda_nstates120_channels.csv
```

Audit tables:

```text
method_summary.csv
tddft100_channel_audit.csv
tda120_channel_audit.csv
tddft100_top_alpha.csv
tda120_top_alpha.csv
tddft100_top_c6.csv
tda120_top_c6.csv
```

## Method Summary

| method | channels | sum osc | alpha0 | C6 | min delta Ha |
|---|---:|---:|---:|---:|---:|
| TDDFT100 | 45 | 3.53813333 | 56.62325050 | 569.96394236 | 0.12716186 |
| TDA120 | 47 | 4.79466186 | 66.36259933 | 831.03236415 | 0.13306059 |

The TDA result is not larger because of a spurious lower-energy transition. Its
minimum excitation energy is slightly higher than TDDFT100. The main difference
is oscillator strength: TDA120 has about 35.5% more summed oscillator strength.

## Dominant Channels

TDDFT100 top alpha/C6 channel:

```text
channel = rpa_022
delta   = 0.242691659878 Ha
osc     = 0.911669722683
alpha fraction = 27.34%
C6 cross-inclusive fraction = 27.66%
```

TDA120 top alpha/C6 channel:

```text
channel = rpa_024
delta   = 0.259339014285 Ha
osc     = 1.244820272470
alpha fraction = 27.89%
C6 cross-inclusive fraction = 28.31%
```

Both methods are dominated by a similar low-frequency valence-response band
around `0.24-0.26 Ha`. TDA assigns substantially larger oscillator strength to
that band, and also increases strength in higher channels near `0.44 Ha`.

## Interpretation

The TDDFT/TDA discrepancy is a real reference-method sensitivity, not an obvious
single bad TDA pole.

Current implication for Sn q4:

```text
TDDFT100 reference -> closure about 96%
TDA120 reference   -> closure about 26%
```

The PSP-RPA side is already stable. The remaining uncertainty is the
all-electron Sn response model. Sn q4 should stay labeled:

```text
promising official-matched candidate, reference-limited
```

It should not be promoted to final benchmark until full TDDFT can be converged
past 100 states or an independently credible all-electron reference is obtained.
