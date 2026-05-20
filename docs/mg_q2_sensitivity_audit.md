# Mg q2 Sensitivity Audit

This audit checks whether the near-complete Mg q2 closure is stable under
basic numerical choices for the all-electron reference and PSP baseline.

## Summary

Total cases: 8
Review cases: 3

## Cases

```text
category,case,quantity,value,reference,delta_pct,status
all_e_nstates,aug-cc-pVQZ_nstates_100,C6_all_e,644.92232877,647.58810041,-0.411646,review
all_e_nstates,aug-cc-pVQZ_nstates_150,C6_all_e,645.37146889,647.58810041,-0.342290,review
all_e_nstates,aug-cc-pVQZ_nstates_200,C6_all_e,647.58810039,647.58810041,-0.000000,pass
all_e_basis,aug-cc-pVTZ_nstates_200,C6_all_e,661.15486234,647.58810041,2.094968,review
all_e_basis,aug-cc-pVQZ_nstates_200,C6_all_e,647.58810040,647.58810041,-0.000000,pass
psp_nstates,nstates_20,C6_psp,638.62015545,638.62015545,0.000000,pass
psp_nstates,nstates_50,C6_psp,638.62015545,638.62015545,0.000000,pass
psp_nstates,nstates_100,C6_psp,638.62015545,638.62015545,0.000000,pass
```

Cases marked `review` are not automatic failures; they identify numerical
choices that shift the benchmark beyond the configured tolerance and should
not be used as the headline value without explanation.
