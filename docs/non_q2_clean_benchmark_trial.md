# Non-q2 Clean Benchmark Trial

Be q2, Be q2 LDA, Kr q8, and Ca q10 test whether the additive dipole correction remains clean
when the corrected shells are absent from the explicit PSP valence space.

```text
case,C6_PSP,C6_PSP+EFT,C6_all_e,closure_pct,residual_C6,audit,status
Be_q2,209.84599821,210.87562603,256.41610946,2.210920,45.54048343,pass,clean_candidate
Be_q2_LDA,216.52565417,217.57226692,262.92582704,2.255623,45.35356012,pass,clean_candidate
Kr_q8,66.35553987,69.27712206,134.64958004,4.277946,65.37245798,pass,clean_candidate
Ca_q10,1424.67308305,1425.30863004,2206.75882541,0.081263,781.45019537,pass,clean_candidate
Ca_q2_PBE_adapted,1496.31224087,1658.03166276,2206.75882541,22.763066,548.72716265,pass,clean_candidate
```

Ca q10 is a deep-core-only diagnostic because 3s/3p/4s are explicit.