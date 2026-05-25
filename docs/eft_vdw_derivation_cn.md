# EFT-vdW 推导说明：从 PRL frozen-core 动力学到 vdW C6

这份文档把 PRL 中的 frozen-core EFT 思想和本项目的 vdW/C6 原型连接起来。

核心结论是：

```text
PRL 的 Eq. (6)/(7) 给出了 frozen core 的动态 Wilson coefficient。
vdW 需要把同一类 core 动态自由度写成 density response / dipole response。
当前项目已经实现了 additive PSP + l=1 core dipole channel 的原子 C6 原型，
但还没有完成 screened W_v 和 full log determinant。
```

## 1. PRL 已经给出的内容

PRL 的起点是把 frozen core 积分掉，得到 valence-only EFT。其核心结构是：

```text
L_val = psi^dagger (g_val)^-1 psi + L_int^val
```

其中静态部分给出普通 pseudopotential，而动态部分给出一个频率依赖 pole：

```text
(delta V_pp)^dyn_st(i omega)
  = sum_c M_sc M_tc^* / (i omega + Delta E_c)
```

在动量空间，它写成可分离形式：

```text
V_dyn(K,K'; i omega)
  = sum_c f_K^c f_K'^c / (i omega + Delta E_c)
```

其中 PRL Eq. (7) 给出 scalar Wilson coefficient：

```text
f_K^c
  = sqrt(4 pi) / K
    integral_0^infty u_c(r) [V_H,c(r) - J_c] sin(K r) dr
```

这里：

- `u_c(r) = r R_c(r)` 是 reduced core radial orbital
- `V_H,c` 是该 core orbital 的 Hartree potential
- `J_c` 是 self-Coulomb integral
- `Delta E_c` 是 core excitation energy

PRL 用这个动态 pole 的低频导数得到 quasiparticle bandwidth renormalization：

```text
z_core^{-1} - 1 ~ sum_c |F_c|^2 / Delta E_c^2
```

这解释了 Na/K/Ca/Mg 等体系的 Kohn-Sham bandwidth narrowing。

## 2. vdW 需要的不是一粒子 self-energy，而是 core response

要算 vdW，不能只看 valence one-particle self-energy。需要把 core 动态 pole 提升为 core density response。

在积分掉 core 之前，引入外源 `phi` 耦合到电子密度。一个离子中心 `A` 看到的总势是：

```text
Phi_A(r,tau)
  = phi_ext(r,tau)
  + integral dr' v(r-r') rho_v(r',tau)
```

对 core 自由度做 cumulant expansion：

```text
S_core^eff[Phi]
  = S_core^0
  + integral n_c^0 Phi
  - (1/2) integral Phi chi_c Phi
  + ...
```

其中：

- 线性项和静态项进入普通 pseudopotential
- 二次项 `chi_c` 是 core 的动态 density response
- vdW 来自不同中心之间的非局域 fluctuation，不是来自重新拟合静态 core energy

## 3. core response 的谱表示

令 `lambda` 表示一个 core excitation，能量为：

```text
Delta_lambda = E_lambda - E_0
```

定义 transition density form factor：

```text
tau_lambda(q) = <0 | rho_c(q) | lambda>
```

则 core density response 是：

```text
chi_c^A(q,q'; i xi)
  = sum_lambda
      2 Delta_lambda tau_lambda(q) tau_lambda^*(q')
      / (Delta_lambda^2 + xi^2)
      exp[-i(q-q') . R_A]
```

这一步回答了 finite imaginary frequency 的问题：

```text
一旦 channel lambda 和 Wilson coefficient 匹配好，
它在 imaginary axis 上的 xi 依赖就是
1 / (Delta_lambda^2 + xi^2)。
```

不需要额外拟合 `xi` 依赖。

## 4. Coulomb-dressed Wilson coefficient 与 dipole coefficient

core density 通过 Coulomb 与 valence density 耦合：

```text
F_lambda(q) = v(q) tau_lambda(q)
            = 4 pi tau_lambda(q) / q^2
```

等价地：

```text
tau_lambda(q) = q^2 F_lambda(q) / (4 pi)
```

注意：PRL 的 `f_K^c` 是 scalar channel，主要用于 bandwidth 的 `z_core`。vdW 的 leading term 是 dipole-dipole fluctuation，所以需要 `l=1` multipole-resolved coefficient。

因此必须区分：

```text
F_c,00(q)        -> bandwidth / z_core
F_lambda,1m(q)  -> vdW dipole polarizability
```

它们属于同一个 integrated-out-core EFT，但不是同一个 Wilson coefficient。

## 5. l=1 multipole Wilson coefficient

采用密度 Fourier convention：

```text
rho(q) = integral d^3r exp(i q.r) rho(r)
```

把 transition density 展开成球谐：

```text
rho_lambda(r)
  = sum_lm rho_lambda,lm(r) Y_lm(rhat)
```

其 Fourier form factor 是：

```text
tau_lambda(q)
  = 4 pi sum_lm i^l Y_lm(qhat)
      integral dr r^2 j_l(qr) rho_lambda,lm(r)
```

对于 `l=1`：

```text
j_1(qr) = qr/3 + O(q^3)
```

因此：

```text
tau_lambda,1m(q)
  = 4 pi i Y_1m(qhat) (q/3)
    integral dr r^3 rho_lambda,1m(r)
  + O(q^3)
```

Coulomb-dressed 形式为：

```text
F_lambda,1m(q)
  = (4 pi / q^2) tau_lambda,1m(q)
  ~ const * Y_1m(qhat) / q
```

这个 `1/q` 奇异性是正常的：dipole transition density 没有 monopole moment，所以 `tau(q) ~ q`，Coulomb dressing 给出 `1/q^2`。

真正进入 polarizability 的是 undressed transition density 的导数：

```text
d_lambda,i^EFT
  = i partial_{q_i} tau_lambda(q)|_{q=0}
```

等价的球张量形式是：

```text
d_lambda,m^EFT
  = sqrt(4 pi / 3)
    integral_0^infty dr r^3 rho_lambda,1m(r)
```

这就是严格的 `l=1` dipole Wilson coefficient。

## 6. 从 dipole Wilson 到 alpha_core(i xi)

core dipole polarizability tensor：

```text
alpha_c,ij(i xi)
  = sum_lambda
      2 Delta_lambda d_lambda,i d_lambda,j^*
      / (Delta_lambda^2 + xi^2)
```

闭壳层各向同性原子：

```text
alpha_c(i xi)
  = (1/3) Tr alpha_c,ij(i xi)
  = sum_lambda f_lambda^EFT / (Delta_lambda^2 + xi^2)
```

其中 oscillator strength 是：

```text
f_lambda^EFT
  = (2/3) Delta_lambda sum_i |d_lambda,i^EFT|^2
```

这正好接上当前代码后端：

```text
delta_Ha, osc -> alpha(i xi) -> C6
```

## 7. PSP + EFT-core 的加性原型

对于 pseudopotential valence response：

```text
alpha_valence^PSP-RPA(i xi)
```

第一版 additive EFT-vdW 原型写成：

```text
alpha_EFT(i xi)
  = alpha_valence^PSP-RPA(i xi)
  + P_frozen alpha_core^EFT(i xi) P_frozen
```

`P_frozen` 是 shell-space projector，实际代码规则是：

```text
explicit_valence_shells ∩ eft_core_shells = empty
```

如果某个 shell 已经显式在 PSP valence 里，就不能再作为 EFT-core correction 加回去。否则就是 double counting。

## 8. alpha_cv 的处理

全电子 response 常写成：

```text
alpha_all = alpha_v + alpha_c + alpha_cv
```

在当前 additive PSP+EFT 实现里：

```text
alpha_valence^PSP-RPA
```

定义为静态 pseudopotential 背景下的 valence response。因此：

```text
alpha_cv 的静态/平均场部分 -> PSP valence orbitals / PSP-RPA baseline
alpha_cv 的动态 collective mixing -> 未来 W_v / vertex stage
```

同时，Casimir-Polder 能量中：

```text
[alpha_valence + alpha_core]^2
```

自动包含能量层面的 cross term：

```text
2 alpha_valence alpha_core
```

但它不包含同一原子内部的 collective valence-core mixed excitation。这部分属于未来 screened/vertex-corrected EFT。

## 9. 当前代码对应关系

当前有三个层级：

```text
EFT_CORE_SCALAR_PROXY
```

只用 PRL scalar `f0/Delta`，是 diagnostic，不是 dipole vdW correction。

```text
EFT_CORE_DIPOLE_WILSON_MO_APPROX
```

使用 neutral atom MO transition dipoles，是真实 `l=1` dipole channel，但不是 core-only transition-density solver。

```text
EFT_CORE_MULTIPOLE_TDENSITY_TDHF
```

使用 isolated core ion TDHF 的 transition density，并显式做：

```text
Gamma_lambda,ia = X_lambda,ia + Y_lambda,ia
d_lambda = 2 sum_ia Gamma_lambda,ia <i|r|a>
```

这是目前代码中最接近严格 `l=1` multipole Wilson 的实现。

## 10. 当前 clean benchmark

最强 clean benchmark 是 Mg q2：

```text
PSP baseline:       Mg GTH-PBE-q2 / TZV2P-MOLOPT-SR-GTH-q2
explicit valence:   3s
EFT shells:         2s,2p
C6_PSP:             638.6202
C6_PSP+dipole_EFT:  647.6079
C6_all-e_PBE_TDDFT: 647.5881
double counting:    clean
```

这说明 additive `l=1` dipole core correction 在一个真正 clean q2 partition 上可以补回 PSP-to-all-electron C6 gap。

Ca q2 adapted 也有希望：

```text
C6_PSP:             1496.3122
C6_PSP+dipole_EFT:  1658.0317
C6_all-e_PBE_TDDFT: 2206.7588
closure:            about 22.8%
```

但它使用：

```text
GTH-PBE-q2 pseudo + TZV2P-MOLOPT-PBE-GTH-q10 basis
```

这是 pseudo-basis mismatch，所以只能作为 diagnostic。

## 11. 还没完成什么

还没完成：

```text
W_v(i xi) = [v^-1 - chi_v^irr(i xi)]^-1
```

还没完成：

```text
screened pairwise vdW with ab initio W_v
full log determinant / MBD-like energy
periodic implementation
forces
```

现在做的是：

```text
atomic long-range C6 additive PSP+EFT-core prototype
```

不是完整 screened EFT-vdW functional。

## 12. 下一步

下一步应该是：

1. 固化 Mg q2 clean benchmark。
2. 找 matched Ca/Sr q2 large-core basis。
3. 用 `EFT_CORE_MULTIPOLE_TDENSITY_TDHF` 替代/对照 `EFT_CORE_DIPOLE_WILSON_MO_APPROX`。
4. 再做 model `W_v`。
5. 最后才做 ab initio screened `W_v` 和 log determinant。
