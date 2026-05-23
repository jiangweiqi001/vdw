# Note: 从 frozen-core EFT 到 KS 带宽重整化与 vdW/TDDFT 响应修正

这份 note 的目标不是写成最终论文，而是把理论逻辑从头到尾打通。核心问题是：

> frozen-core pseudopotential 把芯电子变成静态势；但线性响应、TDDFT、ACFDT、vdW 需要动态密度响应。EFT 要做的事，是把被冻掉的 core dynamic response 作为 Wilson coefficient 匹配回来。

同一个 core dynamic response 有两个投影：

1. 投影到价电子单粒子传播上，给出 EFT-KS 的 $z_{\rm core}$ 和 quasiparticle band narrowing。
2. 投影到 density-density response 上，取偶极长波极限，给出 $\alpha_{\rm core}(i\xi)$、$C_6$、screened vdW / ACFDT correction。

下面从作用量开始推。

---

## 0. 记号和物理图像

把电子场分成 valence 与 core：

$$
\psi = \psi_v + \psi_c.
$$

芯电子激发能记为 $\Delta E_c$，价电子低能尺度记为 $\epsilon_F$。对 alkali / alkaline-earth 这类体系，通常有

$$
\epsilon_F \ll \Delta E_c.
$$

于是 core 可以被积分掉。积分掉后不会只留下普通赝势，还会留下一个动态响应核。普通 pseudopotential 是这个核的静态部分；EFT correction 是它的动态部分。

最短物理图像：

> core electron 像一个高频振子。普通赝势只保留它的平均位置；EFT-KS 看这个高频振子如何拖慢价电子传播；EFT-vdW 看这个高频振子的瞬时偶极涨落如何与远处另一个中心的涨落相关起来。

---

## 1. 从全电子作用量开始

在 Born-Oppenheimer 近似下，固定离子位置 ($\mathbf R_A$)。电子的 Euclidean action 写成

$$
\begin{aligned}
S[\bar\psi,\psi] &= \int_0^\beta d\tau \,
\bar\psi
\left(\partial_\tau - \frac{\nabla^2}{2} + V_{\rm lat}-\mu\right)
\psi \\
&\quad+
\frac12
\int_0^\beta d\tau
\int d\mathbf r d\mathbf r'
\rho(\mathbf r,\tau)
 v(\mathbf r-\mathbf r')
\rho(\mathbf r',\tau).
\end{aligned}
$$

$$
\rho(\mathbf r,\tau)=\bar\psi(\mathbf r,\tau)\psi(\mathbf r,\tau),
\qquad
v(\mathbf q)=\frac{4\pi}{q^2}.
$$

为了研究响应函数，引入外部标量源 $\phi$，耦合到电子密度：

$$
S_\phi = S - \int d\tau d\mathbf r\, \phi(\mathbf r,\tau)\rho(\mathbf r,\tau).
$$

生成泛函为

$$
Z[\phi]=\int D[\bar\psi,\psi]\,e^{-S_\phi}.
$$

密度响应可以从 $W[\phi]=\ln Z[\phi]$ 的二阶变分得到：

$$
\chi(1,2)=\frac{\delta \langle \rho(1)\rangle}{\delta \phi(2)}
=\frac{\delta^2 W[\phi]}{\delta \phi(1)\delta \phi(2)}.
$$

这里 $1\equiv(\mathbf r_1,\tau_1)$。

---

## 2. core/valence 分块

把 Hilbert space 投影到 core 与 valence 子空间：

$$
P_c+P_v=1,
\qquad
\psi_c=P_c\psi,
\qquad
\psi_v=P_v\psi.
$$

形式上，作用量可以写成

$$
S[\psi_v,\psi_c;\phi] =
S_v[\psi_v;\phi]
+S_c[\psi_c;\Phi]
+S_{cv}[\psi_v,\psi_c].
$$

这里最关键的是：core 看到的不是单独的外场 $\phi$，而是外场加上价电子产生的 Coulomb 势。对离子中心 $A$，定义

$$
\Phi_A(\mathbf r,\tau) =
\phi_{\rm ext}(\mathbf r,\tau)
+ \int d\mathbf r'\,
v(\mathbf r-\mathbf r')\rho_v(\mathbf r',\tau).
$$

于是 core 的有效作用量由

$$
 e^{-S_{\rm core}^{\rm eff}[\Phi]}
=\int D[\bar\psi_c,\psi_c]\,e^{-S_c[\psi_c;\Phi]}.
$$

定义。

因为 core excitation 是高能自由度，可以对 $\Phi$ 做 cumulant expansion：

$$
\begin{aligned}
S_{\rm core}^{\rm eff}[\Phi] &=
S_{\rm core}^{(0)}
+ \int n_c^0(1)\Phi(1) \\
&\quad-
\frac12
\int \Phi(1)\chi_c(1,2)\Phi(2)
+O(\Phi^3).
\end{aligned}
$$

这一步非常重要：

* $\int n_c^0\Phi$ 是静态平均 core charge 对价电子的作用。
* $-\frac12\Phi\chi_c\Phi$ 是 core 的动态密度响应。
* 更高阶 cumulant 对应非线性响应、多 core-hole 激发等，在低能展开中被 $\epsilon_F/\Delta E_c$ 压低。

---

## 3. 普通 pseudopotential 从哪里来

把 $\Phi=\phi+v\rho_v$ 代入线性项：

$$
\begin{aligned}
\int n_c^0 \Phi &=
\int n_c^0 \phi
+ \int n_c^0 v \rho_v.
\end{aligned}
$$

第二项就是 core 静态电荷对 valence 的 Coulomb 势。再加上 core-valence exchange、正交化排斥、norm-conserving / Kleinman-Bylander projector 等静态匹配，就得到普通的 static pseudopotential：

$$
\begin{aligned}
V_{\rm PSP}^{\rm static} &=
V_{\rm core\ Hartree}
+V_{\rm core\ exchange}
+V_{\rm orthogonality}
+\cdots.
\end{aligned}
$$

所以 frozen-core pseudopotential 的本质是：

$$
{\text{nucleus + core electrons}}
\longrightarrow
V_{\rm PSP}^{\rm static}.
$$

它是静态低能散射的等效势。它并不自动包含完整的动态响应。

---

## 4. core 动态响应核

二次 cumulant 给出

$$
S_{\rm core}^{(2)}[\Phi] =
- \frac12
\int \Phi\,\chi_c\,\Phi.
$$

展开 $\Phi=\phi+v\rho_v$：

$$
\begin{aligned}
S_{\rm core}^{(2)} &=
- \frac12 \phi\chi_c\phi \\
&\quad- \int \rho_v\, v\chi_c\phi \\
&\quad- \frac12 \rho_v v\chi_c v \rho_v.
\end{aligned}
$$

三项分别表示：

1. core 对外场的直接密度响应；
2. valence-density 与外场通过 core response 混合；
3. valence-valence 相互作用被虚拟 core excitation 动态修正。

这就是 EFT-vdW 和 EFT-KS 的共同起点。

---

## 5. core response 的谱表示

对单个离子中心 $A$，定义 core transition density form factor

$$
\tau_\lambda^A(\mathbf q) =
\langle 0_A|\hat\rho_c(\mathbf q)|\lambda_A\rangle,
$$

其中 $\lambda$ 标记 core excitation channel，包括 shell、角动量、磁量子数等。激发能为

$$
\Delta_\lambda = E_\lambda - E_0 > 0.
$$

在 imaginary frequency 上，core density response 的谱表示是

$$
\chi_c^A(\mathbf q,\mathbf q';i\xi) =
\sum_\lambda
\frac{2\Delta_\lambda\,\tau_\lambda^A(\mathbf q)\tau_\lambda^{A*}(\mathbf q')}
{\xi^2+\Delta_\lambda^2}
 e^{-i(\mathbf q-\mathbf q')\cdot\mathbf R_A}.
$$

这个式子来自两个时间方向的虚激发：

$$
\frac{1}{i\xi+\Delta_\lambda}
+ \frac{1}{-i\xi+\Delta_\lambda}
= \frac{2\Delta_\lambda}{\xi^2+\Delta_\lambda^2}.
$$

所以一粒子 self-energy 中看到的是单边 pole；密度响应中看到的是偶的 Matsubara 组合。

---

## 6. 和 PRL 中 $f_K^c$ 的关系

价电子密度通过 Coulomb 相互作用耦合到 core transition density。因此定义 Coulomb-dressed transition form factor

$$
f_\lambda(\mathbf q) =
\frac{4\pi}{q^2}\tau_\lambda(\mathbf q).
$$

反过来，

$$
\boxed{
\tau_\lambda(\mathbf q) =
\frac{q^2}{4\pi}f_\lambda(\mathbf q)
}
$$

于是 core-induced density kernel 可以写成

$$
v\chi_c v =
\sum_\lambda
\frac{2\Delta_\lambda\,f_\lambda(\mathbf q)f_\lambda^*(\mathbf q')}
{\xi^2+\Delta_\lambda^2}.
$$

EFT-KS 文章中的动态 pseudopotential 写作

$$
V_{\rm dyn}(K,K';i\omega) =
\sum_c
\frac{f_K^c f_{K'}^c}{i\omega+\Delta E_c}.
$$

这就是同一个结构的单粒子投影。它说明：

> $f$ 不是经验 $C_6$ 参数，而是 core transition density 被 Coulomb dressing 后的 Wilson coefficient。

但是要小心：EFT-KS 中用于带宽重整化的 $f_K^c$ 主要是 scalar / $\ell=0$ channel。vdW 的 leading term 是 dipole-dipole fluctuation，需要 $\ell=1$ 的 dipole-resolved form factor。

所以正确说法是：

$$
\text{band narrowing: scalar projection of core pole},
$$

$$
\text{vdW: dipole projection of the same core response spectrum}.
$$

不能机械地把 PRL 的 scalar $f_K^c$ 直接当 $C_6$ 参数。

---

## 7. 单粒子投影：推导 $z_{\rm core}$

先看 EFT-KS 这一路。标准 KS 静态 Hamiltonian 是

$$
\hat H_{\rm KS} =
- \frac{\nabla^2}{2}
+V_{\rm PSP}^{\rm static}
+V_H[n]
+V_{\rm xc}[n].
$$

KS 本征态满足

$$
\hat H_{\rm KS}|\psi_{\nu\mathbf k}\rangle
=\epsilon_{\nu\mathbf k}^{\rm KS}|\psi_{\nu\mathbf k}\rangle.
$$

把 Bloch wave 展成平面波：

$$
\psi_{\nu\mathbf k}(\mathbf r) =
\sum_{\mathbf G}c_{\nu\mathbf k}(\mathbf G)
 e^{i(\mathbf k+\mathbf G)\cdot\mathbf r}.
$$

将 core dynamic self-energy 投影到 KS 态上，定义

$$
F_{\nu\mathbf k}^c =
\sum_{\mathbf G}
 c_{\nu\mathbf k}(\mathbf G)
 f^c_{|\mathbf k+\mathbf G|}.
$$

于是带内近似下，动态 self-energy 为

$$
\Sigma_{\nu\mathbf k}^{\rm dyn}(\omega) =
\sum_c
\frac{|F_{\nu\mathbf k}^c|^2}{\omega+\Delta E_c}.
$$

不过这里有一个 double-counting 问题。$\Sigma^{\rm dyn}$ 含有静态部分，而静态部分已经被 $V_{\rm PSP}^{\rm static}$ 和 $V_{\rm xc}$ 吸收。不能再加一次。

所以真正用于 quasiparticle 方程的是 subtract 之后的 self-energy：

$$
\Sigma_{\rm sub}(\omega) =
\Sigma^{\rm dyn}(\omega)-\Sigma^{\rm dyn}(\epsilon_F).
$$

对一个 core channel，令

$$
D=\epsilon_F+\Delta E_c,
\qquad
\delta=\omega-\epsilon_F.
$$

则

$$
\begin{aligned}
\Sigma_{\rm sub}(\omega) &=
|F|^2
\left[\frac{1}{\omega+\Delta E_c}-\frac{1}{\epsilon_F+\Delta E_c}\right] \\
&= -\frac{|F|^2\delta}{D(D+\delta)}.
\end{aligned}
$$

因为低能 quasiparticle 满足 $|\delta|\ll \Delta E_c$，所以

$$
\Sigma_{\rm sub}(\omega)
\approx
-\frac{|F|^2}{D^2}(\omega-\epsilon_F)
\approx
-\frac{|F|^2}{\Delta E_c^2}(\omega-\epsilon_F).
$$

Dyson 方程为

$$
\omega-\epsilon_{\nu\mathbf k}^{\rm KS}
-\Sigma_{\rm sub}(\omega)=0.
$$

令

$$
s=\epsilon_{\nu\mathbf k}^{\rm KS}-\epsilon_F,
\qquad
\delta=\omega-\epsilon_F.
$$

代入得

$$
\delta-s+\frac{|F|^2}{\Delta E_c^2}\delta=0.
$$

因此

$$
\delta =
\frac{s}{1+|F|^2/\Delta E_c^2}.
$$

多个 core channel 相加，得到

$$
\boxed{
 z_{\nu\mathbf k}^{\rm core} =
\frac{1}{1+\sum_c |F_{\nu\mathbf k}^c|^2/\Delta E_c^2}
}
$$

以及

$$
\boxed{
\epsilon_{\nu\mathbf k}^{\rm QP}
\approx
\epsilon_F
+ z_{\nu\mathbf k}^{\rm core}
\left(\epsilon_{\nu\mathbf k}^{\rm KS}-\epsilon_F\right)
}
$$

这就是 band narrowing：所有能量相对费米能 $\epsilon_F$ 被压缩。

相应 Green's function 可写为

$$
G_{\nu\mathbf k}^{-1}(i\omega)
\approx
(z^{\rm val})^{-1}
\left[
(z_{\nu\mathbf k}^{\rm core})^{-1}i\omega
- (\epsilon_{\nu\mathbf k}^{\rm KS}-\epsilon_F)
\right].
$$

其中 $z^{\rm val}$ 是价电子 Fermi liquid residue。由于 Ward identity，价电子 residue 与 charge vertex 近似抵消，所以 $z^{\rm val}$ 不移动 pole；但是 core 已被积分掉，core polarizability 在低能下被 $\epsilon_F/\Delta E_c$ 压低，因此没有同样的 vertex cancellation，$z^{\rm core}$ 留在 quasiparticle dispersion 里。

---

## 8. 密度响应投影：从 $f$ 到 $\alpha(i\xi)$

vdW 不看单粒子 pole，而看 density fluctuation。对 core transition density $\tau_\lambda(\mathbf q)$，定义 transition dipole：

$$
 d_{\lambda,i} =
\int d\mathbf r\,r_i\rho_\lambda(\mathbf r).
$$

如果

$$
\tau_\lambda(\mathbf q) =
\int d\mathbf r\,e^{-i\mathbf q\cdot\mathbf r}\rho_\lambda(\mathbf r),
$$

则小 $q$ 展开给出

$$
\tau_\lambda(\mathbf q) =
\tau_\lambda(0) - i q_i d_{\lambda,i} + O(q^2).
$$

对中性闭壳层的偶极跃迁，$\tau_\lambda(0)=0$，所以

$$
 d_{\lambda,i}
= i\left.\frac{\partial \tau_\lambda(\mathbf q)}{\partial q_i}\right|_{\mathbf q=0}.
$$

用

$$
\tau_\lambda(\mathbf q) =
\frac{q^2}{4\pi}f_\lambda(\mathbf q),
$$

得到

$$
\boxed{
 d_{\lambda,i} =
 i\left.
\frac{\partial}{\partial q_i}
\left[\frac{q^2}{4\pi}f_\lambda(\mathbf q)\right]
\right|_{\mathbf q=0}
}
$$

这就是从 EFT Wilson coefficient 到偶极矩的桥。

然后 core polarizability tensor 是

$$
\boxed{
\alpha^A_{ij}(i\xi) =
\sum_{\lambda\in {\rm core}}
\frac{2\Delta_\lambda\,d^A_{\lambda,i}d^{A*}_{\lambda,j}}
{\Delta_\lambda^2+\xi^2}
}
$$

如果原子或闭壳层各向同性，则

$$
\boxed{
\alpha_A(i\xi) =
\frac13 {\rm Tr}\,\alpha^A_{ij}(i\xi)
= \frac23
\sum_{\lambda\in {\rm core}}
\frac{\Delta_\lambda |\mathbf d_\lambda^A|^2}
{\Delta_\lambda^2+\xi^2}
}
$$

所以从 EFT-KS 到 EFT-vdW 的核心链条是

$$
\boxed{
 f_\lambda(\mathbf q)
\Rightarrow
\tau_\lambda(\mathbf q)
\Rightarrow
\mathbf d_\lambda
\Rightarrow
\alpha_{\rm core}(i\xi)
\Rightarrow
C_6\ \text{or}\ E_{\rm vdW}
}
$$

---

## 9. 真空中两个中心的 London / Casimir-Polder 公式

若有两个孤立中心 $A,B$，其偶极涨落通过裸 Coulomb 相互作用传播。偶极-偶极张量为

$$
T_{AB,ij} =
\nabla_{R_A,i}\nabla_{R_B,j}\frac{1}{R_{AB}}.
$$

二阶涨落能量为

$$
E_{AB}^{(2)} =
-\frac{1}{2\pi}\int_0^\infty d\xi\,
{\rm Tr}
\left[
\alpha_A(i\xi)T_{AB}
\alpha_B(i\xi)T_{BA}
\right].
$$

如果 $A,B$ 都各向同性，

$$
\alpha_{A,ij}=\alpha_A\delta_{ij},
\qquad
\alpha_{B,ij}=\alpha_B\delta_{ij}.
$$

此时

$$
{\rm Tr}(T_{AB}T_{BA})=\frac{6}{R_{AB}^6}.
$$

于是

$$
E_{AB}^{(2)} =
-\frac{3}{\pi R_{AB}^6}\int_0^\infty d\xi\,
\alpha_A(i\xi)\alpha_B(i\xi).
$$

定义

$$
\boxed{
C_6^{AB} =
\frac{3}{\pi}\int_0^\infty d\xi\,
\alpha_A(i\xi)\alpha_B(i\xi)
}
$$

得到

$$
\boxed{
E_{AB}=-\frac{C_6^{AB}}{R_{AB}^6}
}
$$

这就是普通 vdW tail。

在 EFT-vdW 中，$\alpha_A$ 应该包含 pseudo-valence response 加上 matched core response：

$$
\alpha_A^{\rm target}(i\xi) =
\alpha_{A,v}^{\rm PSP}(i\xi)
+\Delta\alpha_{A,c}^{\rm EFT}(i\xi)
+\Delta\alpha_{A,cv}^{\rm EFT}(i\xi).
$$

最小版本可以先做

$$
\alpha_A^{\rm target}(i\xi) \approx
\alpha_{A,v}^{\rm PSP}(i\xi)
+ \alpha_{A,c}^{\rm EFT}(i\xi).
$$

但完整理论上还要处理 core-valence counterterm，避免 Pauli-forbidden transition 和已被 PSP/NLCC 隐含吸收的部分重复计算。

---

## 10. 固体中必须使用 valence-screened interaction

在材料里不能直接用裸 $v$，因为 core fluctuation 会被显式 valence electrons 屏蔽。

定义 valence irreducible response

$$
\chi_v^{\rm irr}(i\xi),
$$

它由 EFT tree-level KS propagator 和 matched charge vertex 计算。对应 valence-screened Coulomb interaction 为

$$
\boxed{
W_v(i\xi) =
\left[v^{-1}-\chi_v^{\rm irr}(i\xi)\right]^{-1}
}
$$

这一步的含义是：core-core fluctuation 不再通过裸 Coulomb 传播，而是通过被 valence electrons 和 local-field effect 屏蔽后的 $W_v$ 传播。

于是 core-induced nonlocal correlation energy 可以写成 RPA/log 形式：

$$
\boxed{
E_{\rm vdW}^{\rm core} =
\frac{1}{2\pi}
\int_0^\infty d\xi\,
{\rm Tr}_{\rm inter}
\left\{
\ln\left[1-W_v(i\xi)\chi_c^\perp(i\xi)\right]
+W_v(i\xi)\chi_c^\perp(i\xi)
\right\}
}
$$

这里：

* ${\rm Tr}_{\rm inter}$ 表示只取不同离子中心之间的非局域项；
* $\chi_c^\perp$ 表示去掉本地已匹配部分后的 core response；
* 加上线性项 $(+W_v\chi_c^\perp)$ 是为了去掉一阶 Hartree-like / local matching 贡献，使 correction 从二阶涨落开始。

展开 $\ln(1-X)+X$：

$$
\ln(1-X)+X=-\frac12X^2-\frac13X^3-\cdots.
$$

所以 leading 二阶项为

$$
\boxed{
E_{\rm vdW}^{(2)} =
-\frac{1}{4\pi}
\int_0^\infty d\xi\,
{\rm Tr}_{\rm inter}
\left[
W_v\chi_c^\perp W_v\chi_c^\perp
\right]
}
$$

对两个局域中心 $A,B$，做偶极近似，定义 screened dipole tensor

$$
T_{AB,ij}^{\rm scr}(i\xi) =
\nabla_{R_A,i}\nabla_{R_B,j}
W_v(\mathbf R_A,\mathbf R_B;i\xi).
$$

则

$$
\boxed{
E_{AB}^{(2)} =
-\frac{1}{2\pi}\int_0^\infty d\xi\,
{\rm Tr}
\left[
\alpha_A(i\xi)
T_{AB}^{\rm scr}(i\xi)
\alpha_B(i\xi)
T_{BA}^{\rm scr}(i\xi)
\right]
}
$$

在真空远距离极限，$W_v\to v$，恢复普通 $-C_6/R^6$。

---

## 11. double counting 规则

这是整个理论最容易错的地方。EFT-vdW 不是“在 DFT 上随便加一个 core $C_6$”。必须明确哪些东西已经在 PSP/DFT 里，哪些东西缺失。

### 11.1 不能重复加入静态 self-energy

单粒子修正中，动态 self-energy 必须 subtract：

$$
\Sigma_{\rm dyn}^{\rm sub}(i\omega) =
\Sigma_{\rm dyn}(i\omega)
-\Sigma_{\rm dyn}(i\omega_{\rm ref}).
$$

$\omega_{\rm ref}$ 通常取费米能附近。物理含义是：静态部分已经被 $V_{\rm PSP}^{\rm static}$ 和 $V_{\rm xc}$ 匹配掉了。

### 11.2 vdW 只保留 inter-center fluctuation

同一离子内部的 self-polarization、本地 core energy、以及线性 Hartree-like term 都不应该作为 vdW correction 加回来。定义

$$
{\rm Tr}_{\rm inter} =
{\rm Tr}-\sum_A {\rm Tr}_{AA}.
$$

这一步确保补的是 semilocal DFT 缺失的非局域长程 correlation，而不是重新拟合 atom energy 或 static pseudopotential。

### 11.3 显式 valence 中已有的 semicore 必须删掉

如果使用 small-core pseudopotential，把某些 semicore shell 显式放进 valence，那么这些 shell 的响应已经在 $\chi_v^{\rm irr}$ 里。它们不能再放进 $\chi_c$。

用 shell projector 写：

$$
\chi_c^{\rm EFT} =
\sum_{\lambda\in {\rm frozen\ core\ only}}
\chi_{c,\lambda}.
$$

### 11.4 本地匹配部分要投影掉

更一般地，定义本地投影 $P_{\rm loc}$，去掉已被 PSP、NLCC、local XC 或 atomic fitting 吸收的本地部分：

$$
\chi_c^\perp =
(1-P_{\rm loc})\chi_c(1-P_{\rm loc}).
$$

最终 correction 写成

$$
\boxed{
E_{\rm vdW}^{\rm EFT} =
\frac{1}{2\pi}\int_0^\infty d\xi\,
{\rm Tr}_{\rm inter}
\left\{
\ln\left[
1-W_v(i\xi)(1-P_{\rm loc})\chi_c(i\xi)(1-P_{\rm loc})
\right]
+W_v(i\xi)(1-P_{\rm loc})\chi_c(i\xi)(1-P_{\rm loc})
\right\}
}
$$

这就是 double-counting-free 的形式。

---

## 12. 和普通 atomic many-body vdW 的区别

atomic many-body vdW 文献已经会做

$$
\alpha(i\xi)=\alpha_v(i\xi)+\alpha_c(i\xi)+\alpha_{cv}(i\xi),
$$

并用

$$
C_6^{AB}=\frac{3}{\pi}\int_0^\infty d\xi\,
\alpha_A(i\xi)\alpha_B(i\xi).
$$

所以“core polarizability 会影响 $C_6$”不是新发现。

EFT 的新任务是不同的：

1. 给定一个 frozen-core pseudopotential，判断哪些 core physics 已经静态吸收，哪些 dynamic response 还缺失。
2. 从同一个 core Wilson coefficient 推出两个观测量：$z_{\rm core}$ 与 $\alpha_{\rm core}(i\xi)$。
3. 在材料中用 $W_v$ 做 valence-screened response，而不是孤立原子的裸 $C_6$。
4. 给 large-core / small-core / semicore partition 一个可诊断规则。

因此这个方向更准确地叫：

$$
\boxed{
\text{TDEFT for frozen-core pseudopotential response}
}
$$

而不只是“算一个 EFT-vdW $C_6$”。

---

## 13. 最小可计算流程

对每个元素：

1. 解 atomic core orbitals，得到 core shell、$\Delta_\lambda$、transition density。
2. 计算 multipole-resolved form factor $f_{\lambda,\ell m}(\mathbf q)$。
3. 对 $\ell=0$ scalar channel，投影到 KS Bloch 态，得到

$$
F_{\nu\mathbf k}^c =
\sum_{\mathbf G}c_{\nu\mathbf k}(\mathbf G)f^c_{|\mathbf k+\mathbf G|}
$$

和

$$
\boxed{
 z_{\nu\mathbf k}^{\rm core} =
\frac{1}{1+\sum_c |F_{\nu\mathbf k}^c|^2/\Delta E_c^2}
}
$$

4. 对 $\ell=1$ dipole channel，用

$$
 d_{\lambda,i} =
 i\partial_{q_i}\left[\frac{q^2}{4\pi}f_\lambda(\mathbf q)\right]_{q=0}
$$

得到

$$
\alpha_{ij}^{\rm core}(i\xi) =
\sum_\lambda
\frac{2\Delta_\lambda d_{\lambda,i}d_{\lambda,j}^*}
{\Delta_\lambda^2+\xi^2}.
$$

5. 对孤立原子/二聚体，用 Casimir-Polder 积分得到 $C_6$。
6. 对固体/表面/吸附，用 valence response 得到

$$
W_v=(v^{-1}-\chi_v^{\rm irr})^{-1},
$$

再算 screened London 或 full log-RPA correction。
7. 删除所有已显式在 valence 中的 semicore shell，删除 local/self terms，避免 double counting。

---

## 14. 这套理论真正要验证什么

最强验证不是只拟合一个 $C_6$。更强的是同一套 $\Delta_\lambda,f_\lambda$ 同时解释：

1. ARPES band narrowing：通过 $z_{\rm core}$。
2. static core polarizability：通过 $\alpha_{\rm core}(0)$。
3. imaginary-frequency tail：通过 $\alpha_{\rm core}(i\xi)$ 的高频 oscillator strength。
4. atomic / molecular $C_6$：通过 Casimir-Polder integral。
5. material screened vdW：通过 $W_v$-screened RPA/log formula。

如果这五个量能由同一套 core Wilson coefficient 同时给出，就说明 EFT matching 抓住了真正的 frozen-core dynamic physics，而不是又引入了一个经验 correction。

---

## 15. 一句话总结

普通 pseudopotential 完成的是

$$
\text{static matching:}\quad
{\text{nucleus + core}}\to V_{\rm PSP}^{\rm static}.
$$

TDEFT 要补的是

$$
\text{dynamic matching:}\quad
\chi_c(i\xi),\ f_\lambda(\mathbf q),\ \Delta_\lambda.
$$

同一个 dynamic core pole：

$$
\frac{f_\lambda f_\lambda^*}{i\omega+\Delta_\lambda}
$$

在单粒子通道里给

$$
 z_{\rm core},
$$

在密度响应通道里给

$$
\alpha_{\rm core}(i\xi),\quad C_6,\quad E_{\rm vdW}^{\rm ACFDT}.
$$

这就是从 EFT-KS 到 EFT-vdW/TDEFT 的完整理论骨架。
