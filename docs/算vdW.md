下面给一版**可直接补进稿子的 EFT-vdW 推导骨架**。核心点是：现有 Eq. (6)/(7) 已经给出了 core 动力学 pole 和 form factor；要算 vdW，需要把这个 pole 从“valence one-particle self-energy”提升为“core density response / polarizability”的源泛函。

---

## 1. 先把 core integration 写成带外源的形式

在积分掉 core 之前，引入一个标量外源 $\phi$，耦合到电子密度。对某个离子 $A$，core 看到的总势为

$$
\begin{aligned}
\Phi_A(\mathbf r,\tau)
&=
\phi_{\mathrm{ext}}(\mathbf r,\tau)
+
\int d\mathbf r'\,
v(\mathbf r-\mathbf r')\rho_v(\mathbf r',\tau),
\end{aligned}
$$

其中 $\rho_v$ 是 valence 密度，$v(\mathbf q)=4\pi/q^2$。对 core 自由度做 cumulant expansion，有

$$
\begin{aligned}
S_{\mathrm{core}}^{\mathrm{eff}}[\Phi]
&=
S_{\mathrm{core}}^{(0)}
+
\int n_c^0 \Phi
-
\frac{1}{2}
\int
\Phi\,\chi_c\,\Phi
+\cdots .
\end{aligned}
$$

线性项和静态部分给出普通 pseudopotential；二次项给出 core 的动态密度响应。你文中 Eq. (6) 正是这个动态二次核投影到 valence orbital 后的一粒子版本：静态 Hartree–Fock core self-energy 加上

$$
\begin{aligned}
(\delta V_{\mathrm{pp}})^{\mathrm{dyn}}_{st}(i\omega)
&=
\sum_c
\frac{M_{sc}M_{tc}^{\ast}}{i\omega+\Delta E_c},
\end{aligned}
$$

其中 pole 位于 core excitation energy $\Delta E_c$。在 momentum space 里，它被写成可分离 form

$$
\begin{aligned}
V_{\mathrm{dyn}}(K,K';i\omega)
&=
\sum_c
\frac{f^c_K f^c_{K'}}{i\omega+\Delta E_c},
\end{aligned}
$$

而 Eq. (7) 给出 $f^c_K$ 的 atomic core form factor。文中已经说明这个 dynamic term 是 conventional pseudopotential 没有的 core 动力学项。

---

## 2. 从 $f_K^c$ 反推出 core density response

写 core 的 transition density form factor 为

$$
\begin{aligned}
\tau_\lambda(\mathbf q)
&=
\langle 0|\hat{\rho}_c(\mathbf q)|\lambda\rangle ,
\end{aligned}
$$

其中 $\lambda$ 标记 core excitation channel，包括 shell、角动量和磁量子数。core density response 的谱表示是

$$
\begin{aligned}
\chi_c^A(\mathbf q,\mathbf q';i\xi)
&=
\sum_\lambda
\frac{
2\Delta_\lambda\,
\tau_\lambda^A(\mathbf q)
\tau_\lambda^{A\ast}(\mathbf q')
}
{\xi^2+\Delta_\lambda^2}
e^{-i(\mathbf q-\mathbf q')\cdot \mathbf R_A}.
\end{aligned}
$$

valence density 通过 Coulomb 与 core transition density 耦合，因此

$$
\begin{aligned}
f_\lambda(\mathbf q)
&=
v(\mathbf q)\tau_\lambda(\mathbf q),
\qquad
v(\mathbf q)=\frac{4\pi}{q^2}.
\end{aligned}
$$

也就是说

$$
\boxed{
\begin{aligned}
\tau_\lambda(\mathbf q)
&=
\frac{q^2}{4\pi}f_\lambda(\mathbf q)
\end{aligned}
}
$$

在这个约定下，Eq. (6) 的一边 pole

$$
\begin{aligned}
\frac{1}{i\omega+\Delta_\lambda}
\end{aligned}
$$

在响应函数里要换成偶的 Matsubara 组合

$$
\begin{aligned}
\frac{1}{i\xi+\Delta_\lambda}
+
\frac{1}{-i\xi+\Delta_\lambda}
&=
\frac{2\Delta_\lambda}{\xi^2+\Delta_\lambda^2}.
\end{aligned}
$$

于是 core-induced density kernel 是

$$
\begin{aligned}
v\chi_c v
&=
\sum_\lambda
\frac{
2\Delta_\lambda\,
f_\lambda(\mathbf q)f_\lambda^{\ast}(\mathbf q')
}
{\xi^2+\Delta_\lambda^2}.
\end{aligned}
$$

这里有一个重要细节：文中 Eq. (7) 主要是为 bandwidth narrowing 服务的 $s$-like scalar channel。它足以给 $z^{\mathrm{core}}$，但**单独的球对称 $\ell=0$ channel 不产生偶极 polarizability**。vdW 的 leading term 是 dipole–dipole fluctuation，所以需要把 Eq. (7) 推广到 $\ell=1$ 的 dipole-resolved form factor $f_{\lambda,1m}(\mathbf q)$。现有 Eq. (7) 是这个 multipole construction 的 $\ell=0$ 特例。

---

## 3. 从 generalized $f_\lambda(\mathbf q)$ 得到 $\alpha(i\xi)$

core transition dipole 为

$$
\begin{aligned}
d_{\lambda,i}
&=
\int d\mathbf r\, r_i\,\rho_\lambda(\mathbf r)
\\
&=
i
\left.
\frac{\partial \tau_\lambda(\mathbf q)}
{\partial q_i}
\right|_{\mathbf q=0}.
\end{aligned}
$$

代入 $\tau_\lambda=q^2f_\lambda/4\pi$，得到直接由 EFT form factor 给出的公式：

$$
\boxed{
\begin{aligned}
d_{\lambda,i}
&=
i
\left.
\frac{\partial}{\partial q_i}
\left[
\frac{q^2}{4\pi}
f_\lambda(\mathbf q)
\right]
\right|_{\mathbf q=0}
\end{aligned}
}
$$

因此 core 的动态偶极极化率张量为

$$
\boxed{
\begin{aligned}
\alpha^A_{ij}(i\xi)
&=
\sum_{\lambda\in \mathrm{core}}
\frac{
2\Delta_\lambda\,
d_{\lambda,i}^A
d_{\lambda,j}^{A\ast}
}
{\Delta_\lambda^2+\xi^2}
\end{aligned}
}
$$

各向同性原子或闭壳层情形下，

$$
\boxed{
\begin{aligned}
\alpha_A(i\xi)
&=
\frac{1}{3}\mathrm{Tr}\,\alpha^A_{ij}(i\xi)
\\
&=
\frac{2}{3}
\sum_\lambda
\frac{
\Delta_\lambda |\mathbf d_\lambda^A|^2
}
{\Delta_\lambda^2+\xi^2}
\end{aligned}
}
$$

这就是从 Eq. (6)/(7) 到 vdW 所需要的核心桥梁：

$$
\begin{aligned}
f_\lambda(\mathbf q)
&\Rightarrow
\tau_\lambda(\mathbf q)
\Rightarrow
d_\lambda
\Rightarrow
\alpha(i\xi).
\end{aligned}
$$

---

## 4. 接上 valence Coulomb/RPA screening

接下来不能直接用裸 Coulomb $v$，因为固体里 core fluctuation 会被 valence electrons 屏蔽。定义 valence irreducible response

$$
\begin{aligned}
\chi_v^{\mathrm{irr}}(i\xi)
\end{aligned}
$$

由 EFT tree-level KS propagator 和 matched charge vertices 计算。文中 Appendix E 已经把 valence interaction 写成 core-dressed Coulomb vertex，并说明在长程 Coulomb 仍然存在的情况下，金属里的 infrared-sensitive expansion 可以通过 ring diagrams 重组为 screened interaction $W^{\mathrm{val}}$。

因此定义

$$
\boxed{
\begin{aligned}
W_v(i\xi)
&=
\left[
v^{-1}-\chi_v^{\mathrm{irr}}(i\xi)
\right]^{-1}
\end{aligned}
}
$$

这就是被 valence electrons 屏蔽后的 Coulomb interaction。然后 core–core dispersion 的 RPA/log 形式可以写成

$$
\boxed{
\begin{aligned}
E_{\mathrm{vdW}}^{\mathrm{core}}
&=
\frac{1}{2\pi}
\int_0^\infty d\xi\,
\mathrm{Tr}_{\mathrm{inter}}
\left\{
\ln\left[
1-W_v(i\xi)\chi_c^\perp(i\xi)
\right]
+
W_v(i\xi)\chi_c^\perp(i\xi)
\right\}
\end{aligned}
}
$$

其中 $\mathrm{Tr}_{\mathrm{inter}}$ 表示只保留连接不同离子中心的非局域项，$\chi_c^\perp$ 是扣除了本地匹配部分后的 core response。这个式子展开到二阶就是 screened London interaction：

$$
\begin{aligned}
E_{\mathrm{vdW}}^{(2)}
&=
-\frac{1}{4\pi}
\int_0^\infty d\xi\,
\mathrm{Tr}_{\mathrm{inter}}
\left[
W_v\chi_c^\perp W_v\chi_c^\perp
\right].
\end{aligned}
$$

对两个中心 $A,B$，写成偶极形式：

$$
\boxed{
\begin{aligned}
E_{AB}^{(2)}
&=
-\frac{1}{2\pi}
\int_0^\infty d\xi\,
\mathrm{Tr}
\left[
\alpha_A(i\xi)
T_{AB}^{\mathrm{scr}}(i\xi)
\alpha_B(i\xi)
T_{BA}^{\mathrm{scr}}(i\xi)
\right]
\end{aligned}
}
$$

其中

$$
\begin{aligned}
T_{AB,ij}^{\mathrm{scr}}(i\xi)
&=
\nabla_{R_A,i}\nabla_{R_B,j}
W_v(\mathbf R_A,\mathbf R_B;i\xi).
\end{aligned}
$$

在真空、各向同性、远距离极限 $W_v\to v$，恢复

$$
\begin{aligned}
E_{AB}
&=
-\frac{C_6^{AB}}{R_{AB}^6},
\qquad
C_6^{AB}
=
\frac{3}{\pi}
\int_0^\infty d\xi\,
\alpha_A(i\xi)\alpha_B(i\xi).
\end{aligned}
$$

---

## 5. Double counting 的处理

这里的 double counting 要分三层处理。

第一，**不能把 $\Sigma_{\mathrm{dyn}}(0)$ 再加进 KS Hamiltonian**。文中 Appendix D 已经指出，$\Sigma_{\mathrm{dyn}}$ 含有一个静态 piece，这部分已经被 $H_{\mathrm{KS}}$ 中的 $V_{\mathrm{xc}}$ 和 pseudopotential 吸收；如果再放进 Dyson equation 会 double count。
所以一粒子修正必须用 subtracted kernel：

$$
\begin{aligned}
\Sigma_{\mathrm{dyn}}^{\mathrm{sub}}(i\omega)
&=
\Sigma_{\mathrm{dyn}}(i\omega)
-
\Sigma_{\mathrm{dyn}}(i\omega_{\mathrm{ref}}),
\end{aligned}
$$

其中 $\omega_{\mathrm{ref}}$ 取 $\epsilon_F$ 或 $0$，取决于使用的 Matsubara convention。

第二，**vdW energy 只取非局域 inter-center fluctuation**。本征原子 self-polarization、同一离子内的 $A=A$ block、以及线性 Hartree-like term 都应去掉：

$$
\begin{aligned}
\mathrm{Tr}_{\mathrm{inter}}
&=
\mathrm{Tr}
-
\sum_A \mathrm{Tr}_{AA}.
\end{aligned}
$$

这一步保留的是标准 semilocal DFT 缺失的长程 correlation，而不是重新拟合 atomic core energy 或 static pseudopotential。

第三，**如果某个 semicore shell 已经显式放进 valence，就必须从 $\chi_c$ 里删掉这个 shell**。文中在讨论 small-core pseudopotential 时已经采用同样原则：semicore states 显式进入 valence 后，EFT correction 只作用在 conduction band，以避免 double counting。

用 projector 语言，可以把最终修正写成

$$
\boxed{
\begin{aligned}
E_{\mathrm{vdW}}^{\mathrm{EFT}}
&=
\frac{1}{2\pi}
\int_0^\infty d\xi\,
\mathrm{Tr}_{\mathrm{inter}}
\left\{
\ln\left[
1-W_v(i\xi)
(1-P_{\mathrm{loc}})\chi_c(i\xi)(1-P_{\mathrm{loc}})
\right]
+
W_v(i\xi)
(1-P_{\mathrm{loc}})\chi_c(i\xi)(1-P_{\mathrm{loc}})
\right\}
\end{aligned}
}
$$

其中 $P_{\mathrm{loc}}$ 是与文中 $P_2^\Gamma$ 同源的本地匹配 projector：它去掉已经由 $V_H+V_{\mathrm{xc}}+V_{\mathrm{PSP}}$ 吸收的 relevant/marginal local pieces。文中 RPT 正是用 $\delta\Sigma=(1-P_2^\Gamma)\Sigma$ 来定义“KS tree level 已经包含什么、剩下什么作为修正”。

---

## 6. 实际计算流程

最小实现版是：

1. 对每个元素求 atomic core orbitals、$\Delta_\lambda$、以及 generalized multipole form factor $f_{\lambda,\ell m}(\mathbf q)$。
2. 用

   $$
   \begin{aligned}
   d_{\lambda,i}
   &=
   i\partial_{q_i}
   \left[
   q^2 f_\lambda(\mathbf q)/4\pi
   \right]_{\mathbf q=0}
   \end{aligned}
   $$

   得到 $\alpha_A(i\xi)$。
3. 用 KS/EFT valence response 构造

   $$
   \begin{aligned}
   W_v
   &=
   (v^{-1}-\chi_v^{\mathrm{irr}})^{-1}.
   \end{aligned}
   $$

4. 用

   $$
   \begin{aligned}
   E_{AB}^{(2)}
   &=
   -\frac{1}{2\pi}
   \int d\xi\,
   \mathrm{Tr}
   \left[
   \alpha_A T_{AB}^{\mathrm{scr}}
   \alpha_B T_{BA}^{\mathrm{scr}}
   \right]
   \end{aligned}
   $$

   做 pairwise screened vdW；或者用 log formula 做 full RPA many-body dispersion。
5. 删除 $A=A$ self terms、本地 $P_{\mathrm{loc}}$ pieces、以及所有已显式进入 valence 的 semicore channels。

一句话概括：**Eq. (6)/(7) 给的是 core 动力学 Wilson coefficient；vdW 需要把同一个 coefficient 解释为 core density response 的 Coulomb-dressed transition form factor，再取 long-wavelength dipole limit 得到 $\alpha(i\xi)$，最后用 valence-screened $W_v$ 而不是裸 $v$ 做 RPA correlation。**