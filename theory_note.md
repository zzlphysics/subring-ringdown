---
title: "Photon-Ring Subrings as Discrete Ringdown Samplers"
date: 2026-06-11
updated: 2026-06-11
tags: [gr-qc, astro-ph, black-hole, ringdown, photon-ring, polarization, QNM]
status: active-theory
projects: [subring-ringdown]
based_on:
  - "2605.11499 (Huang+2026) - Ringdown in photon polarization swings"
  - "1910.12873 (Gralla+Lupsasca 2020) - Lensing by Kerr Black Holes"
  - "2408.10303 (Zhong+2025 PRL) - Dynamical Lensing Tomography"
  - "2010.03683 (Hadar+2021) - Photon Ring Autocorrelations"
---

# Photon-Ring Subrings as Discrete Ringdown Samplers

## 核心思想

Kerr ringdown 期间的引力波微扰在偏振光中留下 QNM-locked 偏振振荡 [Huang+2026]。
光子环的 subring 层级由普适临界指数 (γ, δ, τ) 控制 [Gralla+Lupsasca 2020]。
**subring 指标 n 天然构成一个由时空几何控制的离散 retarded-time 采样器。**

在 near-critical 极限下，相邻 subring 的复偏振扰动比值收敛到普适复常数。
**存在两种自然定义**，测不同的物理量：

$$
\boxed{
R^{\rm pix}_{n,q}(\varphi) \equiv \frac{\Delta\mathfrak{P}_{n+q}(\varphi, t_o)}{\Delta\mathfrak{P}_n(\varphi, t_o)}
\;\xrightarrow{n\gg 1}\;
e^{-q\gamma}\, e^{i q(\omega\tau - m\delta)}
}
$$

$$
\boxed{
R^{\rm align}_{n,q}(\varphi) \equiv \frac{\Delta\mathfrak{P}_{n+q}(\varphi + q\delta, t_o)}{\Delta\mathfrak{P}_n(\varphi, t_o)}
\;\xrightarrow{n\gg 1}\;
e^{-q\gamma}\, e^{i q \omega\tau}
}
$$

其中：
- $\Delta\mathfrak{P} = \Delta Q + i\Delta U$ 是 subring 横截面积分后的复线偏振扰动
- $\omega = \omega_R - i\omega_I$ 是 QNM 复频率
- $q = 1$：相邻 half-orbit subring；$q = 2$：same-family subring
- $m$ 是 QNM 方位角量子数

**物理含义**：
- $R^{\rm align}$ 通过 pattern-aligned 比较消去了几何旋转 $m\delta$，纯净测量 $\omega_R\tau$（GW 相位）
- $R^{\rm pix}$ 在固定屏幕角比较，测量 $\omega_R\tau - m\delta$（GW+几何相位）
- 振幅 $|R| = e^{-q(\gamma - \omega_I\tau)}$ 在理想 flux-integrated 极限下编码 $\omega_I$

## 与现有文献的本质区别

### vs Huang+2026 (2605.11499)
Huang+2026 用时间域自相关提取 QNM。我们用 subring 指标 n 替代时间滞后 Δt——**空间维度替代时间维度**。

### vs Zhong+2025 PRL (2408.10303, Dynamical Lensing Tomography)
| | DLT | 本文 |
|---|---|---|
| 观测量 | 偏折角 | 复偏振ΔP |
| 采样 | 连续时间序列 | 离散subring指标n |
| 几何 | Geodesic deviation | Photon-ring monodromy |
| 需时间序列? | 是 | 否（单张snapshot） |

### vs Hadar+2021 (2010.03683)
Hadar+2021 做强度自相关，需要时间监测。我们用复偏振比值 R_n，n 替代时间滞后。

## Convention 约定

**为什么有两种 ratio？**

$\Delta\mathfrak{P}_n(\varphi)$ 中角向因子是 $e^{i m (\varphi - n\delta)}$。若分子分母取同一 $\varphi$，$-n\delta$ 和 $-(n+q)\delta$ 的差产生 $e^{-i q m\delta}$。若分子取 $\varphi + q\delta$ 对齐 subring 图案，则角向因子完全相同，$m\delta$ 被吸收。

两者测量不同的量，都是有效的观测量。

**为什么用 $\mathfrak{P}$ 而非 $\mathcal{P}$？**

$e^{-\gamma}$ 是 flux/Jacobian 因子，不是 local surface brightness 因子。
$\mathfrak{P}_n = \int_{\rm subring\;n} d\rho\, W_n(\rho,\varphi)\,[Q+iU](\rho,\varphi)$，
即 subring 横截面积分后的复偏振通量密度。这使得 $|\mathfrak{P}_{n+1}/\mathfrak{P}_n| \to e^{-\gamma}$ 
在几何光学极限下自然成立。

**q 参数**：相邻 half-orbit subring 步长为 (γ,δ,τ)；same-family (same-parity)
步长为 (2γ,2δ,2τ)。引入 q 统一处理。

## 数值结果 (a=0.7, m=2)

| 量 | q=1 | q=2 |
|---|---|---|
| $\arg R^{\rm pix}$ | 325.8° | 291.6° |
| $\arg R^{\rm align}$ | 73.2° | 146.4° |
| $|R|$ | 0.0559 | 0.0031 |

其中 $\omega_R\tau = 73.2°$, $m\delta = 107.4°$ (q=1).

## 验证路线

1. ✅ 解析推导 + 参数扫描
2. ⬜ Kerr geodesic subring classifier (Layer 0)
3. ⬜ Synthetic QNM screen (Layer 1) — 最小可行验证
4. ⬜ Full Teukolsky-CCK polarized ray tracing (Layer 2)

## 论文策略

- **主打相位**：$\arg R^{\rm align} \to q\omega_R\tau$ 最干净、最robust
- **振幅放副结果**：$|R|$ 依赖 flux/Jacobian convention
- **核心卖点**：subring index n 替代时间序列——"single-shot ringdown spectroscopy"
- **与 DLT 的 novelty 对比表**必须在正文中明确

## 项目文件结构

```
subring-ringdown/
  theory_note.md         ← 本文件
  conventions.md         ← 详细 convention 说明
  src/
    subring_core.py      ← 核心计算模块（修正后）
    fig_main.py          ← 主图
  archive/               ← 旧版本（已废弃）
  figures/
    subring_main.png     ← 四面板主图
```

## 参考文献

- [Huang+2026] arXiv:2605.11499
- [Gralla+Lupsasca 2020] PRD 101, 044031 (1910.12873)
- [Zhong+2025] PRL 134, 211402 (2408.10303)
- [Hadar+2021] PRD 103, 104038 (2010.03683)
