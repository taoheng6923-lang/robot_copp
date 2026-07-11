# robot_copp：基于 copp 算法的机器人时间最优轨迹规划模块 —— 设计方案

> **定位**：在 `copp`（Convex-Objective Path Parameterization）时间最优求解内核之上，新增一层**机器人运动指令 + 路径 blending** 前端，形成完整的"指令 → 几何路径 → 时间最优轨迹"规划链路。
>
> **输入**：关节运动 / 直线 / 圆弧三类运动指令序列，相邻指令之间可做 blending 过渡。
> **输出**：满足关节/笛卡尔运动学与动力学约束、按 copp 算法（TOTP-SPLP + 解析插值 + 可选 HLAW）时间最优的轨迹 $q(t),\dot q(t),\ddot q(t),\dddot q(t)$。
>
> **本文档只描述设计方案，不包含可运行代码。**
>
> ### 参考来源
> | 文档 | 贡献 |
> |------|------|
> | [`copp/docs/paper_notes.md`](../0.other_lib_code/copp/copp/docs/paper_notes.md) | TOTP-SPLP 理论：$(a,b,c)$ 状态变换、无损离散化、解析插值、分段线性目标、序列线性化、HLAW 可行性 |
> | [`copp/docs/code_reading_guide.md`](../0.other_lib_code/copp/copp/docs/code_reading_guide.md) | copp 求解器家族 API 与代码↔数学映射（`topp2_ra`/`topp3_lp`/`copp3_socp` 等） |
> | [`opencn-matlab/robot6dof_topp_design.md`](../0.other_lib_code/opencn-matlab/opencn-matlab/robot6dof_topp_design.md) | 六轴机器人离散化 + 关节导数（IK/Jacobian 链式法则）+ G2 Hermite 平滑的工程实现 |
> | [`1.OptimalHermiteInterpolation`](../1.OptimalHermiteInterpolation/) | 最优五次 Hermite 过渡（Frenet 框架 + $\alpha_0,\alpha_1$ 最优化，最小跃度代价） |
>
> ### 两点设计基调（已确认）
> 1. **实现无关**：本文档在数据流、数学公式、模块职责层面描述，不绑定具体语言 API；Rust（`copp` crate 新增 robot 层）、Python（`copp/python` 扩展）、MATLAB（OpenCN）均可据此实现。
> 2. **笛卡尔空间 blending + IK**：相邻指令的最优 Hermite 过渡在**笛卡尔位姿空间**构造（保刀尖几何、zone 距离含义清晰），随后对整条笛卡尔路径做连续解 IK 降维到关节空间，交给 copp 求解。

---

## 目录

1. [定位与范围](#1-定位与范围)
2. [总体架构](#2-总体架构)
3. [指令层：三类运动指令 → 笛卡尔位姿路径段](#3-指令层三类运动指令--笛卡尔位姿路径段)
4. [最优 Hermite Blending 层（笛卡尔空间）](#4-最优-hermite-blending-层笛卡尔空间)
5. [降维层：笛卡尔路径 → 关节空间 q(s)](#5-降维层笛卡尔路径--关节空间-qs)
6. [约束摄入层：物理约束 → 路径域不等式](#6-约束摄入层物理约束--路径域不等式)
7. [copp 求解层（核心）](#7-copp-求解层核心)
8. [HLAW 在线窗口层（长指令序列）](#8-hlaw-在线窗口层长指令序列)
9. [轨迹合成与输出](#9-轨迹合成与输出)
10. [模块划分与数据结构](#10-模块划分与数据结构)
11. [关键公式汇总](#11-关键公式汇总)
12. [设计决策与开放问题](#12-设计决策与开放问题)

---

## 1. 定位与范围

### 1.1 要解决的问题

给定一串机器人运动指令（关节 / 直线 / 圆弧），相邻指令交接处允许 blending 平滑，求一条**时间最优**的关节轨迹，使全程严格满足关节速度/加速度/加加速度/力矩以及笛卡尔 TCP 速度等约束。

`copp` 内核只解决"给定几何路径 $q(s)$，求时间最优/凸目标最优的时间律 $s(t)$"这一层（对应论文 [`paper_notes.md`](../0.other_lib_code/copp/copp/docs/paper_notes.md) 的 §3–§7）。**路径本身如何从指令生成、如何 blending**，不在 copp 范围内——这正是 `robot_copp` 模块要补的两层：**指令层**与 **blending 层**。

### 1.2 与三条已有路线的关系

| 已有方案 | 求解内核 | blending | 本模块的取舍 |
|----------|----------|----------|--------------|
| `robot6dof_topp_design.md` | OpenCN 两阶段 LP（TOPPRA 风格，B 样条决策变量） | G2 Hermite（`calcTransition`） | **复用**其指令→笛卡尔→IK→关节导数管线与 G2 blending；**替换**求解内核为 copp（SPLP+解析插值+HLAW） |
| `copp/python/DESIGN.md` | copp（Python 复刻） | 无（指令直接拼接，段间 C¹） | **复用**指令层三类命令与 PathBuilder 思路；**新增**最优 Hermite blending 层 |
| copp（Rust 库） | copp（原生） | 无 | 作为**求解内核**直接对接（`topp2_ra`→`topp3_lp/socp`/`copp3_socp`） |

一句话：**robot_copp = 指令层（借 python/DESIGN）+ 最优 Hermite blending 层（借 OptimalHermite / robot6dof §2.5）+ IK 降维（借 robot6dof §3）+ copp 时间最优内核（借 paper_notes / code_reading_guide）**。

### 1.3 符号约定

沿用 copp 代码习惯，**路径参数记为 $s$**（论文 `paper_notes.md` 中记作 $u$），状态量与 [`code_reading_guide.md`](../0.other_lib_code/copp/copp/docs/code_reading_guide.md) §2 一致：

| 符号 | 含义 |
|------|------|
| $u\in[0,1]$ | 单条指令的局部参数 |
| $s\in[0,s_f]$ | 全局路径参数（装配后，近似弧长归一化） |
| $\mathbf r(s)=[\mathbf p(s);\boldsymbol\varphi(s)]$ | 笛卡尔位姿路径（位置 + 姿态），$6\times1$ |
| $\mathbf q(s)=\mathrm{IK}(\mathbf r(s))$ | 关节路径，$n\times1$（六轴 $n=6$） |
| $(\cdot)'=d/ds$，$\dot{(\cdot)}=d/dt$ | 路径导数 / 时间导数 |
| $a(s)=\dot s^2,\ b(s)=\ddot s,\ c(s)=\dddot s\,\dot s$ | copp 状态量（对 $t$ 求导，见 §7.1） |
| $\mathbf J(\mathbf q)$ | 几何 Jacobian，$\dot{\mathbf r}=\mathbf J\dot{\mathbf q}$ |

---

## 2. 总体架构

```
┌──────────────────────────────────────────────────────────────────────────┐
│  输入：运动指令序列  [cmd₁, cmd₂, …, cmd_K]  +  机器人参数/约束             │
│        cmd ∈ {JointMove, LinearMove, CircularMove}，可带 blend 半径 zone   │
└───────────────────────────────┬────────────────────────────────────────────┘
                                │
        ═══════════════════════ 离线预处理（几何阶段） ═══════════════════════
                                │
        ┌───────────────────────▼───────────────────────┐
        │  §3 指令层：每条指令 → 笛卡尔位姿路径段 r_i(u)  │
        │       Joint: r=FK(q(u))；Line: 直线+SLERP；    │
        │       Arc: 圆弧几何+SLERP                       │
        └───────────────────────┬───────────────────────┘
                                │  位姿路径段数组（G0/G1 连续）
        ┌───────────────────────▼───────────────────────┐
        │  §4 最优 Hermite Blending（笛卡尔空间）         │
        │       对每对相邻段：zone 截短 → 最优五次 Hermite│
        │       G2 过渡（α0/α1 最小跃度）                 │
        └───────────────────────┬───────────────────────┘
                                │  C²/分段C³ 连续位姿路径 r(s)  ← 满足论文 Assumption 1
        ┌───────────────────────▼───────────────────────┐
        │  §5 降维：自适应离散化 + 连续解 IK + 链式法则   │
        │       → q(s), q'(s), q''(s), q'''(s)  于 {s_k} │
        └───────────────────────┬───────────────────────┘
                                │  PathDerivatives（对齐 Rust Path）
        ┌───────────────────────▼───────────────────────┐
        │  §6 约束摄入：物理约束 → 路径域 1/2/3 阶不等式  │
        │       Robot::with_axial_* / 站点索引存储        │
        └───────────────────────┬───────────────────────┘
                                │
        ═══════════════════════ 求解阶段（时间最优，copp） ═════════════════════
                                │
        ┌───────────────────────▼───────────────────────┐
        │  §7 copp 内核（论文 TOTP-SPLP）                 │
        │   ① topp2_ra 求种子 a⁽⁰⁾（2 阶可达/上界）      │
        │   ② SPLP 迭代（Algorithm 2）：                  │
        │      分段线性目标(PLP,eq.27)+序列线性化(eq.32)  │
        │      落地：topp3_lp+自建 PLP 目标，eq.30 迭代   │
        │            (备选 topp3_socp 精确)；凸目标 copp3_socp │
        │   ③ 解析插值 s↔t（Prop.1/2，无损）            │
        └───────────────────────┬───────────────────────┘
                                │  （长指令序列时外层包 §8 HLAW 三窗口）
        ┌───────────────────────▼───────────────────────┐
        │  §9 轨迹合成：解析插值采样 → q(t),q̇,q̈,q⃛       │
        └────────────────────────────────────────────────┘
```

**离线 / 在线划分**：§3–§6 为几何预处理，一次性完成，无优化循环；§7 求解可整段离线（短指令序列），也可由 §8 HLAW 分窗在线流式求解（长序列 / 边规划边执行）。

**与论文 Assumption 1 的衔接**：copp 内核要求几何路径 **$C^2$、分段 $C^3$ 连续且正则**。§4 的最优 Hermite blending 正是为此服务——它把相邻指令交接处的 G0/G1 间断提升为 G2 连续，使装配后的 $\mathbf r(s)$（进而 $\mathbf q(s)$）满足 copp 的输入前提；未成功 blending 的角点退化为 G1，copp 会在该处自动降速（可行但非全局最优）。

---

## 3. 指令层：三类运动指令 → 笛卡尔位姿路径段

每条指令在局部参数 $u\in[0,1]$ 上采样，输出统一的**笛卡尔位姿路径段** $\mathbf r_i(u)=[\mathbf p_i(u);\boldsymbol\varphi_i(u)]$（借鉴 `python/DESIGN.md` §4 的 `PathSegment`，但降维目标改为笛卡尔位姿以支持 §4 的笛卡尔 blending）。

### 3.1 JointMove（关节运动指令）

- 语义：起止关节角 $\mathbf q_\text{start}\to\mathbf q_\text{end}$，关节空间插值。
- 关节路径 $\mathbf q_i(u)$：单段用**五次多项式**（两端 $\dot q,\ddot q$ 可指定，默认 stationary 即两端 $\dot q=\ddot q=0$），复用 copp `path/spline.rs` 的 Hermite 样条（阶次 $p=5$）。
- 为进入笛卡尔 blending 层，用正运动学给出其位姿表示 $\mathbf r_i(u)=\mathrm{FK}(\mathbf q_i(u))$（仅在过渡区需要，见 §4.5 的快路径说明）。

### 3.2 LinearMove（直线指令）

- 位置：$\mathbf p_i(u)=(1-u)\mathbf p_\text{start}+u\,\mathbf p_\text{end}$。
- 姿态：四元数 SLERP，$\mathbf R_i(u)=\mathrm{Slerp}(\mathbf R_\text{start},\mathbf R_\text{end},u)$。
- 位置/姿态对 $u$ 的一~三阶导均有解析表达（直线导数为常量，SLERP 导数解析），无需数值差分。

### 3.3 CircularMove（圆弧指令）

- 几何：三点定圆（起点/终点/途经点）或显式 $(\text{center},\text{normal})$，求圆心 $\mathbf c$、半径 $r$、法向 $\mathbf n$，平面内正交基 $(\mathbf e_1,\mathbf e_2)$：
  $\mathbf p_i(\theta)=\mathbf c+r(\cos\theta\,\mathbf e_1+\sin\theta\,\mathbf e_2)$，$\theta(u)=\theta_0+u(\theta_1-\theta_0)$。
- 姿态：SLERP，同直线。
- 圆弧位置对 $u$ 的三阶导解析（三角函数），天然适配 copp `Jet3` 三阶前向自动微分。
- 退化情形（起点=终点、共线三点）显式报错。

### 3.4 输出数据结构 `PoseSegment`

```
PoseSegment_i
  .pose(u)      : r_i(u) = [p_i(u); φ_i(u)]，位姿采样/解析求值器
  .dpose_du     : 一~三阶导（供 §4 blending 端点 Frenet 框架、§5 自适应离散化）
  .length_hint  : 段弧长/角度估计（用于全局参数 s 分配与相对权重）
  .native_space : {"joint","cartesian"}（决定 §4.5 过渡桥接策略）
  .kin_branch   : 该段 IK 解分支/构型标记（保证 §5 连续解选择）
  .constraint_overrides : 该段专属限速/限力矩（可选，映射到论文分段约束 Assumption 2）
```

---

## 4. 最优 Hermite Blending 层（笛卡尔空间）

### 4.1 动机：G2 间断的危害

输入路径由多段基元拼接，交接处往往只有 G0（位置）或 G1（切向）连续。G2（曲率）间断会导致 $\mathbf r''(s_j^-)\ne\mathbf r''(s_j^+)$，经链式法则（§5）传到 $\mathbf q''$ 也跳变，copp 被迫在角点强制降速，无法真正时间最优（详见 `robot6dof_topp_design.md` §2.5）。blending 层把角点替换为一段 G2 连续的过渡曲线。

### 4.2 zone 距离（CutOff）与截短

对应工业机器人的 zone 参数（ABB z5/z10、KUKA C_DIS）。交接点两侧各截去弧长 $\Delta_j$：

$$
\Delta_j=\min\!\left(\Delta_\text{max},\ \tfrac{L_1}{3},\ \tfrac{L_2}{3}\right)
$$

其中 $\Delta_\text{max}$ 为用户设定最大 zone，$L_1,L_2$ 为相邻两段弧长，$L/3$ 上限防止相邻过渡区重叠。截短后在两个新端点间插入过渡曲线。典型取值 $1\sim10\,\text{mm}$（精加工取小、搬运取大）。

### 4.3 最优五次 Hermite 过渡（G2，最小跃度）

在截短端点的 **Frenet 框架**（切向 $\mathbf t$、法向 $\mathbf n$、曲率 $\kappa$）下，构造满足两端 位置+切向+曲率（G2）的五次多项式：

$$
\mathbf p_5(u)=\mathbf r_L h_{00}+\alpha_0\mathbf t_0 h_{10}+(\beta_0\mathbf t_0+\alpha_0^2\kappa_0\mathbf n_0)h_{20}
+\mathbf r_R h_{01}+\alpha_1\mathbf t_1 h_{11}+(\beta_1\mathbf t_1+\alpha_1^2\kappa_1\mathbf n_1)h_{21}
$$

自由参数为切向拉伸 $\alpha_0,\alpha_1>0$（$\beta_0,\beta_1$ 由其线性确定）。按两端曲率是否为零分四种情形求解（直线-直线 → 2×2 线性；直线-曲线/曲线-直线 → 三次多项式；曲线-曲线 → 结式消元至九次），多根时按**最小化三阶导范数积分（跃度代价 $\int\|\mathbf p_5'''\|^2$）**选最优解。此即 [`1.OptimalHermiteInterpolation`](../1.OptimalHermiteInterpolation/) 与 `robot6dof_topp_design.md` §2.5.2 的"最优 Hermite 过渡"，直接复用其已有实现（`G2_Hermite_Interpolation_nAxis`）。

> **为何是"最优"而非任意 G2 过渡**：$\alpha_0,\alpha_1$ 的自由度让过渡曲线在满足 G2 的前提下最小化跃度积分——这与 copp 后续的 jerk 时间最优目标同向，过渡段本身就"好走"，减少 copp 在过渡区的降速。

### 4.4 位姿联合 blending

`G2_Hermite_Interpolation_nAxis` 天然支持 N 维。6D 位姿 $\mathbf r=[p_x;p_y;p_z;\varphi_x;\varphi_y;\varphi_z]$ 中位置轴（mm）与姿态轴（rad）**同时**纳入 Frenet 框架，共用同一 zone 距离 $\Delta_j$；用特征长度 $D$ 做量纲归一化，防止 mm 与 rad 混合的数值病态。姿态过渡在旋转矢量/四元数切空间上处理以保证流形一致性。

### 4.5 混合指令的过渡桥接

三类指令的过渡在**笛卡尔空间**统一处理，但需注意 JointMove 的 native space 是关节空间：

| 相邻组合 | 过渡构造 |
|----------|----------|
| Cartesian–Cartesian（Line/Arc 之间） | 直接在笛卡尔位姿空间截短 + 最优 Hermite（主路径，最干净） |
| Joint–Cartesian / Cartesian–Joint | 对 Joint 段截短端点用 FK 求其位姿 + 位姿导数（$\dot{\mathbf r}=\mathbf J\dot{\mathbf q}$ 等），再与 Cartesian 端点在笛卡尔空间做 Hermite 过渡；IK 回关节时以 Joint 段端点已知 $\mathbf q$ 为 seed 保证分支连续 |
| Joint–Joint | 若两段均纯关节且无笛卡尔精度要求，可选**关节空间**直接做最优 Hermite（$n$ 维），跳过 FK/IK 往返（快路径）；否则同上走笛卡尔 |

> **快路径**：纯 JointMove 段的**内部**（非过渡区）$\mathbf q(u)$ 已解析已知，无需 IK；只有引入了笛卡尔过渡曲线的 **zone 区间**才需要 §5 的 IK 采样。降维层据 `native_space` 标记跳过不必要的 IK。

### 4.6 输出

装配后得到一条 **$C^2$、分段 $C^3$ 连续**的笛卡尔位姿路径 $\mathbf r(s)$（原始截短段 + TransP5 过渡段交替），并记录每段的全局参数区间 $[s_i,s_{i+1}]$、过渡段高曲率标记（供 §5 自适应加密）、以及约束域边界（供 §6/§7 网格对齐论文 §3.1.2 要求）。

---

## 5. 降维层：笛卡尔路径 → 关节空间 q(s)

本层将 §4 的 $\mathbf r(s)$ 转成 copp 需要的关节路径导数，**完全复用** `robot6dof_topp_design.md` §二、§三的方法，此处只给要点与差异。

### 5.1 自适应离散化

曲率驱动步长（位置曲率 $\kappa_\text{pos}$ + 姿态曲率 $\kappa_\text{ori}$）：

$$
\Delta s_m=\min\!\left(\Delta s_\text{max},\ \sqrt{\tfrac{8\varepsilon_\text{pos}}{\kappa_\text{pos}(s_m)}},\ \sqrt{\tfrac{8\varepsilon_\text{ori}}{\kappa_\text{ori}(s_m)}}\right)
$$

过渡段曲率高会自动加密。**关键补充（对齐论文 §3.1.2）**：离散网格 $\{s_k\}$ **必须包含所有约束域边界**（不同工艺段/不同坐标系约束的分界 $s$），以保证每个区间内约束光滑，这是论文离散化误差 $O(\Delta^2)$ 有界（Theorem 1）的前提。

### 5.2 连续解 IK

逐点 $\mathbf q_m=\mathrm{IK}(\mathbf r(s_m))$，按"最近解 → 构型一致 → 关节限位合法"优先级选分支，相邻点以上一点解为 seed，避免关节跳变（`robot6dof` §3.1）。JointMove 段快路径直接用解析 $\mathbf q(u)$。

### 5.3 关节导数（Jacobian 链式法则，解析）

$$
\mathbf q'_m=\mathbf J_m^{-1}\mathbf r'_m,\quad
\mathbf q''_m=\mathbf J_m^{-1}\bigl(\mathbf r''_m-\mathbf J'_m\mathbf q'_m\bigr),\quad
\mathbf q'''_m=\mathbf J_m^{-1}\bigl(\mathbf r'''_m-2\,\mathbf J'_m\mathbf q''_m-\mathbf J''_m\mathbf q'_m\bigr)
$$

> 注意三阶式中 $\mathbf J'\mathbf q''$ 的系数是 **2**（对 $\mathbf J'\mathbf q'+\mathbf J\mathbf q''=\mathbf r''$
> 再求导，乘积法则产生两个 $\mathbf J'\mathbf q''$ 项）。本文档 v0.5 及之前版本漏写了该系数，
> 已由实现（`path/lowering/derivatives.py`）的有限差分交叉验证确认并修正。

$\mathbf J',\mathbf J''$ 沿路径方向前向/中心差分（`robot6dof` §3.2–3.4；实现取方向差分
$\mathbf J'=\mathrm DJ[\mathbf q']$、$\mathbf J''=\mathrm D^2J[\mathbf q',\mathbf q']+\mathrm DJ[\mathbf q'']$）。
精确 $\mathbf q'''$ 是保证 jerk 约束精度 $<2\%$ 的关键。

### 5.4 奇异处理

$\mathrm{rcond}(\mathbf J_m)<10^{-6}$ 时切阻尼最小二乘逆（DLS，$\lambda=0.05$），并局部加密采样。奇异区间打标记，供 §6 施加更保守约束或 §7 降速。

### 5.5 输出 `PathDerivatives`

与 Rust `Path::evaluate_up_to_3rd` 同结构：$\mathbf q,\mathbf q',\mathbf q'',\mathbf q'''$，形状 $(n,N)$，外加 $\{s_k\}$、段/约束域索引、奇异标记。直接可被 §6 的 `Robot` 摄入。

---

## 6. 约束摄入层：物理约束 → 路径域不等式

将用户物理约束按论文附录 A.2 映射为路径域 1/2/3 阶不等式（对应 copp `Robot::with_axial_velocity/acceleration/jerk`、`with_axial_torque` 及 `Constraints`）。

### 6.1 支持的约束类型

| 类型 | 阶数 | 路径域形式 |
|------|------|-----------|
| 轴向速度 $\dot{\mathbf q}=\mathbf q'\dot s$ | 1 | $\dot s\le v_\alpha(s)$，即 $a\le\bar a_\alpha$（线性上界） |
| 轴向加速度 $\ddot{\mathbf q}=\mathbf q''\dot s^2+\mathbf q'\ddot s$ | 2 | $\mathbf n_\beta a+\mathbf m_\beta b\le\mathbf g_\beta$（40b，逐轴线性） |
| 轴向 jerk $\dddot{\mathbf q}=\mathbf q'''\dot s^3+3\mathbf q''\dot s\ddot s+\mathbf q'\dddot s$ | 3 | $\mathbf r_\gamma a+\mathbf s_\gamma b+\mathbf t_\gamma c+\mathbf h_\gamma\le\mathbf f_\gamma a^{-1/2}$（40c/6f） |
| 关节力矩 / 力矩率 | 2/3 | 论文 eq.43–45；copp `RobotTorque` |
| TCP 位置速度模长 $\|\dot{\mathbf p}\|$ | 1 | $\|\mathbf p'(s)\|^2\,a\le v_\text{tcp,max}^2$（线性上界；$\mathbf p'=0$ 时平凡满足，退化安全） |
| TCP 姿态角速度模长 $\|\boldsymbol\omega\|$ | 1 | $\|\mathbf J_\omega(\mathbf q)\,\mathbf q'(s)\|^2\,a\le\omega_\text{tcp,max}^2$（线性上界） |
| 弦高误差 / 轮廓误差（工艺约束） | 1 | 论文 eq.46/47；转成 $a$ 上界 |

> **TCP 约束仅限速度层**：本版 TCP 只限制**位置速度的模** $\|\dot{\mathbf p}\|$ 与**姿态角速度的模** $\|\boldsymbol\omega\|$，暂不引入 TCP 加速度/jerk 约束。两者都是对 $a=\dot s^2$ 的**线性上界**（与轴向速度同构，逐点取 min 合并进 $\bar a(s)$），不产生二阶锥，LP 版求解器即可胜任。系数 $\|\mathbf p'(s)\|^2$、$\|\mathbf J_\omega\mathbf q'(s)\|^2$ 在 §5 降维阶段随 $\mathbf q',\mathbf q''$ 一并预计算。

> **为何不采用 CNC 的切向/法向分解**（$v_t,a_t,j_t$ / $a_n,j_n$，论文 eq.41/42）——这是机器人相对 CNC 的一个关键差异：
> - **关节侧不需要**：机器人电机限位是**逐轴**的（轴向约束已覆盖），关节路径上的"切向速度" $\|\mathbf q'\|\dot s$ 不对应任何物理限位。
> - **笛卡尔侧会退化**：切向/法向分解依赖刀尖位置恒动（$\mathbf p'\ne0$ 正则），CNC 满足；但机器人存在**纯姿态调整**（TCP 位置不动、仅换姿），此时 $\mathbf p'\to0$，切向/法向方向未定义、系数 $\dfrac{\mathbf p'^\top\mathbf p''}{\|\mathbf p'\|}$ 发散，无法形成良态约束。
> - **替代方案（本版采用）**：笛卡尔层只保留**速度模长**约束——TCP 位置速度模长与姿态角速度模长，二者都是 $a$ 的**线性上界**（无 $\|\mathbf p'\|$ 相除、无除零退化），对应真实工艺速度窗口；TCP 加速度/jerk 暂不约束（见下表说明）。
>
> 需要说明的是：在正则段（$\mathbf p'\ne0$）切向/法向加速度**确实**是 $(a,b)$ 的线性表达式（CNC 正据此处理），并非全局非线性；但上述"纯姿态退化"与"无对应物理限位"两点，使其不宜作为机器人约束。故本版从 §6.1 移除切向/法向行，Cartesian 约束一律走模长形式。

### 6.2 关键性质（继承论文能力）

- **非对称、参数变**约束：上下界可随 $s$ 变、可不对称（论文 Assumption 2）。
- **跨坐标系**：同一路径可在工件坐标系与刀具/关节坐标系**同时**施加约束（论文 §6.2.2 已证不同坐标系约束不等价）。
- **分段约束**：不同指令段/工艺段用不同约束索引 $\alpha,\beta,\gamma$，站点索引循环存储（copp `constraints.rs`），网格已在 §5.1 对齐约束域边界。

---

## 7. copp 求解层（核心）

本层的**算法方案完全采用论文的 TOTP-SPLP**（分段线性目标 + 序列线性化 + 无损离散化 + 解析插值），理论细节见 [`paper_notes.md`](../0.other_lib_code/copp/copp/docs/paper_notes.md) §4–§7。copp Rust 库作为**可选后端**：其无损离散化、jerk 线性化、解析插值可直接复用；但**目标函数与迭代需按 §7.2.1 补齐到论文方案**——现成的 `topp3_lp` 只是论文的 TOTP-LP 基线（$\max\int a$、单次求解），并非本方案。

### 7.1 状态变换与无损离散化

构造状态 $a(s)=\dot s^2,\ b(s)=\ddot s,\ c(s)=\dddot s\,\dot s$，动力学线性化 $a'=2b,\ b'=c$，目标 $J=\int_0^{s_f}\!ds/\sqrt{a}$。区间内按 Prop.1（非静止点 $c$ 零阶保持）/ Prop.2（静止点 $\dddot s$ 零阶保持）离散，配合解析插值实现**无损离散化**，区间违约 $O(\Delta^2)$ 有界（Theorem 1）。这些已由 copp 内核实现，本层只需提供 §5/§6 的 `PathDerivatives` 与约束。

### 7.2 TOTP-SPLP：分段线性目标 + 序列线性化（论文方案）

本层严格按论文 Algorithm 2 的 TOTP-SPLP 求解，而非 $\max\int a$ 的 LP 基线。

**① 分段线性目标（PLP）**：真实时间被积函数为 $1/\sqrt{a_k}$，用若干采样点 $0<\delta_{k,0}<\delta_{k,1}<\dots$ 的**割线上包络**逼近（论文 eq.27），引入辅助变量 $J_k$ 与线性约束（eq.29d）把 $\max_i(\cdot)$ 摊平为线性。相比 LP 基线 $\max\int a$，PLP 逼近真实时间惩罚，**近最优**；再加下界 $a_k\ge\delta_{k,0}>0$（论文 Prop.3，$\delta_{k,0}$ 足够小不影响最优），等效在 $a_k\to0^+$ 施加无穷惩罚，**根除零进给奇异**。

**② jerk 凹约束线性化**：3 阶凹约束 $a^{-1/2}$ 在 $a^{(p-1)}$ 处取切线（论文 eq.32，非保守，优于 pseudo-jerk）：
$$\Big(r+\tfrac{f}{2\,a^{(p-1)\,3/2}}\Big)a+s\,b+t\,c\le\tfrac{3f}{2\sqrt{a^{(p-1)}}}-h$$
线性可行域是原凹域的线性子集。此式即 copp `build_with_linearization` 已实现者（`formulation.rs:246`），可直接复用。

**③ 序列迭代（SPLP 的"S"）**：从种子 $a^{(0)}$ 出发，第 $p$ 次在 $a^{(p-1)}$ 处线性化 + 解一次子问题得 $a^{(p)}$，重复至停止准则（论文 eq.30）：
$$\big|t_f^{(p)}-t_f^{(p-1)}\big|<\varepsilon_t \quad\text{或}\quad \big\|a^{(p)}-a^{(p-1)}\big\|<\varepsilon_a \quad\text{或}\quad p\ge N_\text{iter}$$
实践中 $N_\text{iter}$ 取小（2~3）即得可行近优解（论文 §4.1）。

**④ 种子与可行性**：种子 $a^{(0)}$ 取 2 阶问题最优解（`topp2_ra`），**不必**是可行轨迹，只需保证 $p=1$ 的线性化子问题可行（论文 §7.3）。由此 **Theorem 2** 保证 SPLP 全程可行并收敛到（PLP 目标 + 线性化约束下的）KKT 解。

#### 7.2.1 与 copp Rust 库的落地映射（重要）

现成 `topp3_lp` **不等于**本方案：它用 LP 目标 $\max\int a\,ds$ 且**单次**求解（= 论文 TOTP-LP 基线，`topp3_lp.rs:24/241`）。要落地论文 TOTP-SPLP 有两条路线：

| 路线 | 做法 | 取舍 |
|------|------|------|
| **(a) 忠实 PLP + LP（推荐，算力最省）** | 扩展 `topp3_lp`：决策向量增补辅助变量 $J_k$，目标改为 $\min\sum(u_{k+1}-u_{k-1})J_k$（eq.27 割线上包络），约束矩阵加 eq.29d 行与下界 $a_k\ge\delta_{k,0}$；把"`build_with_linearization`(重线性化) + LP 求解"包成 Algorithm 2 的**迭代循环** | 即论文原版 TOTP-SPLP，纯 LP、**算力最省**；需自建 PLP 目标向量、$J_k$ 变量与 eq.29d 约束行 |
| (b) 精确目标 + SOCP（备选） | 直接用 `topp3_socp`：以 $\eta_k=1/\sqrt{a_k}$ + 二阶锥**精确**建模真实时间目标（`topp3_socp.rs:284`），外层同样按 eq.30 迭代重线性化 | 目标精确（非近似），实现最省（无需自建目标）；代价是 SOCP 比 LP 略贵、更适合追求数值鲁棒时 |

两条路线都复用同一套 jerk 线性化（eq.32）、无损离散化与解析插值；差别只在**目标**（PLP 近似走 LP vs 精确走 SOCP）。**默认走 (a)**：PLP+LP 是论文原版 TOTP-SPLP，单次迭代计算量最低（LP ≪ SOCP），最契合"在线、算力受限"的机器人场景；论文也表明 PLP 已足够接近精确目标的时间最优（`paper_notes.md` §7.1 / 附录 B.2 中 C4≈C5）。仅当需要精确目标或更强数值鲁棒时再用 (b)。凸目标（时间+热能等）走 `copp3_socp`。

> **单次求解 = 论文一次迭代**：无论 (a)/(b)，只调一次 `build_with_linearization + solve` 相当于论文 SPLP 的 $p=1$（即可行初值 $a^{(1)}$，对应 HLAW 的"可行窗"）。时间最优性来自后续迭代（对应"最优窗"，§8）。

### 7.3 求解器家族选择

| 需求 | copp 求解器 | 与论文方案的关系 |
|------|-------------|------------------|
| 种子 $a^{(0)}$ / 2 阶时间最优 | `topp2_ra` | 提供 SPLP 线性化种子（论文 §5.1 种子窗的 2 阶问题） |
| 可行性 / 可达域分析 | `reach_set2` | 后向/双向可达集 |
| 2 阶凸目标 | `copp2_socp` | — |
| **3 阶时间最优（本方案，默认）** | 迭代 `topp3_lp` + 自建 PLP 目标（推荐，算力最省）／ 迭代 `topp3_socp`（精确目标，备选） | 落地论文 TOTP-SPLP，见 §7.2.1 |
| 3 阶凸目标（时间+热能等） | `copp3_socp` | `CoppObjective::{Time,ThermalEnergy,TotalVariationTorque,Linear}` 线性组合 |

标准 3 阶调用序列：`topp2_ra` 得 $a^{(0)}$ → **循环** { `Topp3ProblemBuilder::new(&mut robot, idx_s_start, a^{(p-1)}, a_bnd, b_bnd).build_with_linearization()` → `topp3_lp`+PLP 目标（或备选 `topp3_socp`）得 $a^{(p)}$ } 直到 eq.30 → 解析插值。

> **目标决定求解器**：本版运动学约束全为逐轴/箱式线性，故求解器由**目标**决定——默认用论文 PLP 分段线性目标走 **LP（`topp3_lp` 循环 + 自建 PLP 目标向量，算力最省）**；仅当需要精确时间目标或更强数值鲁棒时改走 SOCP（`topp3_socp`）；凸目标（热能等）走 `copp3_socp`。TCP 两个速度模长约束只是 $a$ 的线性上界，不改变求解器选择。

### 7.4 解析插值 s↔t

用 copp 的 `s_to_t_topp3` / `t_to_s_topp3`（论文 Prop.1/2 的闭式 $\Phi_k,\Phi_k^{-1}$）把解 $(a,b)$ 转成时间律与轨迹采样，**无插值误差**，无需弧长参数化。2 阶方法用 `*_topp2` 对应版本。

### 7.5 边界条件

- 默认两端静止 $a=b=0$（工业"静止起、静止止"），也是 HLAW 可行性定理（Theorem 3）的充分条件。
- 非静止边界用于**多窗口衔接**（§8）：过渡段两端的 $a,b$ 由相邻窗口解提供，通过 `a_boundary,b_boundary` 传入 builder。

---

## 8. HLAW 在线窗口层（长指令序列）

短指令序列可整段离线求解（§7 一次）。长序列（网格点上万乃至百万，或边规划边执行）需分窗在线，并保证**每个子问题可行**——这是 copp 内核之上、robot_copp 需要新增的调度层，直接落地论文 §5 的 HLAW。

### 8.1 三层前瞻窗口

三个窗口沿路径层级式前移，始终 $k^r_\text{seed}<k^l_\text{fea}$、$k^r_\text{fea}<k^l_\text{opt}$：

| 窗口 | 求解 | copp 调用 |
|------|------|-----------|
| 种子窗 | 2 阶 LP 得 $a^{(0)}$（终止态置 0） | `topp2_ra` / 2 阶 LP |
| 可行窗 | 在 $a^{(0)}$ 线性化的 3 阶 LP（$p=1$）得可行初值 | `topp3_lp` 单次 |
| 最优窗 | SPLP 迭代（$p\ge2$）得近最优 | `topp3_lp/socp` 迭代 |

### 8.2 可行性保证（Theorem 3）

**核心机制**：可行窗与最优窗都在**同一种子窗生成的 $a^{(0)}$** 上线性化——线性化点**跨窗口一致**，配合每窗尾部可补的静止轨迹（Assumption 3），归纳证明所有窗口可行。实测长刀路 **100% 可行**（论文 §6.3.2，把已有框架上千次不可行降为 0）。

### 8.3 与滑动窗口 MPC 的区别

`robot6dof_topp_design.md` §5.5 的滑动窗口 MPC 属于**非分层**（LAW/PW 类）：逐窗依赖上一窗最优解，线性化点不一致 → 存在"人为不可行"。HLAW 用**层级单向依赖 + 逐层求解**消除该不一致，是本模块相对 robot6dof 方案的**关键升级**（用可证可行性换取工程可靠性）。

### 8.4 窗口参数

窗宽与重叠不影响可行性，但影响时间最优性与算力（论文 §5.4 复杂度 $O(K\beta_L N_\text{win}^3 N)$，每区间成本与总长无关）。重叠过短 → 连接处未消除减速；过长 → 冗余计算。调优准则：终止时间随窗宽/重叠增大不再显著下降即够（论文附录 A.3）。

---

## 9. 轨迹合成与输出

对 §7/§8 的解 $(a,b,c)$ 用 §7.4 解析插值，按伺服周期 $T_s$ 采样：

- $\mathbf q(t)=\mathbf q(s(t))$
- $\dot{\mathbf q}(t)=\mathbf q'(s)\dot s$，$\dot s=\sqrt{a}$
- $\ddot{\mathbf q}(t)=\mathbf q''(s)a+\mathbf q'(s)b$
- $\dddot{\mathbf q}(t)=\mathbf q'''(s)\,a^{3/2}+3\mathbf q''(s)\,\sqrt a\,b+\mathbf q'(s)\,c$（由 $c=\dddot s\dot s$ 反解）

输出 $\mathbf q(t),\dot{\mathbf q}(t),\ddot{\mathbf q}(t),\dddot{\mathbf q}(t)$，全程 $C^2$（分段 $C^3$）。**验证指标**（论文 §6.1.2）：超限率 $R_v$、超限时长比 $D_v$，应 $<0.1\%$；可在区间内多点评估约束（论文 eq.26）抑制进给率振荡。

---

## 10. 模块划分与数据结构

> **与实际实现的差异**：下面这棵树是"实现无关"的职责划分（本文档 §12.1 的立场），
> 实际 Python 实现把它拆到了顶层 `robot/`（本节没有单列的机器人本体，运动学/动力学
> 现由顶层独立包 `robot/ur5.py` 提供）+ `trajectory-planning/{copp,path,planner}`
> 三个包：`commands/blending/lowering` → `trajectory-planning/path/`；
> `constraints/solve` → `trajectory-planning/copp/`；`hlaw/synth` 及本文档未单列的
> 门面 → `trajectory-planning/planner/`。当前只有 `copp/`（对应本节 §6/§7）与顶层
> `robot/` 已实现，其余仍是占位目录。逐模块的真实路径、实现状态见
> [`README_M1.md`](./README_M1.md)；本节的职责划分与算法归属依旧是设计权威。

```
robot_copp/
│
├── commands/                     ── §3 指令层
│   ├── joint_move                    JointMove → q(u) 五次多项式 → r=FK
│   ├── linear_move                   LinearMove → 直线 + SLERP
│   ├── circular_move                 CircularMove → 圆弧几何 + SLERP
│   └── pose_segment                  PoseSegment 数据结构
│
├── blending/                     ── §4 最优 Hermite blending（笛卡尔）
│   ├── zone_trim                     zone 距离 Δ_j + 截短
│   ├── optimal_hermite               [复用] G2_Hermite（Frenet + α0/α1 最优）
│   ├── pose_blend                    位姿联合过渡 + 量纲归一化 D
│   └── junction_bridge               §4.5 混合指令桥接 / 快路径判定
│
├── lowering/                     ── §5 降维到关节空间（复用 robot6dof §2-3）
│   ├── adaptive_sample               曲率驱动步长 + 约束域边界对齐
│   ├── continuous_ik                 连续解 IK
│   ├── joint_derivatives             q',q'',q''' 链式法则
│   └── singularity                   DLS + 局部加密
│
├── constraints/                  ── §6 约束摄入（对接 copp Robot/Constraints）
│   └── robot_constraint_ingest       物理约束 → 路径域 1/2/3 阶不等式
│
├── solve/                        ── §7 copp 内核对接（论文 TOTP-SPLP）
│   ├── seed                          topp2_ra 种子 a⁽⁰⁾
│   ├── splp_loop                     SPLP 迭代(Algorithm 2)：重线性化+topp3_lp+自建PLP目标(默认) / topp3_socp(备选)；凸目标 copp3_socp
│   ├── plp_objective                 PLP 目标向量(eq.27)+辅助变量 J_k+eq.29d 约束行+下界 a_k≥δ_0
│   └── analytic_interp               s_to_t / t_to_s（Prop.1/2）
│
├── hlaw/                         ── §8 分层前瞻窗口（长序列，可选）
│   ├── window_schedule               种子/可行/最优 三窗调度
│   └── boundary_relay                跨窗口边界条件传递
│
└── synth/                        ── §9 轨迹合成与验证
    ├── resample                      解析插值采样 → q(t),q̇,q̈,q⃛
    └── verify                        R_v / D_v 超限指标
```

**核心数据结构**（实现无关的字段视图）：

```
PoseSegment       : pose(u), dpose_du(1~3), length_hint, native_space, kin_branch, overrides
BlendedPath       : r(s)（截短段 + TransP5 段），段区间 [s_i,s_{i+1}]，约束域边界，过渡高曲率标记
PathDerivatives   : s_grid[N], q[n×N], dq[n×N], ddq[n×N], dddq[n×N], singular[N]（对齐 Rust Path）
CoppProblem       : PathDerivatives + Constraints + 边界(a,b) + 目标(CoppObjective) + 求解器选择
TrajectoryResult  : s_grid, a_profile, b_profile, num_stationary, t_final, t_s, q/qd/qdd/qddd(t)
```

---

## 11. 关键公式汇总

**zone 距离**：$\Delta_j=\min(\Delta_\text{max},L_1/3,L_2/3)$

**最优五次 Hermite（G2）**：$\mathbf p_5(u)$ 见 §4.3，$\alpha_0,\alpha_1$ 由 $\min\int\|\mathbf p_5'''\|^2$ 选取

**关节导数（链式法则）**：$\mathbf q'=\mathbf J^{-1}\mathbf r'$，$\mathbf q''=\mathbf J^{-1}(\mathbf r''-\mathbf J'\mathbf q')$，$\mathbf q'''=\mathbf J^{-1}(\mathbf r'''-2\mathbf J'\mathbf q''-\mathbf J''\mathbf q')$（三阶式系数 2，见 §5.3 注）

**copp 状态与动力学**：$a=\dot s^2,\ b=\ddot s,\ c=\dddot s\dot s$；$a'=2b,\ b'=c$；$J=\int ds/\sqrt a$

**分段线性目标（SPLP）**：$J_k=\max_i\dfrac{\delta_{k,i-1}+\sqrt{\delta_{k,i-1}\delta_{k,i}}+\delta_{k,i}-a_k}{(\sqrt{\delta_{k,i-1}}+\sqrt{\delta_{k,i}})\sqrt{\delta_{k,i-1}\delta_{k,i}}}$

**3 阶凹约束线性化**：$a_k^{-1/2}\ge\dfrac{3a^{(p-1)}_k-a_k}{2(a^{(p-1)}_k)^{3/2}}$

**离散化误差界（Theorem 1）**：约束违约 $\le C\Delta^2=O(\Delta^2)$

**轨迹合成**：$\dot{\mathbf q}=\mathbf q'\sqrt a$，$\ddot{\mathbf q}=\mathbf q''a+\mathbf q'b$，$\dddot{\mathbf q}=\mathbf q'''a^{3/2}+3\mathbf q''\sqrt a\,b+\mathbf q'c$

---

## 12. 设计决策与开放问题

### 12.1 已定设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 文档载体 | 实现无关（算法/架构层） | Rust/Python/MATLAB 均可据此实现；接口示例仅作参照 |
| blending 空间 | 笛卡尔位姿空间 + IK 降维 | 保刀尖几何、zone 含义清晰，与 robot6dof 管线一致；JointMove 经 FK 桥接或走关节快路径 |
| 过渡曲线 | 最优五次 Hermite（G2，最小跃度） | 复用现成实现；与 copp jerk 时间最优目标同向 |
| 求解内核 | 论文 TOTP-SPLP（PLP 分段线性目标 + 序列线性化 + 解析插值） | 近最优、避奇异、约束严格、误差 $O(\Delta^2)$ 有界；落地默认 `topp3_lp`+自建 PLP 目标（LP，算力最省），精确目标备选 `topp3_socp`（§7.2.1） |
| 长序列可行性 | HLAW 三窗（替代非分层 MPC） | 线性化点跨窗一致 → 可证可行，实测 100% |
| 默认阶数 | 3 阶（jerk） | 抑制振动、提升表面质量；2 阶作种子 |

### 12.2 开放问题 / 后续扩展点

1. **HLAW 尚未在 copp Rust 库内实现**：copp 目前聚焦单段求解（`code_reading_guide.md` / paper_notes §12）。§8 的三窗调度需在 robot_copp 层新建，或反哺进 copp crate。**建议先实现整段离线（§7），HLAW 作为第二阶段。**
2. **JointMove 与 Cartesian 指令交接的 IK 分支一致性**：FK→blend→IK 往返可能选到不同构型；需在 §4.5 桥接处固定 seed 并校验 $\|\Delta\mathbf q\|$ 跳变阈值（robot6dof §7.2）。
3. **姿态过渡的流形处理**：SLERP 与五次 Hermite 在旋转流形上的联合最优仍需细化（切空间线性化 vs 单位四元数上直接优化）。
4. **约束域边界与自适应网格的强一致**：论文 Theorem 1 要求网格含所有约束域边界；需在 §5.1 与 §6 之间加一致性校验，否则误差界不成立。
5. **非静止边界的可行性**：HLAW Theorem 3 以静止边界为充分条件；首尾窗口/在线插入指令的非静止衔接可能需要论文 [32] 的增量线性化，属边界工程。
6. **实现载体最终选型**：三条并行路线（Rust copp / Python copp / MATLAB OpenCN）中落哪一条，决定后需据 §10 模块树细化具体接口——本文档保持中立。
7. **copp Rust 库需补 PLP 目标与 SPLP 循环**：现成 `topp3_lp` 是论文 TOTP-LP 基线（$\max\int a$、单次求解），并非本方案。落地论文 TOTP-SPLP 的**默认路线**（§7.2.1(a)）：扩展 `topp3_lp` 的决策向量与约束矩阵——增补辅助变量 $J_k$、把目标换成 PLP 割线上包络（eq.27 → $\min\sum(u_{k+1}-u_{k-1})J_k$）、加 eq.29d 约束行与下界 $a_k\ge\delta_{k,0}$，再外套 Algorithm 2 迭代（eq.30）。这是纯 LP、**算力最省**的实现，最契合在线机器人场景。`topp3_socp`（精确目标）作为备选。jerk 线性化（eq.32）与解析插值均可直接复用 copp 现成实现。这是本模块相对现成 copp API 的**主要新增工作量**。

---

*文档版本：v0.5（设计方案）｜作者视角：robot_copp 模块｜依赖：copp 内核 + 最优 Hermite 过渡*
*v0.2：§6.1 移除 CNC 式切向/法向约束（关节侧无物理意义、笛卡尔侧纯姿态运动退化）。*
*v0.3：TCP 约束收敛为仅两项——位置速度模长 + 姿态角速度模长，均为 $a$ 的线性上界，暂不引入 TCP 加速度/jerk。*
*v0.4：§7 求解层改为忠实论文 TOTP-SPLP——分段线性目标(PLP,eq.27)+序列线性化(Algorithm 2,eq.30)；新增 §7.2.1 落地映射。*
*v0.5：默认落地路线改为 (a) `topp3_lp`+自建 PLP 目标（纯 LP、算力最省，即论文原版 TOTP-SPLP），`topp3_socp` 精确目标降为备选；§2/§7/§10/§12 同步。*
