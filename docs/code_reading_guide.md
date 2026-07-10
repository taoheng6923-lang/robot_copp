# COPP 项目代码阅读指南

> 适用版本：`copp v0.1.0`  
> 参考论文：Wang et al., *Online time-optimal trajectory planning along parametric toolpaths with strict constraint satisfaction and certifiable feasibility guarantee*, IJMTM, 2026.

---

## 目录

1. [项目概述](#1-项目概述)
2. [核心数学变量速查](#2-核心数学变量速查)
3. [整体目录结构](#3-整体目录结构)
4. [模块依赖与调用层次](#4-模块依赖与调用层次)
5. [路径层（path/）](#5-路径层-path)
6. [机器人与约束层（robot/ + constraints.rs）](#6-机器人与约束层-robot--constraintsrs)
7. [2阶求解器（copp2/）](#7-2阶求解器-copp2)
8. [3阶求解器（copp3/）](#8-3阶求解器-copp3)
9. [目标函数（objectives.rs）](#9-目标函数-objectivesrs)
10. [Clarabel 后端（clarabel_backend.rs）](#10-clarabel-后端-clarabel_backendrs)
11. [完整端到端调用流程](#11-完整端到端调用流程)
12. [约束映射公式速查](#12-约束映射公式速查)
13. [关键设计决策说明](#13-关键设计决策说明)
14. [诊断层（diag/）](#14-诊断层-diag)
15. [数学内核（math/）](#15-数学内核-math)
16. [样条路径（path/spline.rs）](#16-样条路径-pathspliners)
17. [Clarabel 约束矩阵组装](#17-clarabel-约束矩阵组装-copp2--copp3)

---

## 1. 项目概述

**COPP（Convex-Objective Path Parameterization）** 是一个 Rust 轨迹规划库，解决以下问题：

> 给定几何路径 `q(s)`（关节角度沿路径参数的函数），在满足速度/加速度/急动度/力矩约束的前提下，找到最优时间律 `s(t)`，使代价函数最小化。

### 求解器家族

| 求解器 | 文件 | 约束阶数 | 目标函数 | 算法 |
|--------|------|----------|----------|------|
| `topp2_ra` | `copp2/dp2/topp2_ra.rs` | 2阶 | 时间最优 | DP 可达集 + 前向贪心 |
| `copp2_socp` | `copp2/opt2/copp2_socp.rs` | 2阶 | 凸目标 | Clarabel SOCP |
| `reach_set2_*` | `copp2/dp2/reach_set2.rs` | 2阶 | 可行性分析 | DP 可达集 |
| `topp3_lp` | `copp3/opt3/topp3_lp.rs` | 3阶 | 时间最优 | SCP + Clarabel LP |
| `topp3_socp` | `copp3/opt3/topp3_socp.rs` | 3阶 | 时间最优 | SCP + Clarabel SOCP |
| `copp3_socp` | `copp3/opt3/copp3_socp.rs` | 3阶 | 凸目标 | SCP + Clarabel SOCP |

---

## 2. 核心数学变量速查

| 数学符号 | 含义 | 代码变量 | 单位 |
|----------|------|----------|------|
| $s$ | 路径参数（弧长归一化） | `s_grid: Vec<f64>` | 无量纲 |
| $\mathbf{q}(s)$ | 关节位置沿路径的函数 | `derivs.q` (dim×N) | rad |
| $\mathbf{q}'(s) = d\mathbf{q}/ds$ | 路径一阶导 | `derivs.dq` | rad |
| $\mathbf{q}''(s) = d^2\mathbf{q}/ds^2$ | 路径二阶导 | `derivs.ddq` | rad |
| $\mathbf{q}'''(s) = d^3\mathbf{q}/ds^3$ | 路径三阶导 | `derivs.dddq` | rad |
| $a(s) = \dot{s}^2$ | 路径速度平方 | `a_profile: Vec<f64>` | $\text{s}^{-2}$ |
| $b(s) = \ddot{s}$ | 路径加速度 | `b_profile: Vec<f64>` | $\text{s}^{-2}$ |
| $c(s) = \dddot{s}/\dot{s}$ | 路径归一化急动度 | LP/SOCP 决策变量 | $\text{s}^{-2}$ |
| $\dot{s} = \sqrt{a}$ | 路径速度 | — | $\text{s}^{-1}$ |

### 关节空间与路径空间的关系

$$
\dot{\mathbf{q}} = \mathbf{q}' \cdot \dot{s} = \mathbf{q}' \cdot \sqrt{a}
$$

$$
\ddot{\mathbf{q}} = \mathbf{q}'' \cdot a + \mathbf{q}' \cdot b
$$

$$
\dddot{\mathbf{q}} = \mathbf{q}''' \cdot a\dot{s} + 3\mathbf{q}'' \cdot b\dot{s} + \mathbf{q}' \cdot c\dot{s}
= \sqrt{a}\left(\mathbf{q}'''a + 3\mathbf{q}''b + \mathbf{q}'c\right)
$$

---

## 3. 整体目录结构

```
src/
├── lib.rs                      公开 API 门面，定义 solver 命名空间与 prelude
│
├── path/                       路径层：将几何路径解析为导数矩阵
│   ├── autodiff.rs             Jet3：三阶前向自动微分标量
│   ├── path_core.rs            Path 结构体（参数化 / 样条统一接口）
│   └── spline.rs               Hermite 样条（O(n) 块三对角，Rayon 并行）
│
├── robot/                      机器人层：物理约束 → 路径域不等式
│   └── robot_core.rs           Robot<M> 包装器 + 约束摄入 API
│
├── copp/                       核心求解命名空间（crate 内部 pub(crate)）
│   ├── mod.rs                  内部子模块路由
│   ├── constraints.rs          Constraints：环形缓冲约束存储（1/2/3阶）
│   ├── objectives.rs           CoppObjective 目标枚举
│   ├── general.rs              InterpolationMode + 浮点近似比较工具
│   ├── clarabel_backend.rs     Clarabel 求解器选项封装与解提取
│   │
│   ├── copp2/                  2阶求解器
│   │   ├── formulation.rs      Topp2Problem / Copp2Problem 构造器
│   │   ├── interpolation.rs    a↔b 转换，s↔t 映射（后处理）
│   │   ├── dp2/                动态规划后端
│   │   │   ├── reach_set2.rs   后向/双向可达集（DP 核心）
│   │   │   └── topp2_ra.rs     RA 前向贪心求解
│   │   └── opt2/               SOCP 优化后端
│   │       ├── clarabel_constraints.rs  约束矩阵组装
│   │       └── copp2_socp.rs   COPP2 Clarabel 接口
│   │
│   └── copp3/                  3阶求解器
│       ├── formulation.rs      Topp3Problem / Copp3Problem（含线性化）
│       ├── interpolation.rs    3阶 s↔t 映射
│       └── opt3/               LP / SOCP 优化后端
│           ├── clarabel_constraints.rs  约束矩阵组装
│           ├── topp3_lp.rs     TOPP3 LP 接口
│           ├── topp3_socp.rs   TOPP3 SOCP 接口
│           └── copp3_socp.rs   COPP3 SOCP 接口
│
├── math/                       数值内核
│   └── numerical/
│       ├── lp.rs               1D/2D LP 求解器（RA 用）
│       └── general.rs          辅助数学函数
│
└── diag/                       诊断层
    ├── error.rs                CoppError + ConstraintError + PathError
    └── diagnostics.rs          Verbosity + 日志宏
```

---

## 4. 模块依赖与调用层次

```
lib.rs（公开 API）
  │
  ├─ path/          ─────────────────────────────────┐
  │   autodiff.rs → path_core.rs → PathDerivatives   │
  │                                                   │ q/dq/ddq/dddq 矩阵
  ├─ robot/         ─────────────────────────────────┘
  │   robot_core.rs: Robot<M>
  │     with_s / with_q / with_axial_velocity / _acceleration / _jerk / _torque
  │     → 写入 Constraints（环形缓冲）
  │                       │
  ├─ copp/constraints.rs  │ 被 Topp2/3ProblemBuilder 引用
  │                       ↓
  ├─ copp2/formulation.rs: Topp2Problem / Copp2Problem
  │       ↓                      ↓
  │   dp2/reach_set2.rs    opt2/copp2_socp.rs
  │   dp2/topp2_ra.rs       → clarabel_backend.rs → Clarabel 求解器
  │       ↓
  │   copp2/interpolation.rs: s_to_t / t_to_s
  │
  └─ copp3/formulation.rs: Topp3Problem / Copp3Problem
          ↓
      opt3/topp3_lp.rs / topp3_socp.rs / copp3_socp.rs
          → clarabel_backend.rs → Clarabel 求解器
          ↓
      copp3/interpolation.rs: s_to_t / t_to_s
```

---

## 5. 路径层（path/）

### 5.1 自动微分 `Jet3`（autodiff.rs）

**职责**：将路径函数 `q(s)` 的闭包表达式自动求出三阶导数，无需手动推导。

```rust
// 用户仅需写出 q(s) 的解析形式：
Path::from_parametric(|s: Jet3| vec![
    sin(2.0 * PI * s),   // 关节 0
    sin(3.0 * PI * s),   // 关节 1
], 0.0, 1.0)
```

**调用流程**：

```
Path::from_parametric(closure)
  └─ eval_parametric()  对每个 s[j] 调用：
       vals = closure(Jet3::seed(s[j]))  ← 种子：f=s, f'=1, f''=0, f'''=0
       ↓ 闭包内算术运算自动传播导数（链式法则）
       jet.v  → q_col[i]     位置
       jet.d1 → dq_col[i]    一阶导 dq/ds
       jet.d2 → ddq_col[i]   二阶导
       jet.d3 → dddq_col[i]  三阶导
```

**关键公式（以 sin 为例）**，设 $u = u(s)$：

$$
f = \sin u, \quad f' = \cos u \cdot u'
$$
$$
f'' = -\sin u \cdot (u')^2 + \cos u \cdot u''
$$
$$
f''' = -\cos u \cdot (u')^3 - 3\sin u \cdot u' u'' + \cos u \cdot u'''
$$

所有基本函数（sin/cos/exp/ln/sqrt/powi）均用同一思路实现链式法则，
运算符重载（Add/Sub/Mul/Div）也各自实现了相应的导数传播。

### 5.2 路径核心 `Path`（path_core.rs）

**两种后端**：

| 后端 | 构建方式 | 求导方式 |
|------|----------|----------|
| `Parametric` | `from_parametric(closure)` | `Jet3` 前向 AD |
| `Spline` | `from_waypoints(&waypoints, cfg)` | 多项式系数直接求导 |

**求值方法层次**（按计算量递增）：

```rust
path.evaluate_q(&s)           // 仅位置，dq=ddq=dddq=None
path.evaluate_up_to_2nd(&s)   // 位置+速度+加速度导数（TOPP2 所需）
path.evaluate_up_to_3rd(&s)   // 全部三阶（TOPP3/COPP3 所需）
```

所有方法底层均调用 `evaluate_impl(s, Order)`，再通过 Rayon 并行分发。

---

## 6. 机器人与约束层（robot/ + constraints.rs）

### 6.1 Robot<M> 包装器（robot_core.rs）

**职责**：将机器人关节空间的物理约束（速度/加速度/急动度/力矩）
转换为路径域的不等式约束，写入 `Constraints` 环形缓冲区。

**类型参数**：
- `M: RobotBasic`：只需 `dim()` 方法，适用于 TOPP 问题
- `M: RobotTorque`：额外需要逆动力学 `inverse_dynamics(q, dq, ddq, tau)`，适用于 COPP 问题

**约束转换公式**：

#### 速度约束 → 1阶约束

关节速度约束 $v_{\min} \le \dot{q}_i \le v_{\max}$，由 $\dot{q}_i = q'_i \cdot \dot{s} = q'_i\sqrt{a}$ 映射为：

$$
a \le a_{\max,i} = \begin{cases} (v_{\max}/q'_i)^2 & q'_i > 0 \\ (v_{\min}/q'_i)^2 & q'_i < 0 \\ +\infty & q'_i = 0 \end{cases}
$$

#### 加速度约束 → 2阶约束

关节加速度约束 $\alpha_{\min} \le \ddot{q}_i \le \alpha_{\max}$，由 $\ddot{q}_i = q''_i a + q'_i b$ 映射为：

$$
q''_i \cdot a + q'_i \cdot b \le \alpha_{\max,i}
$$
$$
-q''_i \cdot a - q'_i \cdot b \le -\alpha_{\min,i}
$$

代码符号：`acc_a = q''`，`acc_b = q'`。

#### 急动度约束 → 3阶非线性约束

关节急动度约束 $j_{\min} \le \dddot{q}_i \le j_{\max}$，由

$$
\dddot{q}_i = \sqrt{a}\left(q'''_i a + 3q''_i b + q'_i c\right)
$$

映射为：

$$
\sqrt{a}\left(\text{jerk\_a}_i \cdot a + \text{jerk\_b}_i \cdot b + \text{jerk\_c}_i \cdot c\right) \le j_{\max,i}
$$

其中 $\text{jerk\_a} = q'''$，$\text{jerk\_b} = 3q''$，$\text{jerk\_c} = q'$（代码中通过 `jerk_b_new.scale_mut(3.0)` 应用系数 3）。

#### 力矩约束 → 2阶约束（需逆动力学）

力矩分解：

$$
\tau = M(\mathbf{q})(\mathbf{q}''a + \mathbf{q}'b) + C(\mathbf{q}, \mathbf{q}'\sqrt{a})\mathbf{q}'\sqrt{a} + \mathbf{g}(\mathbf{q})
\approx \text{coeff\_a} \cdot a + \text{coeff\_b} \cdot b + \text{coeff\_g}
$$

通过三次逆动力学调用计算系数：

$$
\text{coeff\_g} = \tau(\mathbf{q}, \mathbf{0}, \mathbf{0}), \quad
\text{coeff\_b} = \tau(\mathbf{q}, \mathbf{0}, \mathbf{q}') - \text{coeff\_g}, \quad
\text{coeff\_a} = \tau(\mathbf{q}, \mathbf{q}', \mathbf{q}'') - \text{coeff\_g}
$$

约束行：$\text{coeff\_a} \cdot a + \text{coeff\_b} \cdot b \le \tau_{\max} - \text{coeff\_g}$。

### 6.2 约束存储 Constraints（constraints.rs）

**存储结构**：环形缓冲列矩阵，物理列索引 = `(head_col + 逻辑索引) % capacity_col`。

**约束家族**：

| 名称 | 数学形式 | 缓冲区字段 |
|------|----------|------------|
| 1阶上界 | $0 \le a_k \le a_{\max,k}$ | `amax` |
| 2阶行 | $f_a \cdot a_k + f_b \cdot b_k \le f_{\max}$ | `acc_a, acc_b, acc_max` |
| 3阶非线性行 | $\sqrt{a_k}(g_a a_k + g_b b_k + g_c c_k + g_d) \le g_{\max}$ | `jerk_{a,b,c,d,max}` |
| 3阶线性化行 | $h_a a_k + h_b b_k + h_c c_k \le h_{\max}$ | `jerk_a_linear, jerk_max_linear` |

---

## 7. 2阶求解器（copp2/）

### 7.1 TOPP2-RA 完整调用流程

```
用户调用 topp2_ra(problem, options)
  │
  ├─ 按 verbosity 分发到 topp2_ra_core(...)
  │
  ├─ 步骤1：reach_set2_backward(problem, options)
  │    ├─ 初始化：$a_{\max}[n]=a_{\text{final}}$，$a_{\max}[0..n-1]=+\infty$
  │    │
  │    └─ 后向 DP 循环 $k = n-1 \to 0$：
  │         fill_acc_topp2 收集段 $[s_k, s_{k+1}]$ 上的约束行：
  │           $f_a \cdot a_{k+1} + f_b \cdot a_k \le f_{\max}$
  │         lp_2d_incre_max_y: 在 $a_{k+1}\in[a_{\min},a_{\max}]$ 下求 $a_k$ 的最大值（2D LP）
  │         裁剪：$a_{\max}[k] = \min(\text{LP解},\; a_{\max,k})$
  │         返回 ReachSet2 { a_max, a_min }
  │
  └─ 步骤2：前向贪心遍历 $k = 1,\ldots,n$
       ① 收集段 $[s_{k-1}, s_k]$ 的约束行
       ② 代入 $a_{k-1}=a_{\text{prev}}$，化为 1D LP 求 $a_k$ 上界
       ③ 与后向区间取交：$a_{\max,\text{curr}} = \min(\text{LP},\; a_{\max}[k])$
       ④ 贪心选择（时间最优关键步）：
            $a_k = a_{\max,\text{curr}}$（最大化 $\dot{s}$，最小化 $dt = ds/\dot{s}$）

  返回 $a$: `Vec<f64>`（长度 $n+1$，$a_k = \dot{s}_k^2$）
```

### 7.2 TOPP2 后处理（interpolation.rs）

$a$ profile 的后处理链：

**`a_to_b_topp2`**（有限差分）：
$$b_k = \frac{a_{k+1} - a_k}{2(s_{k+1}-s_k)}, \quad k=0,\ldots,n-2$$

**`s_to_t_topp2`**（梯形时间积分）：
$$\Delta t_k = \frac{2(s_{k+1}-s_k)}{\sqrt{a_k}+\sqrt{a_{k+1}}}, \quad t_k = t_0 + \sum_{j<k}\Delta t_j$$

**`t_to_s_topp2`**（段内解析反演，$a(s) = c_0 + c_1 x$ 线性假设，$x$ 为段内偏移）：
$$x_{\text{right}} = \frac{\left(\sqrt{c_0 + c_1 x_{\text{left}}} + \frac{c_1 \Delta t}{2}\right)^2 - c_0}{c_1}$$

### 7.3 COPP2-SOCP（copp2_socp.rs）

将问题转化为 Clarabel SOCP 后求解，目标函数由 `CoppObjective` 列表决定：

```
copp2_socp(problem, options)
  ├─ clarabel_standard_capacity_copp2()   预估约束矩阵规模
  ├─ clarabel_standard_constraint_copp2() 组装稀疏矩阵 A, P, q
  │    └─ 时间目标线性化（SOC辅助变量引入 1/√a 的锥重表达）
  ├─ DefaultSolver::new().solve()          Clarabel 内点法
  └─ clarabel_to_copp2_solution()          从 x 向量提取 a profile
```

---

## 8. 3阶求解器（copp3/）

### 8.1 核心：非线性约束线性化（formulation.rs）

3阶求解器需要**先运行 TOPP2-RA 得到参考点 `a_ref`**，再执行 SCP 迭代：

**原始非线性 jerk 约束**（约束存储中的 `jerk_*` 字段）：

$$
\sqrt{a}\left(g_a a + g_b b + g_c c + g_d\right) \le g_{\max}
$$

**在 $a_{\text{lin}}$ 处的线性化**（写入 `jerk_a_linear` / `jerk_max_linear`）：

利用 $\sqrt{a} \approx \sqrt{a_{\text{lin}}} + \dfrac{a - a_{\text{lin}}}{2\sqrt{a_{\text{lin}}}}$，展开整理后：

$$
\underbrace{\left(g_a + \frac{g_{\max}}{2\,a_{\text{lin}}^{3/2}}\right)}_{h_a} a
+ g_b b + g_c c
\le \underbrace{\frac{3}{2}\frac{g_{\max}}{\sqrt{a_{\text{lin}}}} - g_d}_{h_{\max}}
$$

其中 $a_{\text{lin,eff}} = \max(a_{\text{lin}},\; a_{\text{linearization\_floor}})$ 防止除零（静止边界附近）。

**两次 SCP 迭代过程**：

```
【第0步】TOPP2-RA → a_ref（2阶时间最优解）
  ↓
【第1次 SCP】
  robot.constraints.amax_substitute(&a_ref, 0)   用 a_ref 收紧速度上界
  Topp3ProblemBuilder::new(&mut robot, 0, &a_ref, ...)
    .build_with_linearization()
      → linearize_constraint_3order_with_floor(&a_ref, ...)  写入线性化系数
      → 返回 Topp3Problem（只读，线性化已就位）
  topp3_lp / topp3_socp → 求解 → (a1, b1, num_stat1)
  ↓
【第2次 SCP】（用 a1 替换 a_ref）
  Topp3ProblemBuilder::new(&mut robot, 0, &a1, ...)
    .build_with_linearization()  重新线性化
  topp3_lp / topp3_socp → 求解 → (a2, b2, num_stat2)  更优解
```

### 8.2 静止边界建模（num_stationary）

当边界条件 $a=b=0$（起点或终点静止），常规 $c=\dddot{s}/\dot{s}$ 公式在 $\dot{s}\to 0$ 时退化。
解决方案：在边界附近加入若干"常急动度（jerk-constant）"段作为静止缓冲区：

$$\text{num\_stationary} = (\text{start},\; \text{end})$$

各端的缓冲站数由 `determine_num_stationary_pair` 自动确定：
- 若 $a_{\text{boundary}} \approx 0$ 且 $b_{\text{boundary}} \approx 0$ → 使用 `num_stationary_max`
- 否则 → 缓冲站数为 0

### 8.3 TOPP3 后处理（copp3/interpolation.rs）

原理与 2 阶类似，但节点 $b_k$ 为节点值（非段值），时间积分更精确：

**`s_to_t_topp3`**：对每个非静止段用含 $b$ 的局部模型积分倒二次根积分：

$$
\Delta t = \int_{s_k}^{s_{k+1}} \frac{ds}{\sqrt{c_0 + c_1 x + c_2 x^2}}, \quad
c_0 = a_k,\; c_1 = 2b_k,\; c_2 = \frac{b_{k+1}-b_k}{s_{k+1}-s_k}
$$

对静止段用匀急动度（jerk-constant）三次律 $s(t) = s_0 + \frac{1}{6}\dddot{s}_0 t^3$ 积分。

**`t_to_s_topp3`**：均匀时间采样，通过反函数 `inverse_rsrqp` 段内解析求逆。

---

## 9. 目标函数（objectives.rs）

`CoppObjective` 枚举仅在 COPP 系列（非 TOPP）中使用，通过 `Copp2/3ProblemBuilder` 传入：

```rust
let objectives = [
    CoppObjective::Time(1.0),                      // 最小化时间，权重 1.0
    CoppObjective::ThermalEnergy(0.1, &weights),   // 最小化热能，权重 0.1
];
```

| 枚举变体 | 连续形式 | 适用场景 |
|----------|----------|----------|
| `Time(w)` | $w\displaystyle\int \frac{ds}{\sqrt{a}}$ | 基础时间最优 |
| `ThermalEnergy(w, ν)` | $w\displaystyle\int \frac{\sum_i(\tau_i\nu_i)^2}{\sqrt{a}}\,ds$ | 减少电机发热 |
| `TotalVariationTorque(w, ν)` | $w\displaystyle\sum_i\int\left|\frac{d\tau_i}{ds}\right|\nu_i\,ds$ | 平滑力矩曲线 |
| `Linear(w, α, β)` | $w\displaystyle\int(\alpha a+\beta b)\,ds$ | 用户自定义代价 |

---

## 10. Clarabel 后端（clarabel_backend.rs）

**职责**：统一封装 Clarabel 求解器的配置、求解调用和结果提取。

**调用位置**：所有 SOCP/LP 求解器（`copp2_socp`, `topp3_lp`, `topp3_socp`, `copp3_socp`）

```
ClarabelOptionsBuilder::new()
  .allow_almost_solved(true)   // 接受"近似求解"状态（生产推荐）
  .allow_max_iterations(true)  // 接受迭代次数达上限的解
  .build()
  → ClarabelOptions { clarabel_settings, allow_*, verbosity }

求解器调用：
  DefaultSolver::new(&P, &q, &A, &b, &cones, &settings).solve()
  → DefaultSolution<f64> { x, status, ... }
  → options.is_allow(status)  决定是否接受此解
  → clarabel_to_copp2_solution(&solution, ...)  提取 a profile
  → clarabel_to_copp3_solution(&solution, ...)  提取 (a, b) + 静止段填充
```

---

## 11. 完整端到端调用流程

以 TOPP2-RA + COPP3-SOCP（两次 SCP）为例：

```
━━━━━━━━━━━ 步骤 1：路径构建 ━━━━━━━━━━━
Path::from_parametric(|s: Jet3| vec![sin(2πs), sin(3πs), sin(5πs)], 0.0, 1.0)
  └─ 存储闭包，记录 dim=3

━━━━━━━━━━━ 步骤 2：路径求值 ━━━━━━━━━━━
path.evaluate_up_to_3rd(&s_grid)
  └─ 对每个 s[j]，调用 closure(Jet3::seed(s[j]))
  └─ 提取 jet.v/d1/d2/d3 → PathDerivatives { q(3×n), dq, ddq, dddq }

━━━━━━━━━━━ 步骤 3：约束摄入 ━━━━━━━━━━━
let mut robot = Robot::with_capacity(3usize, n);
robot.with_s(&s_grid)
  → Constraints::with_s()  写入路径参数
robot.with_q(&q, &dq, &ddq, Some(&dddq), 0)
  → Constraints::with_q()  写入路径导数矩阵
robot.with_axial_velocity((vel_max, n), (vel_min, n), 0)
  → 计算 amax = (v/q')²，写入 Constraints::amax
robot.with_axial_acceleration((acc_max, n), (acc_min, n), 0)
  → acc_a=q'', acc_b=q', 写入 acc_{a,b,max}
robot.with_axial_jerk((jerk_max, n), (jerk_min, n), 0)
  → jerk_a=q''', jerk_b=3q'', jerk_c=q', 写入 jerk_{a,b,c,d,max}

━━━━━━━━━━━ 步骤 4：TOPP2-RA 求解 ━━━━━━━━━━━
let topp2_problem = Topp2ProblemBuilder::new(&robot, (0, n-1), (0.0, 0.0)).build()?;
let opts_ra = ReachSet2OptionsBuilder::new().build()?;
let a_ref = topp2_ra(&topp2_problem, &opts_ra)?;
  → 后向 DP 可达集 + 前向贪心 → a_ref: Vec<f64>（2阶最优 a profile）

━━━━━━━━━━━ 步骤 5：收紧速度上界 ━━━━━━━━━━━
robot.constraints.amax_substitute(&a_ref, 0)?;
  → 每站 amax[k] = min(amax[k], a_ref[k])
  → 为 3 阶线性化提供更紧的参考域

━━━━━━━━━━━ 步骤 6：COPP3 第1次 SCP ━━━━━━━━━━━
let copp_obj = [CoppObjective::Time(1.0), CoppObjective::ThermalEnergy(0.1, &weights)];
let problem_c3 = Copp3ProblemBuilder::new(&mut robot, &copp_obj, 0, &a_ref, (0.0,0.0), (0.0,0.0))
  .build_with_linearization()?;
  → linearize_constraint_3order_with_floor(&a_ref, ...)
       将 jerk 约束在 a_ref 处线性化 → 写入 jerk_a_linear / jerk_max_linear
  → 返回 Copp3Problem { robot, objectives, constraints（含线性化） }
let (a1, b1, num_stat1) = copp3_socp(&problem_c3, &clarabel_opts)?;
  → 组装 Clarabel SOCP → 求解 → 提取 (a1, b1)

━━━━━━━━━━━ 步骤 7：COPP3 第2次 SCP ━━━━━━━━━━━
let problem_c3_2 = Copp3ProblemBuilder::new(&mut robot, &copp_obj, 0, &a1, (0.0,0.0), (0.0,0.0))
  .build_with_linearization()?;
  → 用 a1 重新线性化（更接近真实约束）
let (a2, b2, num_stat2) = copp3_socp(&problem_c3_2, &clarabel_opts)?;

━━━━━━━━━━━ 步骤 8：转换为时间域轨迹 ━━━━━━━━━━━
let (t_final, t_s) = s_to_t_topp3(&s_grid, &a2, &b2, num_stat2, 0.0);
let s_t = t_to_s_topp3(&s_grid, &a2, &b2, num_stat2, &t_s,
                         InterpolationMode::UniformTimeGrid(0.0, 1e-3, true));
  → s_t[i] = s(t_i)，均匀时间采样

━━━━━━━━━━━ 步骤 9：关节空间重建（Python 脚本） ━━━━━━━━━━━
q(t)  = interp(s(t), s_grid, q_grid)          位置
q̇(t)  = q'(s(t)) · √a(s(t))                   速度
q̈(t)  = q''(s(t)) · a(s(t)) + q'(s(t)) · b(t) 加速度（b 从 profile 或差分得到）
```

---

## 12. 约束映射公式速查

| 约束类型 | 关节空间 | 路径域约束行 | 代码 API |
|----------|----------|-------------|----------|
| 速度 | $v_{\min} \le \dot{q} \le v_{\max}$ | $a \le (v_{\max}/q')^2$ | `with_axial_velocity` |
| 加速度 | $\alpha_{\min} \le \ddot{q} \le \alpha_{\max}$ | $q'' a + q' b \le \alpha$ | `with_axial_acceleration` |
| 急动度 | $j_{\min} \le \dddot{q} \le j_{\max}$ | $\sqrt{a}(q'''a+3q''b+q'c) \le j$ | `with_axial_jerk` |
| 力矩 | $\tau_{\min} \le \tau \le \tau_{\max}$ | $c_a a + c_b b \le \tau_{\max} - g(\mathbf{q})$ | `with_axial_torque` |

---

## 13. 关键设计决策说明

### 为什么 TOPP3 需要 `&mut robot`？

`Topp3ProblemBuilder::new(&mut robot, ...)` 调用 `.build_with_linearization()` 时，
会将线性化后的 jerk 系数**原地写入** `robot.constraints`（`jerk_a_linear` / `jerk_max_linear` 字段），
而 2阶构造器只需只读引用 `&robot`。这是接口可变性差异的唯一原因。

### 为什么用环形缓冲区存储约束？

`Constraints` 使用环形缓冲列矩阵，支持在线/窗口化规划：
- `pop_front()` / `pop_back()` 以 O(1) 移除旧站点
- `expand_capacity()` 动态扩容（类似 Vec 的分摊 O(1)）
- 适合滑动窗口场景（如 CNC 加工中每个刀路段分批处理）

### 为什么 TOPP2-RA 在理论上保证全局最优？

2阶约束形如 $f_a a_{k+1} + f_b a_k \le f_{\max}$，是关于 $(a_k, a_{k+1})$ 的线性约束。
动态规划的**最优子结构**成立：若 $a_0,\ldots,a_k$ 已最优且 $a_k$ 最大，则在 $a_k$ 条件下 $a_{k+1}$ 贪心最优。
因此，后向可达集 + 前向贪心构成完整的 DP 策略，保证全局时间最优。

### SOCP 求解器选择：Clarabel

项目使用 [Clarabel](https://clarabel.org/)（Rust 原生内点法 SOCP/LP 求解器），
相比 ECOS/MOSEK 的优势：
- 纯 Rust 无外部依赖，适合嵌入工业控制器
- 支持锥约束（SOC + 非负锥），覆盖所有 TOPP/COPP 问题形式
- `ClarabelOptionsBuilder` 提供 `allow_almost_solved` 等策略，可在精度和速度间权衡

---

## 14. 诊断层（diag/）

### 14.1 错误类型（error.rs）

库统一使用 `CoppError` 作为顶层错误类型，内部子错误通过 `From` 自动转换：

```
CoppError
  ├─ ConstraintError  约束存储/访问错误（outOfBounds, 维度不匹配等）
  ├─ PathError        路径参数越界、样条次数无效等
  ├─ InvalidInput     用户输入参数校验失败（如权重为负）
  ├─ ClarabelSolverStatus  Clarabel 返回非成功状态（Infeasible 等）
  └─ Other            其他运行时错误
```

**使用模式**：
```rust
fn my_func() -> Result<(), CoppError> {
    let s = constraints.get_s(idx)?;  // ConstraintError 自动 → CoppError
    // ...
}
```

### 14.2 诊断工具（diagnostics.rs）

**Verbosity 四级**：

| 级别 | 信息量 | 典型用途 |
|------|--------|----------|
| `Silent` | 无输出 | 生产/嵌入式 |
| `Summary` | 关键里程碑 + 总耗时 | 日常调试 |
| `Debug` | 矩阵规模、各阶段行数 | 算法验证 |
| `Trace` | 逐步差量、解向量头部 | 精细诊断 |

**日志目标 `VerbosityOutput`**：
- `Println` — 控制台（默认）
- `Log` — `log` crate facade（需先初始化全局 logger）
- `File(path)` — 追加到日志文件

**Verboser trait（编译期分发）**：

求解核心函数如 `topp3_lp_core` 通过泛型参数 `impl Verboser` 接受具体实现，
在 `Silent` 分支下编译器可以消除所有日志分支，零运行时开销。

---

## 15. 数学内核（math/）

### 15.1 增量 LP 求解器（math/numerical/lp.rs）

TOPP2-RA 的后向 DP 步骤需要在每个路径站求解：

> 在所有约束 $a_i x + b_i y \le c_i$ 的交集中，找到最大化 $y$ 的点。

这是一个**2D LP 问题**，本模块用 **Seidel 增量算法**（随机化，期望 $O(m)$）求解：

**算法思路**：
1. 初始解：`(x*, y*) = (0, +∞)`（无约束上界）
2. 对每条新约束 `a_i·x + b_i·y ≤ c_i`：
   - 若当前解满足，跳过
   - 否则当前最优一定在新约束边界 `a_i·x + b_i·y = c_i` 上 → 退化为 1D LP

**关键常数**：
- `EPS_ZERO = 1e-9` — 近零分支决策阈值
- `LP_BOUND = 1e6` — 无界方向的默认盒子约束
- `EPS_SCALE = 10.0` — 维度规约时的容差缩放

### 15.2 辅助函数（math/numerical/general.rs）

- `cross_product_2d(x1, x2)` — 2D 叉积，被 LP 几何内核用于判断半空间方向
- `solve_2x2(A, b)` — 2×2 线性系统，用于 `force_positive_a` 的局部矫正

---

## 16. 样条路径（path/spline.rs）

当用户调用 `Path::from_waypoints(waypoints, config)` 时，内部调用此模块：

**支持阶次（奇次，`p = 2m+1`）**：

| 阶次 p | m | 连续性 | 边界条件 |
|--------|---|--------|----------|
| 3 | 1 | C² | 两端速度 |
| 5 | 2 | C⁴ | 两端速度+加速度 |
| 7 | 3 | C⁶ | 两端速度+加速度+急动度 |

**求解算法（O(N) 块 Thomas）**：

```
Hermite 参数化 + C^{m+1}..C^{2m} 连续性
  → m×m 块三对角线性系统
  → 前向消元（块 LU 分解）
  → 回代求各段多项式系数
```

**求值**：多项式系数直接对路径参数 s 求导，速度比 Jet3 AD 快（无重复计算）。

---

## 17. Clarabel 约束矩阵组装（copp2/ + copp3/）

### 17.1 TOPP2 约束矩阵

`clarabel_standard_constraint_topp2()` 组装决策变量 `x = [a[0..=n]]`：

| 约束块 | 锥类型 | 数学含义 |
|--------|--------|----------|
| 边界等式 | `ZeroConeT` | $a_0 = a_{\text{start}},\; a_n = a_{\text{final}}$ |
| 一阶不等式 | `NonnegativeConeT` | $a_{\max,k} - a_k \ge 0$ |
| 二阶不等式 | `NonnegativeConeT` | $f_{\max} - f_a a_{k+1} - f_b a_k \ge 0$ |

**关键**：TOPP2 中 $b_k = (a_{k+1}-a_k)/(2\Delta s)$ 通过约束行的 `acc_a, acc_b` 系数隐式编码，
$b$ 不是独立决策变量。

### 17.2 TOPP3 约束矩阵

`clarabel_standard_constraint_topp3()` 组装决策变量 `x = [a[0..=n], b[0..=n]]`：

| 约束块 | 锥类型 | 数学含义 |
|--------|--------|----------|
| 边界等式（4行） | `ZeroConeT` | $a_0,b_0,a_n,b_n$ 固定为边界值 |
| 动力学等式 | `ZeroConeT` | $a_{k+1} = a_k + 2\Delta s_k b_k$（离散导数连续性） |
| 一阶不等式 | `NonnegativeConeT` | $a_{\max,k} - a_k \ge 0$ |
| 二阶不等式 | `NonnegativeConeT` | $f_{\max} - f_a a_k - f_b b_k \ge 0$ |
| 线性化三阶不等式 | `NonnegativeConeT` | $h_{\max} - h_a a_k - h_b b_k - h_c c_k \ge 0$ |

**与 TOPP2 的关键区别**：
- 动力学等式约束 $a_{k+1} = a_k + 2\Delta s_k b_k$ 将 $a, b$ 耦合
- $c_k$ 通过 $b$ 的差分隐式表示：$c_k \approx (b_{k+1}-b_{k-1})/(2\Delta s)$（约束行中已展开）
- 静止边界段用特殊的 `set_ab_stationary_topp3()` 前置处理

---

*文档生成于项目 `copp v0.1.0`，如代码有更新请对照源文件核实。*
