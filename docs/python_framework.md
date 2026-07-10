# robot_copp：Python 算法框架（代码骨架 + 模块算法原理）

> **配套**：本文件是 [`robot_copp_design.md`](./robot_copp_design.md) 的 Python 实现框架。设计文档回答"**做什么**"（数据流、数学、模块职责），本文件回答"**怎么组织代码**"（依赖、目录、每模块算法原理、接口签名骨架）。
>
> **范围**：只给出**算法框架与软件骨架**——库依赖、目录层级、模块划分、类/函数签名（`...` 占位）、关键算法伪代码。**不含具体实现**。
>
> **符号**：与设计文档一致——路径参数 $s$，copp 状态 $a=\dot s^2,\ b=\ddot s,\ c=\dddot s\dot s$。默认求解路线为设计文档 §7.2.1(a)：**PLP 分段线性目标 + LP**（论文 TOTP-SPLP，算力最省）。

---

## 目录

1. [依赖库](#1-依赖库)
2. [顶层包结构（目录层级）](#2-顶层包结构目录层级)
3. [数据流与模块映射](#3-数据流与模块映射)
4. [核心数据类型](#4-核心数据类型)
5. [各子包算法原理与接口骨架](#5-各子包算法原理与接口骨架)
   - [5.1 robot/ 运动学-动力学后端](#51-robot-运动学-动力学后端)
   - [5.2 commands/ 指令层](#52-commands-指令层)
   - [5.3 blending/ 最优 Hermite 过渡](#53-blending-最优-hermite-过渡)
   - [5.4 lowering/ 降维到关节空间](#54-lowering-降维到关节空间)
   - [5.5 constraints/ 约束摄入](#55-constraints-约束摄入)
   - [5.6 solve/ copp 求解层（SPLP 核心）](#56-solve-copp-求解层splp-核心)
   - [5.7 hlaw/ 分层前瞻窗口](#57-hlaw-分层前瞻窗口)
   - [5.8 synth/ 轨迹合成与验证](#58-synth-轨迹合成与验证)
   - [5.9 backend/ 求解器后端抽象](#59-backend-求解器后端抽象)
   - [5.10 planner.py 门面 + 全局设施](#510-plannerpy-门面--全局设施)
6. [关键算法伪代码](#6-关键算法伪代码)
7. [端到端调用示例](#7-端到端调用示例)
8. [测试策略](#8-测试策略)
9. [实现里程碑建议](#9-实现里程碑建议)

---

## 1. 依赖库

分层依赖（沿用 [`copp/python/DESIGN.md`](../0.other_lib_code/copp/copp/python/DESIGN.md) §2 的选型）：

| 层 | 用途 | 库 | 选型说明 |
|----|------|----|----------|
| 数值基础 | 数组/线代 | `numpy` | 全项目数据交换格式（`ndarray`） |
| 数值基础 | 样条/积分/差分/LP 兜底 | `scipy` | `interpolate`（Hermite）、`optimize.linprog`（LP 兜底）、`spatial.transform`（`Rotation`/`Slerp` 姿态插值） |
| 凸优化建模 | LP/SOCP 统一建模 | `cvxpy` | solve 层组装 SPLP 的 LP（默认）与 SOCP（备选） |
| 凸优化求解 | LP/SOCP 后端 | `clarabel` | 与 Rust copp 同核，数值可对照；`glpk`/`highs`/`ecos`/`scs` 作 cvxpy 兜底 |
| 运动学 | FK/IK/Jacobian | `roboticstoolbox-python`（Peter Corke）／`ikpy`／自实现 DH-POE | 经 `KinematicsModel` 协议解耦；默认自实现 DH-POE，rtb 作可插拔适配器 |
| 动力学（可选） | 逆动力学 RNE（力矩约束/热能目标） | `roboticstoolbox-python`（`rne`）／`pinocchio`（`rnea`，可选） | 经 `DynamicsModel` 协议解耦；无力矩约束时不装 |
| 姿态几何 | 四元数/SLERP、圆弧 | `scipy.spatial.transform` | 避免自写四元数数值细节 |
| 测试 | 单元/集成 | `pytest` | — |
| 属性测试 | 随机路径批量验证 | `hypothesis` | 随机指令/约束组合，捕获边界 |
| 可视化（非本体） | 轨迹绘图 | `matplotlib` | 复用 copp `scripts/plot_*.py` 风格 |

版本建议：Python `>=3.10`（`match`/`Protocol`/`Literal`）、`numpy>=1.26`、`scipy>=1.11`、`cvxpy>=1.4`、`clarabel>=0.6`。

extras 分组（核心可最小安装）：

```toml
[project.optional-dependencies]
kinematics = ["roboticstoolbox-python"]
dynamics   = ["roboticstoolbox-python"]   # 或 pinocchio
socp       = ["clarabel"]                 # 备选精确目标路线
viz        = ["matplotlib"]
test       = ["pytest", "hypothesis"]
```

> **与现成 copp 的关系**：solve 层默认**在 Python 内用 cvxpy 自建 SPLP**（对应设计文档 §7.2.1(a)：`topp3_lp` + 自建 PLP 目标）。若已装 copp 的 Rust/Python 绑定，可通过 `backend/copp_backend.py` 直接调用其 `topp2_ra`/`topp3_lp`/`build_with_linearization`，复用其离散化与 jerk 线性化（见 §5.9）。

---

## 2. 顶层包结构（目录层级）

```
robot_copp/
├── pyproject.toml
├── README.md
├── docs/
│   ├── robot_copp_design.md          # 设计文档（已存在）
│   └── python_framework.md           # 本文档
│
├── copp/                             # 包本体
│   ├── __init__.py
│   ├── types.py                      # 全局数据类型（PoseSegment / BlendedPath / PathDerivatives / ...）
│   ├── options.py                    # BlendOptions / SampleOptions / SolveOptions / PlannerOptions
│   ├── errors.py                     # 异常层次（RobotCoppError ...）
│   ├── diagnostics.py                # Verbosity + 日志
│   │
│   ├── robot/                        # 机器人本体：运动学/动力学后端     §5.1
│   │   ├── base.py                   #   KinematicsModel / DynamicsModel Protocol
│   │   ├── synthetic.py              #   SyntheticRobotModel（M1 stand-in）
│   │   ├── dh_poe.py                 #   自实现 DH/POE 最小实现
│   │   └── rtb_adapter.py            #   roboticstoolbox 适配器
│   │
│   ├── commands/                     # 指令层                          §5.2 / 设计 §3
│   │   ├── base.py                   #   MotionCommand Protocol + PoseSegment 工厂
│   │   ├── joint_move.py             #   JointMoveCommand
│   │   ├── linear_move.py            #   LinearMoveCommand
│   │   └── circular_move.py          #   CircularMoveCommand
│   │
│   ├── blending/                     # 最优 Hermite 过渡（笛卡尔）      §5.3 / 设计 §4
│   │   ├── zone.py                   #   zone 距离 Δ_j + 截短
│   │   ├── frenet.py                 #   端点 Frenet 框架（t,n,κ）
│   │   ├── optimal_hermite.py        #   G2 五次 Hermite + α0/α1 最优化
│   │   ├── pose_blend.py             #   位姿联合过渡 + 量纲归一化 D
│   │   └── junction.py               #   混合指令桥接 / 快路径判定
│   │
│   ├── lowering/                     # 降维到关节空间                  §5.4 / 设计 §5
│   │   ├── sampling.py               #   曲率驱动自适应离散化
│   │   ├── ik.py                     #   连续解 IK
│   │   ├── derivatives.py            #   q',q'',q''' 链式法则
│   │   └── singularity.py            #   DLS + 局部加密
│   │
│   ├── constraints/                  # 约束摄入                        §5.5 / 设计 §6
│   │   ├── model.py                  #   ConstraintSet（站点索引存储）
│   │   └── ingest.py                 #   物理约束 → 路径域 1/2/3 阶不等式
│   │
│   ├── solve/                        # copp 求解层（SPLP）★核心        §5.6 / 设计 §7
│   │   ├── state.py                  #   (a,b,c) 状态 + 无损离散化系数 + 静止段
│   │   ├── seed.py                   #   种子 a⁽⁰⁾（2 阶可达/LP，= topp2_ra）
│   │   ├── linearize.py              #   jerk 凹约束切线线性化（eq.32）
│   │   ├── plp_objective.py          #   PLP 目标（eq.27 割线 + J_k + eq.29d + 下界 a_k≥δ_0）
│   │   ├── lp_problem.py             #   单次 LP 组装（cvxpy）→ (a,b,c)
│   │   ├── splp.py                   #   Algorithm 2 迭代循环（eq.30 停止）
│   │   └── interp.py                 #   解析插值 s↔t（Prop.1/2 闭式）
│   │
│   ├── hlaw/                         # 分层前瞻窗口（长序列，可选）     §5.7 / 设计 §8
│   │   ├── windows.py                #   种子/可行/最优 三窗调度
│   │   └── relay.py                  #   跨窗口边界条件传递
│   │
│   ├── synth/                        # 轨迹合成与验证                  §5.8 / 设计 §9
│   │   ├── resample.py               #   解析插值采样 → q(t),q̇,q̈,q⃛
│   │   └── verify.py                 #   R_v / D_v 超限指标
│   │
│   ├── backend/                      # 求解器后端抽象                  §5.9
│   │   ├── base.py                   #   SolverBackend Protocol
│   │   ├── cvxpy_backend.py          #   cvxpy + clarabel（默认）
│   │   └── copp_backend.py           #   可选：调用现成 copp（Rust/py 绑定）
│   │
│   └── planner.py                    # TrajectoryPlanner 门面           §5.10
│
└── tests/
    ├── conftest.py
    ├── unit/                         # 逐模块（blending / lowering / solve / interp ...）
    ├── integration/                  # 端到端（joint-only / mixed / objectives）
    └── benchmark/                    # 随机路径批量（@pytest.mark.slow）
```

**依赖方向**（单向无环）：

```
robot     ─┐
commands  ──┼─→ blending ─→ lowering ─→ constraints ─→ solve ─→ synth
            │                                            ↑
            └───────────────────────── hlaw ────────────┘（长序列时编排 solve）
planner.py 编排以上全部；backend/ 被 solve/ 依赖注入
```

---

## 3. 数据流与模块映射

| 阶段 | 输入 → 输出 | 子包 | 设计文档 |
|------|-------------|------|----------|
| 指令解析 | `list[MotionCommand]` → `list[PoseSegment]`（笛卡尔位姿段） | `commands/` | §3 |
| 过渡平滑 | `list[PoseSegment]` → `BlendedPath`（C²/分段C³ 位姿路径 r(s)） | `blending/` | §4 |
| 降维 | `BlendedPath` → `PathDerivatives`（q,q',q'',q''' 于 {s_k}） | `lowering/` | §5 |
| 约束摄入 | 物理限值 + `PathDerivatives` → `ConstraintSet` | `constraints/` | §6 |
| 求解 | `PathDerivatives`+`ConstraintSet` → `Profile(a,b,c)` | `solve/`（长序列外套 `hlaw/`） | §7 / §8 |
| 合成 | `Profile` → `TrajectoryResult`（q(t),q̇,q̈,q⃛） | `synth/` | §9 |

---

## 4. 核心数据类型

`copp/types.py`（全部为 `@dataclass`，字段视图；无实现）：

```python
from dataclasses import dataclass
from typing import Literal, Protocol
import numpy as np

Pose = np.ndarray            # (7,) = [x,y,z, qw,qx,qy,qz]  位置+单位四元数
NativeSpace = Literal["joint", "cartesian"]

@dataclass
class PoseSegment:
    """单条指令在局部参数 u∈[0,1] 上的笛卡尔位姿段（设计 §3.4）。"""
    pose_of_u:  "Callable[[np.ndarray], np.ndarray]"   # u → (7,·) 位姿求值器
    dpose_of_u: "Callable[[np.ndarray, int], np.ndarray]"  # (u, 阶数1~3) → 导数
    length_hint: float                                  # 弧长/角度估计
    native_space: NativeSpace
    kin_branch: "IKBranch | None"                       # IK 分支/构型标记
    constraint_overrides: "ConstraintOverrides | None"  # 该段专属限值（可选）

@dataclass
class BlendedPath:
    """截短原始段 + TransP5 过渡段交替的 C²/分段C³ 位姿路径（设计 §4.6）。"""
    r_of_s:  "Callable[[np.ndarray], np.ndarray]"       # s → 位姿
    dr_of_s: "Callable[[np.ndarray, int], np.ndarray]"  # (s, 阶数1~3) → 笛卡尔导数
    s_breaks: np.ndarray                                # 段边界（含约束域边界）
    s_total: float
    high_curvature_marks: np.ndarray                    # 过渡段高曲率标记（供加密）

@dataclass
class PathDerivatives:
    """关节路径导数，喂给 solve/（对齐 Rust copp Path，设计 §5.5）。"""
    s_grid: np.ndarray          # (N,)
    q:    np.ndarray            # (n, N)
    dq:   np.ndarray            # (n, N)   q'(s)
    ddq:  np.ndarray            # (n, N)   q''(s)
    dddq: np.ndarray            # (n, N)   q'''(s)
    singular: np.ndarray        # (N,) bool
    seg_index: np.ndarray       # (N,) 每站点所属指令/约束段
    cart: "CartCoeffs | None"   # 预计算 ||p'||², ||J_ω q'||²（TCP 速度约束系数）

@dataclass
class Profile:
    """求解结果 (a,b,c) 剖面（设计 §7；3 阶方法）。"""
    a: np.ndarray               # (N,)  ṡ²
    b: np.ndarray               # (N,)  s̈
    c: np.ndarray | None        # (N,)  s⃛ṡ（分段常值）
    num_stationary: tuple[int, int]

@dataclass
class TrajectoryResult:
    """最终轨迹（设计 §9）。"""
    s_grid: np.ndarray
    profile: Profile
    t_s: np.ndarray             # 到达每个 s_k 的时间
    t_final: float
    t: np.ndarray               # 均匀时间栅格
    q: np.ndarray; qd: np.ndarray; qdd: np.ndarray; qddd: np.ndarray   # (n, T)
    metrics: "VerifyMetrics | None"
```

---

## 5. 各子包算法原理与接口骨架

> 每个模块给出：**职责 / 算法原理（引设计文档与论文）/ 关键接口签名骨架**。所有函数体以 `...` 占位。

### 5.1 robot/ 运动学-动力学后端

**职责**：把机器人本体（FK/IK/Jacobian/逆动力学）经协议解耦，供 blending（FK 桥接）、lowering（IK + Jacobian 链）、constraints（逆动力学）注入。

**算法原理**：纯适配层，不含规划算法。IK 要求**支持 seed**（返回离 seed 最近的解分支）以保证 lowering 的连续解选择（设计 §5.2）。`jacobian_derivative` 用沿路径方向的有限差分（设计 §5.3 的 $\mathbf J',\mathbf J''$）。

```python
# robot/base.py
class KinematicsModel(Protocol):
    dof: int
    def fk(self, q: np.ndarray) -> Pose: ...
    def ik(self, pose: Pose, seed: np.ndarray) -> np.ndarray: ...        # 返回离 seed 最近解
    def jacobian(self, q: np.ndarray) -> np.ndarray: ...                 # (6, n) 几何 Jacobian
    def jacobian_omega(self, q: np.ndarray) -> np.ndarray: ...           # (3, n) 角速度分量 J_ω
    def cond(self, q: np.ndarray) -> float: ...                         # 条件数（奇异检测）

class DynamicsModel(Protocol):                                          # 可选（力矩约束/热能目标）
    def inverse_dynamics(self, q, dq, ddq) -> np.ndarray: ...           # RNE → 关节力矩
    def inertia(self, q) -> np.ndarray: ...

# robot/synthetic.py   : class SyntheticRobotModel                     M1 stand-in（解析 TCP + 对角惯性）
# robot/dh_poe.py      : class DhPoeKinematics(KinematicsModel)        自实现最小 DH/POE
# robot/rtb_adapter.py : class RtbKinematics(KinematicsModel, DynamicsModel)  包装 roboticstoolbox
```

### 5.2 commands/ 指令层

**职责**：三类运动指令 → 统一 `PoseSegment`（笛卡尔位姿路径段，含解析导数）。

**算法原理**（设计 §3）：
- `JointMoveCommand`：关节空间五次多项式 $q(u)$（两端 $\dot q,\ddot q$ 可指定，默认 stationary）；位姿表示 $r(u)=\mathrm{FK}(q(u))$（仅过渡区需要，`native_space="joint"` 标记走快路径）。
- `LinearMoveCommand`：位置线性插值 + 姿态 `Slerp`；导数解析。
- `CircularMoveCommand`：三点/圆心法定圆 → $p(\theta)=c+r(\cos\theta\,e_1+\sin\theta\,e_2)$，$\theta(u)$ 线性；姿态 `Slerp`；导数解析（三角函数）。退化情形抛错。

```python
# commands/base.py
class MotionCommand(Protocol):
    def to_segment(self, kin: KinematicsModel) -> PoseSegment: ...

@dataclass
class JointMoveCommand(MotionCommand):
    q_start: np.ndarray; q_end: np.ndarray
    boundary: Literal["stationary", "continuous"] = "stationary"
    def to_segment(self, kin): ...        # 五次多项式 q(u) → PoseSegment(native="joint")

@dataclass
class LinearMoveCommand(MotionCommand):
    pose_start: Pose; pose_end: Pose
    def to_segment(self, kin): ...        # 位置直线 + Slerp → PoseSegment(native="cartesian")

@dataclass
class CircularMoveCommand(MotionCommand):
    pose_start: Pose; pose_end: Pose
    aux: "Pose | tuple[np.ndarray, np.ndarray]"    # via 点 或 (center, normal)
    direction: Literal["shortest", "ccw", "cw"] = "shortest"
    def to_segment(self, kin): ...        # 圆弧几何 + Slerp → PoseSegment(native="cartesian")
```

### 5.3 blending/ 最优 Hermite 过渡

**职责**：相邻 `PoseSegment` 交接处在笛卡尔空间做 G2 五次 Hermite 过渡，装配出 C²/分段C³ 的 `BlendedPath`。

**算法原理**（设计 §4）：
- `zone.py`：$\Delta_j=\min(\Delta_\text{max},L_1/3,L_2/3)$，两侧各截 $\Delta_j$。
- `frenet.py`：在截短端点提取 Frenet 框架 $(\mathbf t,\mathbf n,\kappa)$。
- `optimal_hermite.py`：解 G2 五次 Hermite，自由参数 $\alpha_0,\alpha_1>0$ 按两端曲率是否为零分四情形（线性/三次/九次），多根按**最小跃度代价** $\min\int\|\mathbf p_5'''\|^2$ 选；复用 `1.OptimalHermiteInterpolation` 已有算法。
- `pose_blend.py`：位置+姿态同纳 Frenet，量纲归一化 $D$；姿态在四元数切空间处理。
- `junction.py`：混合指令桥接——Joint 段截短端用 FK 求位姿导数后与 Cartesian 段过渡；Joint–Joint 可走关节空间快路径。

```python
# blending/zone.py
def zone_distance(L_prev: float, L_next: float, delta_max: float) -> float: ...
def trim_segment(seg: PoseSegment, delta: float, side: Literal["head","tail"]) -> PoseSegment: ...

# blending/frenet.py
def frenet_frame(pose, dpose_1, dpose_2) -> "FrenetFrame":  ...   # (t, n, κ)

# blending/optimal_hermite.py
def solve_g2_quintic(left: "FrenetFrame", right: "FrenetFrame",
                     dim_weights: np.ndarray) -> "QuinticTransition":
    """G2 五次 Hermite：解 α0,α1（分四情形），最小跃度选根。返回多项式系数。"""
    ...
def jerk_cost(coeffs: np.ndarray) -> float: ...   # ∫||p'''||² 代价（选根用）

# blending/pose_blend.py
def blend_pair(seg_prev, seg_next, opts: "BlendOptions", kin) -> "list[PoseSegment]":
    """截短两段 + 插入 TransP5，返回 [prev_cut, trans, next_cut]（失败则退回 G1 角点）。"""
    ...

# blending/junction.py
def assemble(segments: list[PoseSegment], opts, kin) -> BlendedPath:
    """对每对相邻段调 blend_pair，拼成 C²/分段C³ 的 BlendedPath，记录 s_breaks。"""
    ...
```

### 5.4 lowering/ 降维到关节空间

**职责**：把笛卡尔 `BlendedPath` 采样并 IK 降维为关节 `PathDerivatives`。

**算法原理**（设计 §5，复用 `robot6dof_topp_design.md` §2–3）：
- `sampling.py`：曲率驱动步长 $\Delta s_m=\min(\Delta s_\text{max},\sqrt{8\varepsilon_\text{pos}/\kappa_\text{pos}},\sqrt{8\varepsilon_\text{ori}/\kappa_\text{ori}})$；**网格强制含所有 `s_breaks`（约束域边界）**——保证论文 Theorem 1 的 $O(\Delta^2)$ 误差界成立。
- `ik.py`：逐点 `kin.ik(pose, seed=q_prev)` 连续解；检测 $\|\Delta q\|$ 跳变。
- `derivatives.py`：$\mathbf q'=\mathbf J^{-1}\mathbf r'$，$\mathbf q''=\mathbf J^{-1}\mathbf r''-\mathbf J^{-1}\mathbf J'\mathbf q'$，$\mathbf q'''=\mathbf J^{-1}\mathbf r'''-\mathbf J^{-1}\mathbf J'\mathbf q''-\mathbf J^{-1}\mathbf J''\mathbf q'$；$\mathbf J',\mathbf J''$ 前向/中心差分。同时预计算 TCP 约束系数 $\|\mathbf p'\|^2$、$\|\mathbf J_\omega\mathbf q'\|^2$。
- `singularity.py`：$\mathrm{cond}(\mathbf J)$ 超阈值 → DLS 逆 + 局部加密。

```python
# lowering/sampling.py
def adaptive_sample(path: BlendedPath, opts: "SampleOptions") -> np.ndarray:
    """曲率驱动采样，返回 s_grid（含 path.s_breaks）。"""
    ...

# lowering/ik.py
def solve_ik_sequence(path, s_grid, kin, seed0) -> tuple[np.ndarray, np.ndarray]:
    """逐点连续解 IK，返回 q[n,N] 与 singular[N]。"""
    ...

# lowering/derivatives.py
def joint_derivatives(path, s_grid, q, kin) -> PathDerivatives:
    """链式法则求 q',q'',q'''（含 J',J'' 差分）+ 预计算 TCP 系数。"""
    ...

# lowering/singularity.py
def damped_inverse(J: np.ndarray, lam: float = 0.05) -> np.ndarray: ...
def refine_local(s_grid, idx, factor: int = 4) -> np.ndarray: ...
```

### 5.5 constraints/ 约束摄入

**职责**：物理限值 → 站点索引的路径域 1/2/3 阶不等式（对齐 copp `Robot`/`Constraints`）。

**算法原理**（设计 §6，论文附录 A.2）：
- 轴向速度/加速度/jerk → eq.40；力矩/力矩率 → eq.43–45（需 `DynamicsModel`）。
- **TCP 仅两项速度模长**（设计 v0.3）：位置速度模长 $\|\mathbf p'\|^2 a\le v_\text{tcp,max}^2$、姿态角速度模长 $\|\mathbf J_\omega\mathbf q'\|^2 a\le\omega_\text{tcp,max}^2$——均为 $a$ 的**线性上界**，逐点取 min 合并进 $\bar a(s)$。
- 支持**非对称、参数变、分段、跨坐标系**约束（论文 Assumption 2）；`ConstraintSet` 按站点区间存不同上下界。

```python
# constraints/model.py
@dataclass
class ConstraintSet:
    a_max: np.ndarray               # (N,)   a 上界（含轴向速度+TCP 速度模长，取 min 后）
    acc_n: np.ndarray; acc_m: np.ndarray; acc_g: np.ndarray   # 2 阶：n·a+m·b≤g（逐轴）
    jerk_r: ...; jerk_s: ...; jerk_t: ...; jerk_h: ...; jerk_f: ...   # 3 阶：r·a+s·b+t·c+h ≤ f·a^{-1/2}
    torque: "TorqueRows | None"     # 力矩（可选）
    # 允许每站点不同 → 参数变/非对称/分段

# constraints/ingest.py
def ingest(pd: PathDerivatives, limits: "RobotLimits",
           dyn: DynamicsModel | None) -> ConstraintSet:
    """把物理限值映射为路径域系数矩阵；TCP 速度模长并入 a_max。"""
    ...
```

### 5.6 solve/ copp 求解层（SPLP 核心）

**职责**：在 `PathDerivatives`+`ConstraintSet` 上按论文 **TOTP-SPLP** 求时间最优 `Profile(a,b,c)`。**默认路线 §7.2.1(a)：PLP 分段线性目标 + LP + Algorithm 2 迭代。**

**算法原理**（设计 §7 / [`paper_notes.md`](../0.other_lib_code/copp/copp/docs/paper_notes.md) §4–§7）：

**`state.py` — 无损离散化**
- 状态 $a=\dot s^2,b=\ddot s,c=\dddot s\dot s$；系统动力学离散 $a_k=a_{k-1}+2b_{k-1}\Delta_k+c_k\Delta_k^2,\ b_k=b_{k-1}+c_k\Delta_k$（Prop.1）。
- 静止段（$a_s/a_f=0$）用 $\dddot s$ 常值模型：$a_k\propto(\Delta s)^{4/3}$（Prop.2 / eq.20）。返回等式约束系数与静止段长度。

```python
def build_dynamics_rows(s_grid, num_stationary) -> "LinearSystem": ...   # a'=2b,b'=c 的离散等式
def stationary_profile(s_grid, a_stationary, num_stationary): ...         # eq.20
```

**`seed.py` — 种子 a⁽⁰⁾（= topp2_ra）**
- 2 阶问题（忽略 3 阶约束）求 $a^{(0)}$：可用 DP 可达集（前/后向）或 2 阶 LP $\max\sum(u_{k+1}-u_{k-1})a_k$。只需保证首次线性化 LP 可行（论文 §5.1 / §7.3）。

```python
def compute_seed(pd, cons: ConstraintSet, boundary) -> np.ndarray: ...    # → a0[N]
```

**`linearize.py` — jerk 凹约束线性化（eq.32）**
- 在参考 $a^{(p-1)}$ 处对 $a^{-1/2}$ 取切线，凹约束变仿射：
  $(r+\tfrac{f}{2a_{\text{lin}}^{3/2}})a+s\,b+t\,c\le\tfrac{3f}{2\sqrt{a_\text{lin}}}-h$。等价于 copp `build_with_linearization`。

```python
def linearize_jerk(cons: ConstraintSet, a_lin: np.ndarray, floor=1e-10) -> "AffineJerkRows": ...
```

**`plp_objective.py` — PLP 分段线性目标（本方案核心，eq.27）**
- 每个 $a_k$ 选采样点 $0<\delta_{k,0}<\dots<\delta_{k,P_k}$，用割线上包络逼近 $1/\sqrt{a_k}$；引入辅助变量 $J_k$ 与约束（eq.29d）
  $J_k\ge\dfrac{\delta_{k,i-1}+\sqrt{\delta_{k,i-1}\delta_{k,i}}+\delta_{k,i}-a_k}{(\sqrt{\delta_{k,i-1}}+\sqrt{\delta_{k,i}})\sqrt{\delta_{k,i-1}\delta_{k,i}}}$；
  目标 $\min\sum(u_{k+1}-u_{k-1})J_k$；再加下界 $a_k\ge\delta_{k,0}$（Prop.3，根除零进给奇异）。

```python
def build_plp_objective(s_grid, num_stationary,
                        deltas: "DeltaSamples") -> "PlpObjective":
    """返回 J_k 辅助变量声明、eq.29d 割线约束行、目标权重、a_k≥δ0 下界。"""
    ...
def default_delta_samples(a_seed: np.ndarray, levels=(1e-4,1e-3,1e-2,1e-1)) -> "DeltaSamples":
    """δ_{k,l}=10^{l-4}·a_seed[k]（对齐论文实验设置）。"""
    ...
```

**`lp_problem.py` — 单次 LP 组装（cvxpy）**
- 决策变量 $x=[a,b,c,J]$；等式=离散动力学（state.py）+ 边界；不等式=速度上界 $a\le\bar a$、加速度箱式、线性化 jerk（linearize.py）、PLP 割线（plp_objective.py）、下界 $a\ge\delta_0$；目标=PLP。经 `backend` 求解。

```python
def build_and_solve_lp(pd, cons, dyn_rows, jerk_rows, plp, boundary,
                       backend: "SolverBackend") -> Profile:
    """组装一次 LP（对应论文式 29）并求解，返回 (a,b,c)。"""
    ...
```

**`splp.py` — Algorithm 2 迭代循环**
- 从 $a^{(0)}$ 起，第 $p$ 次：linearize_jerk(a^{(p-1)}) → build_and_solve_lp → a^{(p)}；停止准则 eq.30。

```python
def solve_splp(pd, cons, boundary, opts: "SolveOptions",
               backend, seed: np.ndarray | None = None) -> Profile:
    """论文 TOTP-SPLP（Algorithm 2）。opts.n_iter 默认 2~3；opts.mode∈{plp_lp, socp}。"""
    ...
```

**`interp.py` — 解析插值（Prop.1/2）**
- 闭式 $\Phi_k,\Phi_k^{-1}$（按 $c_k>0/<0/=0$ 分支）；`s_to_t` 算区间时长，`t_to_s` 时间域采样。**无插值误差**。

```python
def s_to_t(s_grid, profile: Profile) -> tuple[float, np.ndarray]: ...     # → (t_final, t_s)
def t_to_s(s_grid, profile, t_s, dt: float) -> np.ndarray: ...            # 均匀时间栅格 → s(t)
```

> **备选路线 (b) SOCP**：`opts.mode="socp"` 时，`lp_problem.py` 改用 SOCP 建模——引入 $\eta_k=1/\sqrt{a_k}$ 二阶锥、目标 $\min\sum w_k\eta_k$（精确时间），其余（jerk 线性化、动力学、约束、Algorithm 2 循环）不变。对应设计 §7.2.1(b)。

### 5.7 hlaw/ 分层前瞻窗口

**职责**：长指令序列/流式规划时，分窗调度并**保证每窗可行**（论文 §5 HLAW）。

**算法原理**（设计 §8）：三窗（种子/可行/最优）层级前移，**可行窗与最优窗都在同一种子窗的 $a^{(0)}$ 上线性化**——线性化点跨窗一致 → Theorem 3 可证可行。种子窗跑 `seed.compute_seed`；可行窗跑 `solve_splp(n_iter=1)`；最优窗跑 `solve_splp(n_iter≥2)`。

```python
# hlaw/windows.py
def plan_windowed(pd, cons, boundary, win: "WindowSpec",
                  opts, backend) -> Profile:
    """三窗调度：seed → feasibility(p=1) → optimality(p≥2)，逐层 level-by-level。"""
    ...

# hlaw/relay.py
def relay_boundary(profile_prev: Profile, k_left, k_right) -> "BoundaryCond":
    """从上一窗解提取下一窗 (a,b) 边界（非静止衔接）。"""
    ...
```

### 5.8 synth/ 轨迹合成与验证

**职责**：`Profile` → 时间域轨迹 + 约束满足度校验。

**算法原理**（设计 §9）：解析插值采样（interp.py）后
$\dot q=q'\sqrt a$，$\ddot q=q''a+q'b$，$\dddot q=q'''a^{3/2}+3q''\sqrt a\,b+q'c$；校验超限率 $R_v$、超限时长比 $D_v$（论文 §6.1.2）。

```python
# synth/resample.py
def synthesize(pd, profile, dt: float) -> TrajectoryResult: ...

# synth/verify.py
def verify_limits(result: TrajectoryResult, limits) -> "VerifyMetrics": ...   # R_v / D_v
```

### 5.9 backend/ 求解器后端抽象

**职责**：把"组装好的凸问题 → 数值解"与建模解耦，便于 LP/SOCP 切换与对照现成 copp。

**算法原理**：`cvxpy_backend` 默认（LP 走 clarabel/glpk/highs，SOCP 走 clarabel）；`copp_backend` 可选——直接调用现成 copp 的 `topp2_ra`/`build_with_linearization`/`topp3_lp`/`topp3_socp`（Rust 经 pyo3 或 copp_py），复用其离散化与 eq.32 线性化以获得更高性能与数值一致性。

```python
# backend/base.py
class SolverBackend(Protocol):
    def solve_lp(self, problem: "LpSpec") -> "Solution": ...
    def solve_socp(self, problem: "SocpSpec") -> "Solution": ...

# backend/cvxpy_backend.py : class CvxpyBackend(SolverBackend)     默认
# backend/copp_backend.py  : class CoppBackend(SolverBackend)      可选，委托现成 copp
```

### 5.10 planner.py 门面 + 全局设施

**职责**：串起全流程的用户入口；`options/errors/diagnostics` 为横切设施。

```python
# planner.py
class TrajectoryPlanner:
    def __init__(self, kin: KinematicsModel, dyn: DynamicsModel | None = None,
                 backend: SolverBackend | None = None): ...
    def add_command(self, cmd: MotionCommand) -> "TrajectoryPlanner": ...     # 链式累积
    def set_limits(self, **limits) -> "TrajectoryPlanner": ...                # v/a/j/τ/tcp
    def plan(self, opts: "PlannerOptions | None" = None) -> TrajectoryResult:
        """编排：commands→blending→lowering→constraints→(hlaw?)solve→synth。"""
        ...

# options.py     : BlendOptions / SampleOptions / SolveOptions(mode,n_iter,eps_t,eps_a) / WindowSpec / PlannerOptions
# errors.py      : RobotCoppError ← {DegenerateArc, KinematicSingularity, BlendFailed, InfeasibleWindow, SolverStatus}
# diagnostics.py : Verbosity(Silent/Summary/Debug/Trace) + 计时日志（对齐 copp Verboser 风格）
```

---

## 6. 关键算法伪代码

> 伪代码描述控制流，非实现。

**(A) 全流程编排（`TrajectoryPlanner.plan`）**
```
segments = [cmd.to_segment(kin) for cmd in commands]           # §5.2
path     = blending.assemble(segments, blend_opts, kin)        # §5.3  → C²/分段C³ r(s)
s_grid   = lowering.adaptive_sample(path, sample_opts)         # §5.4  含 path.s_breaks
q, sing  = lowering.solve_ik_sequence(path, s_grid, kin, seed0)
pd       = lowering.joint_derivatives(path, s_grid, q, kin)
cons     = constraints.ingest(pd, limits, dyn)                 # §5.5
if long_sequence:  profile = hlaw.plan_windowed(pd, cons, bnd, win, solve_opts, backend)   # §5.7
else:              profile = solve.solve_splp(pd, cons, bnd, solve_opts, backend)          # §5.6
return synth.synthesize(pd, profile, dt)                       # §5.8
```

**(B) SPLP 迭代（`solve.solve_splp`，论文 Algorithm 2）**
```
a_prev = seed or seed.compute_seed(pd, cons, bnd)              # a⁽⁰⁾
dyn_rows = state.build_dynamics_rows(pd.s_grid, num_stat)
plp      = plp_objective.build_plp_objective(pd.s_grid, num_stat,
                    plp_objective.default_delta_samples(a_prev))
for p in 1..n_iter:
    jerk_rows = linearize.linearize_jerk(cons, a_prev)         # eq.32 在 a_prev 处
    profile   = lp_problem.build_and_solve_lp(pd, cons, dyn_rows,
                    jerk_rows, plp, bnd, backend)               # 一次 LP（式29）
    if |t_f(profile) - t_f_prev| < eps_t or ||a - a_prev|| < eps_a:  break   # eq.30
    a_prev = profile.a
return profile
```

**(C) HLAW 三窗（`hlaw.plan_windowed`）**
```
for each window j (level-by-level, 单向依赖):
    seed窗:   a0_j   = seed.compute_seed(window_j)                     # 2 阶
    可行窗:   feas_j = solve_splp(window_{j}, seed=a0_j, n_iter=1)     # p=1，在 a0_j 线性化
    最优窗:   opt_j  = solve_splp(window_{j}, seed=a0_j, n_iter≥2)     # 仍在 a0_j 线性化 → 跨窗一致
    relay.relay_boundary(opt_j, ...) → 下一窗边界
拼接各窗 opt_j（重叠区取后窗）→ 全局 Profile
```

**(D) 最优 Hermite 过渡（`blending.blend_pair`）**
```
Δ = zone.zone_distance(L_prev, L_next, delta_max)
prev_cut = zone.trim_segment(prev, Δ, "tail");  next_cut = zone.trim_segment(next, Δ, "head")
fL = frenet.frenet_frame(prev_cut@u=1);  fR = frenet.frenet_frame(next_cut@u=0)
trans = optimal_hermite.solve_g2_quintic(fL, fR, dim_weights)   # α0,α1 分四情形，min jerk_cost 选根
return [prev_cut, trans, next_cut]  (若 solve 失败 → 保留 G1 角点)
```

---

## 7. 端到端调用示例

> API 使用示意（非实现）。

```python
from copp import TrajectoryPlanner
from copp.robot import RtbKinematics
from copp.commands import JointMoveCommand, LinearMoveCommand, CircularMoveCommand
from copp.options import PlannerOptions, SolveOptions

kin = RtbKinematics.from_urdf("ur5.urdf")
planner = (TrajectoryPlanner(kin)
    .add_command(JointMoveCommand(q_start=q0, q_end=q1))
    .add_command(LinearMoveCommand(pose_start=p1, pose_end=p2))
    .add_command(CircularMoveCommand(pose_start=p2, pose_end=p3, aux=p_via))
    .set_limits(qd_max=..., qdd_max=..., qddd_max=..., v_tcp_max=..., w_tcp_max=...))

result = planner.plan(PlannerOptions(
    blend=dict(delta_max=0.01),
    solve=SolveOptions(mode="plp_lp", n_iter=3),     # 默认 PLP+LP（算力最省）
    dt=1e-3,
))
# result.q / result.qd / result.qdd / result.qddd, result.t_final, result.metrics
```

---

## 8. 测试策略

| 层 | 内容 | 依赖 |
|----|------|------|
| unit | 各模块隔离：blending（G2 连续性/最小跃度）、lowering（q',q'',q''' 对解析解）、solve（单次 LP 可行性、PLP 割线上包络单调性）、interp（$\Phi_k\circ\Phi_k^{-1}=\mathrm{id}$） | pytest |
| integration | 端到端：joint-only（无需 IK）、mixed（joint+line+arc）、objectives（time vs 备选 SOCP） | pytest |
| benchmark | 随机路径批量（对照论文 $R_v,D_v<0.1\%$、$t_f$ 单调不增），`@pytest.mark.slow` | pytest + hypothesis |

关键**性质断言**（对齐论文）：SPLP 迭代 $t_f^{(p)}$ 单调不增；约束超限 $R_v,D_v<0.1\%$；HLAW 全窗 0 不可行；解析插值零误差。

---

## 9. 实现里程碑建议

按依赖顺序、可独立验证的粒度推进：

1. **M1 数值内核**：`types` + `solve/`（state / seed / linearize / plp_objective / lp_problem / splp / interp）+ `backend/cvxpy_backend`。用**合成解析路径**（如正弦 $q(s)$，绕开 IK）跑通 SPLP，对照 copp Rust `topp3_lp` 测试数据。← 先落地论文 TOTP-SPLP。
2. **M2 指令+降维**：`robot/`（先 DH-POE）+ `commands/` + `lowering/`。joint-only 与 line/arc（无 blending）端到端。
3. **M3 blending**：`blending/`（可先移植 `1.OptimalHermiteInterpolation` 算法）+ `junction` 混合桥接。
4. **M4 约束扩展**：`constraints/` 力矩（接 `DynamicsModel`）+ TCP 速度模长；备选 SOCP 路线。
5. **M5 HLAW**：`hlaw/` 三窗调度（长序列/流式）——设计文档 §12.2 标注的主要新增工作量之一。
6. **M6 性能**：`backend/copp_backend`（委托 Rust copp）+ 直接稀疏组装替代 cvxpy 热路径。

---

*文档版本：v0.1（Python 框架）｜配套设计文档 robot_copp_design.md v0.5｜仅框架，无实现*
