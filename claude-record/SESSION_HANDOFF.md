# Session Handoff — COPP-Python 轨迹规划框架设计

> 导出时间：2026-07-10
> 来源仓库：`e:\CodePrj\CppPri-VScode\0.other_lib_code\copp\copp`（Rust `copp` 库，分支 `feature/wc/add_save_data_plot_script`）
> 用途：把本次会话的目标、决策与产出打包成一份自包含文档，方便在**另一个项目**中新开会话时粘贴/上传作为上下文，继续后续工作（例如：把设计落地为实际 Python 代码）。

---

## 1. 背景：Rust `copp` 库是什么

`copp`（Convex-Objective Path Parameterization）是一个 Rust 轨迹规划库，解决：给定几何路径 `q(s)`，在速度/加速度/急动度/力矩约束下，求满足约束且按凸目标最优的时间律 `s(t)`。

求解器家族：

| 阶数 | 时间最优 | 一般凸目标 |
|------|----------|------------|
| 2 阶（速度/加速度/力矩） | TOPP2-RA | COPP2-SOCP |
| 3 阶（+ 急动度） | TOPP3-LP / TOPP3-SOCP | COPP3-SOCP |

Rust 侧代码在 `src/`（`path/`、`robot/`、`copp/{copp2,copp3}`、`math/`、`diag/`），公开 API 见 `src/lib.rs`。项目仓库根目录已有一份详细的中文代码阅读指南：`docs/code_reading_guide.md`（约 700 行，覆盖每个模块的数学公式、调用链、约束映射公式、Clarabel 后端封装等）。**如果继续在原仓库工作，这份文件是最权威的算法参考**；如果换到新项目且没有原仓库访问权限，下面第 3 节已经把关键内容摘录进本文档。

---

## 2. 本次会话的任务与结论

### 2.1 原始需求

> 根据 `code_reading_guide.md` 方案，设计一个 Python 软件框架，包含环境依赖（数学库、LP 求解库、机器人运动学动力学库等）、基于 robot 的 COPP 算法，指令为关节运动指令、直线指令、圆弧指令，仅需要设计出轨迹规划模块及相关测试用例。

### 2.2 通过 AskUserQuestion 确认的范围（重要，决定了产出形态）

- **交付形式**：仅设计文档，不写代码。
- **算法实现深度的描述口径**：按 `code_reading_guide.md` 的算法描述来（即文档中要求 Python 复刻要与 Rust 版本算法逻辑 1:1 对应，而不是简化版）。
- **代码位置**：新建 `python/` 子目录（与 Rust crate 同仓库并存）。

### 2.3 产出

在仓库中新建了两个文件（**设计文档，无可运行代码**）：

- `python/DESIGN.md` —— 完整架构设计文档（第 3 节完整复制于此）
- `python/README.md` —— 指向 `DESIGN.md` 的简短索引

这两个文件目前**还没有被 git commit**（用户没有要求提交，遵循"仅在用户明确要求时才 commit"的规则）。如果继续在原仓库工作，记得检查 `git status` 确认这两个新文件的状态。

---

## 3. 设计文档全文（`python/DESIGN.md`）

> 以下为该文件在导出时的完整内容，供在没有原仓库访问权限的新项目中直接使用。

<!-- BEGIN DESIGN.md -->

# COPP-Python 轨迹规划框架设计文档

> 定位：以 `copp`（Rust）库的算法体系为基础，面向机器人应用层设计的 Python 轨迹规划框架。
> 输入为**指令序列**（关节运动 / 直线 / 圆弧），输出为满足运动学约束、按凸目标最优的时间律轨迹。
> 本文档仅描述架构设计，**不包含可运行代码**；核心求解算法的数学形式与 `docs/code_reading_guide.md` 保持一致。

## 目录

1. 设计目标与范围
2. 环境依赖
3. 总体架构
4. 指令体系（commands/）
5. 轨迹规划模块（planning/）——核心设计
6. COPP 求解器复刻（copp/）
7. 求解器后端抽象（solver_backend/）
8. 端到端示例（伪代码）
9. 测试用例设计（tests/）
10. 与 Rust 版本的关系与演进路径
11. 待决问题 / 后续扩展点

## 1. 设计目标与范围

Rust 版 `copp` 只解决"给定几何路径 `q(s)`，求时间最优/凸目标最优的 `s(t)`"这一层问题，路径本身如何生成不在其范围内。

本框架在此基础上新增**指令层**（Rust 版没有），把工业机器人/CNC 常见的三类运动指令翻译为几何路径，再复用 COPP 的约束摄入与求解流程：

```
指令序列（关节运动 / 直线 / 圆弧）
        │  PathBuilder：翻译 + 拼接 + 重参数化
        ▼
   几何路径 q(s), s ∈ [0, s_f]
        │  Robot：物理约束 → 路径域不等式
        ▼
   COPP 求解（TOPP2-RA / COPP2-SOCP / TOPP3-LP / TOPP3-SOCP / COPP3-SOCP）
        │  s_to_t / t_to_s 后处理
        ▼
   时间域轨迹 q(t), q̇(t), q̈(t)
```

**范围声明**：本次设计仅覆盖 **轨迹规划模块（`planning/`）及其直接依赖（`commands/`、`copp/`、`robot/`、`path/`）与相关测试用例**。不涉及底层伺服通信、可视化 UI、在线重规划调度等应用层功能，这些留作后续扩展点（见第 11 节）。

## 2. 环境依赖

### 2.1 依赖分层

| 层次 | 用途 | 候选库 | 选型建议 |
|------|------|--------|----------|
| 数值基础 | 向量化数值计算 | `numpy` | 必选，全项目的数据交换格式（`ndarray`） |
| 数值基础 | 稀疏矩阵/数值积分/样条 | `scipy` | 必选；`scipy.interpolate`（Hermite/B样条）、`scipy.optimize.linprog`（LP 兜底）、`scipy.spatial.transform.Rotation`（姿态插值） |
| 自动微分 | 路径三阶导数（对应 Rust `Jet3`） | 自实现前向 AD 标量类 / `jax.grad`+`jax.jacfwd` | **优先自实现**轻量三阶前向 AD 类（见 4.1 节 `autodiff.py` 设计），避免引入 `jax` 的安装体积与 GPU 依赖；`jax` 作为可选加速后端 |
| 凸优化建模 | SOCP/LP/QP 统一建模 | `cvxpy` | 必选；对应 Rust 中手工组装 Clarabel 稀疏矩阵的角色，`cvxpy` 负责把锥约束/目标翻译为标准形式 |
| 凸优化求解器 | SOCP/LP 后端求解 | `clarabel`（PyPI 有官方 Python 绑定）、`ecos`、`scs`、`osqp` | 优先 `clarabel`（与 Rust 版本算法行为一致，便于数值对照）；`ecos`/`scs` 作为 `cvxpy` 默认可用的兜底后端 |
| 2 阶 DP 可达集 | TOPP2-RA 的 2D LP 增量求解 | 自实现 Seidel 增量 2D LP（对应 `math/numerical/lp.rs`）；兜底用 `scipy.optimize.linprog` | 自实现以保证与 Rust 版相同的 `O(m)` 期望复杂度；`linprog` 仅用于单元测试中的交叉验证 |
| 机器人运动学 | 正逆运动学（直线/圆弧指令需要 IK） | `roboticstoolbox-python`（Peter Corke）、`ikpy`、自实现 DH/POE | 设计上通过 `KinematicsModel` 协议解耦，默认提供自实现 DH/POE 最小实现，`roboticstoolbox-python` 作为可插拔适配器（支持解析式/数值 IK、雅可比） |
| 机器人动力学 | 力矩约束、热能目标（`with_axial_torque`、`ThermalEnergy`） | `roboticstoolbox-python`（内建 RNE 逆动力学）、`pinocchio`（若可用） | 通过 `DynamicsModel` 协议解耦；`pinocchio` 性能更优但在 Windows 上安装成本高，默认建议 `roboticstoolbox-python`，`pinocchio` 作为可选后端 |
| 姿态与几何 | 四元数/旋转插值、圆弧几何 | `scipy.spatial.transform`（`Rotation`, `Slerp`） | 必选，避免自行实现四元数 SLERP 数值细节 |
| 测试 | 单元/集成测试 | `pytest` | 必选 |
| 属性测试 | 随机样条路径批量验证（对应 Rust `tests/test_random_spline.rs`） | `hypothesis` | 建议，用于生成随机约束/随机样条组合，捕获边界情况 |
| 可视化（供人工核查，非规划模块本体） | 轨迹曲线绘制 | `matplotlib` | 建议，复用现有 `scripts/plot_*.py` 的风格 |
| 打包 | 依赖与元数据管理 | `pyproject.toml` + `hatchling`（或 `setuptools`） | 建议用 `uv`/`pip` 管理虚拟环境 |

### 2.2 版本与约束建议

- Python `>=3.10`（使用 `match` 语句表达指令类型分发、`typing` 的 `Literal`/`Protocol` 描述接口）。
- `numpy>=1.26`、`scipy>=1.11`、`cvxpy>=1.4`、`clarabel>=0.6`（Python 绑定）。
- 可选依赖以 extras 形式声明，例如：

  ```toml
  [project.optional-dependencies]
  kinematics = ["roboticstoolbox-python"]
  jax = ["jax"]
  viz = ["matplotlib"]
  test = ["pytest", "hypothesis"]
  ```

  这样核心 `planning/` 模块在只安装最小依赖集时即可完成"关节运动指令 + TOPP2-RA/COPP2-SOCP"这类不需要 IK/动力学的场景。

### 2.3 为什么不直接绑定 Rust 求解器？

技术上可行（`pyo3` + `maturin` 生成 Python 扩展），且是长期最优路径（见第 10 节），但作为**框架设计的第一版**，先用纯 Python + `cvxpy`/`clarabel-python` 复刻算法，原因：

- 便于阅读、调试、教学，不需要为贡献者引入 Rust 工具链；
- `clarabel` 的 Python 绑定与 Rust crate 是同一核心求解器，数值行为一致，性能损失主要来自 Python 层的约束矩阵组装（可接受，规划模块非高频热路径）；
- 保留一个 `solver_backend/rust_ffi_backend.py` 扩展点，性能敏感场景可直接切换到 Rust 后端（见第 7 节）。

## 3. 总体架构

```
python/
├── pyproject.toml
├── README.md
├── docs/
│   └── DESIGN.md                      # 本文档
├── copp_py/                           # 包名
│   ├── __init__.py
│   │
│   ├── path/                          # 对应 Rust path/
│   │   ├── autodiff.py                # Jet3 等价：三阶前向自动微分标量
│   │   ├── path_core.py               # Path：parametric / spline 统一接口
│   │   └── spline.py                  # Hermite 样条封装（基于 scipy）
│   │
│   ├── robot/                         # 对应 Rust robot/
│   │   ├── robot_core.py              # Robot：约束摄入 API（with_axial_*）
│   │   ├── kinematics.py              # KinematicsModel 协议 + 默认 DH/POE 实现
│   │   └── dynamics.py                # DynamicsModel 协议（逆动力学，供力矩约束/热能目标）
│   │
│   ├── commands/                      # ★ 指令层（Rust 版没有，本框架新增）
│   │   ├── base.py                    # MotionCommand 抽象基类 + PathSegment
│   │   ├── joint_move.py              # JointMoveCommand（关节运动指令）
│   │   ├── linear_move.py             # LinearMoveCommand（直线指令）
│   │   └── circular_move.py           # CircularMoveCommand（圆弧指令）
│   │
│   ├── copp/                          # 核心求解命名空间（复刻 src/copp/）
│   │   ├── constraints.py             # Constraints：站点索引约束存储
│   │   ├── objectives.py              # CoppObjective 枚举等价（Time/ThermalEnergy/...)
│   │   ├── copp2/
│   │   │   ├── formulation.py         # Topp2Problem / Copp2Problem 构造
│   │   │   ├── interpolation.py       # a_to_b / s_to_t / t_to_s（2阶）
│   │   │   ├── dp2/
│   │   │   │   ├── reach_set2.py      # 后向/双向可达集 DP
│   │   │   │   └── topp2_ra.py        # RA 前向贪心
│   │   │   └── opt2/
│   │   │       └── copp2_socp.py      # COPP2 cvxpy/SOCP 接口
│   │   └── copp3/
│   │       ├── formulation.py         # Topp3Problem / Copp3Problem（含线性化）
│   │       ├── interpolation.py       # s_to_t / t_to_s（3阶）
│   │       └── opt3/
│   │           ├── topp3_lp.py
│   │           ├── topp3_socp.py
│   │           └── copp3_socp.py
│   │
│   ├── planning/                      # ★★ 核心：轨迹规划模块
│   │   ├── __init__.py
│   │   ├── trajectory_planner.py      # TrajectoryPlanner 门面类
│   │   ├── path_builder.py            # 指令序列 → 统一几何路径 q(s)
│   │   ├── segment.py                 # PathSegment / SegmentBoundary 数据结构
│   │   ├── options.py                 # PlannerOptions（方法选择、边界条件、Verbosity）
│   │   └── result.py                  # TrajectoryResult 数据结构
│   │
│   ├── solver_backend/                # 求解器后端抽象
│   │   ├── base.py                    # SolverBackend 协议
│   │   ├── cvxpy_backend.py           # cvxpy + clarabel/ecos/scs
│   │   └── rust_ffi_backend.py        # 可选：PyO3 调用本仓库 Rust 实现
│   │
│   └── diag/
│       ├── errors.py                  # CoppError 层次（对应 diag/error.rs）
│       └── diagnostics.py             # Verbosity + 日志
│
└── tests/
    ├── conftest.py
    ├── unit/
    │   ├── test_path_autodiff.py
    │   ├── test_path_spline.py
    │   ├── test_robot_constraints.py
    │   ├── test_commands_joint_move.py
    │   ├── test_commands_linear_move.py
    │   ├── test_commands_circular_move.py
    │   ├── test_topp2_ra.py
    │   ├── test_reach_set2.py
    │   ├── test_copp2_socp.py
    │   ├── test_topp3_lp.py
    │   ├── test_topp3_socp.py
    │   └── test_copp3_socp.py
    ├── integration/
    │   ├── test_trajectory_planner_joint_only.py
    │   ├── test_trajectory_planner_mixed_commands.py
    │   └── test_trajectory_planner_objectives.py
    └── benchmark/
        └── test_random_spline_benchmark.py   # @pytest.mark.slow，对应 Rust ignored 基准测试
```

**模块依赖方向**（与 Rust 版一致，单向、无环）：

```
commands/  ──┐
             ├─→ planning/path_builder.py ─→ path/ ─→ robot/ ─→ copp/{copp2,copp3} ─→ solver_backend/
robot/kinematics.py, robot/dynamics.py ──┘
                                                        planning/trajectory_planner.py 编排以上全部
```

## 4. 指令体系（`commands/`）

这是 Python 框架相对 Rust 版本**新增的一层**，目的是把"做什么运动"的用户意图翻译成"路径长什么样"的几何描述，翻译结果统一交给 `planning/path_builder.py`。

### 4.1 `MotionCommand` 抽象基类（`commands/base.py`）

```python
class MotionCommand(Protocol):
    def to_segment(self, kin: KinematicsModel, u_samples: np.ndarray) -> PathSegment:
        """将指令在局部参数 u∈[0,1] 上采样，返回关节空间路径段。"""

@dataclass
class PathSegment:
    q: np.ndarray            # (dim, n) 关节位置采样
    dq_du: np.ndarray        # 对局部参数 u 的一阶导（供拼接时估计边界导数）
    ddq_du: np.ndarray       # 二阶导
    dddq_du: np.ndarray | None
    length_hint: float       # 段的弧长/角度估计，用于跨段重参数化时的相对权重
    boundary_velocity: tuple[np.ndarray, np.ndarray] | None = None   # 段起止建议速度方向（用于段间 C¹ 连接）
    local_constraint_overrides: ConstraintOverrides | None = None    # 该段专属限速/限力矩（可选）
```

`PathSegment` 的导数是对**局部参数 u**（每条指令各自的 [0,1]）求的，而不是全局弧长参数 `s`；跨段拼接统一到 `s` 时由 `path_builder.py` 按链式法则换算（对应 Rust `Path` 中弧长重参数化的思路）。

### 4.2 `JointMoveCommand`（关节运动指令）

```python
@dataclass
class JointMoveCommand(MotionCommand):
    q_start: np.ndarray            # (dim,) 起始关节角
    q_end: np.ndarray              # (dim,) 终止关节角
    boundary: Literal["stationary", "continuous"] = "stationary"
```

- 最简单情形：不需要逆运动学，直接在关节空间做插值。
- 单段内部用**五次多项式**（两端速度、加速度可指定，默认 stationary 即两端 `q̇=q̈=0`）参数化 `q(u)`，与 Rust `path/spline.rs` 的 Hermite 样条（阶次 `p=5`）复用同一实现。
- `boundary="continuous"` 用于指令序列中间段，此时边界导数由前一/后一指令的采样切向量决定，保证多指令拼接后的路径在关节空间 C¹（若约束允许，C²）连续。

### 4.3 `LinearMoveCommand`（直线指令）

```python
@dataclass
class LinearMoveCommand(MotionCommand):
    pose_start: Pose             # 位置 (3,) + 姿态四元数/旋转矩阵
    pose_end: Pose
    frame: Literal["base", "tool"] = "base"
    ik_seed: np.ndarray | None = None    # IK 数值解的初始关节角猜测
```

翻译流程：

1. 位置：`p(u) = (1-u)·p_start + u·p_end`（Cartesian 直线）。
2. 姿态：`scipy.spatial.transform.Slerp` 在 `p_start`、`p_end` 姿态间做四元数球面插值 `R(u)`。
3. 对采样得到的 `u_samples` 逐点调用 `kin.inverse_kinematics(pose(u), seed=...)`，相邻点用上一点的解作为下一点 IK 的 `seed`（保证解分支连续，避免关节跳变）。
4. 对得到的 `q(u)` 序列做有限差分或局部多项式拟合估计 `dq/du, d²q/du², d³q/du³`（IK 本身不提供解析导数时的兜底方案）；若 `KinematicsModel` 提供解析雅可比 `J(q)` 与其导数，则优先用 `q̇ = J⁻¹ ẋ` 解析计算，精度更高、也是后续优化的自然切入点。
5. 奇异位形检测：若 `J(q)` 条件数超过阈值，抛出 `CoppError.KinematicSingularity`，并在异常信息中给出发生的 `u` 值，便于用户调整路径或改用关节运动指令绕过奇异点。

### 4.4 `CircularMoveCommand`（圆弧指令）

```python
@dataclass
class CircularMoveCommand(MotionCommand):
    pose_start: Pose
    pose_end: Pose
    aux: Pose | tuple[np.ndarray, np.ndarray]
    # aux 二选一：
    #   - 空间中的第三点 pose_via（三点定圆弧，工业机器人常见写法）
    #   - (center, normal)：显式给出圆心与圆弧所在平面法向量
    direction: Literal["shortest", "ccw", "cw"] = "shortest"
```

翻译流程：

1. 若给定三点（起点、终点、途经点），先求解圆心 `c`、半径 `r`、法向量 `n`（三点确定唯一圆，退化共线情形需抛出 `CoppError.DegenerateArc`）。
2. 在圆弧平面内建立正交基 `(e1, e2)`（`e1` 指向起点方向，`e2 = n × e1`），弧上任一点：
   `p(θ) = c + r·(cos θ · e1 + sin θ · e2)`，`θ` 从起点角 `θ0` 单调过渡到终点角 `θ1`（按 `direction` 决定走劣弧/优弧/指定旋向）。
3. 局部参数 `u ↦ θ(u) = θ0 + u·(θ1-θ0)` 线性映射，位置对 `u` 的导数解析可得（三角函数求导，无需数值差分），这一点上圆弧指令天然适合复用 `path/autodiff.py` 的 `Jet3` 机制：把 `θ(u)` 送入 `Jet3` 通道，`p(θ)` 的三阶导直接解析获得。
4. 姿态插值同直线指令（SLERP），随后同样逐点 IK。
5. 圆弧半径、圆心退化情况（起点=终点、共线三点）作为单元测试的显式边界用例（见第 9 节）。

### 4.5 指令序列 → 统一路径

`planning/path_builder.py::build(commands: list[MotionCommand]) -> Path` 的职责：

1. 为每条指令分配全局弧长区间 `[s_i, s_{i+1}]`（区间长度按 `PathSegment.length_hint` 归一化，保证 `s∈[0,1]` 整体覆盖）。
2. 在段边界处按 `boundary` 策略拼接局部导数到全局 `s` 参数下（链式法则：`dq/ds = (dq/du)·(du/ds)`，`du/ds` 为该段的局部→全局缩放常数，二阶、三阶导同理逐级展开，公式形式与 `code_reading_guide.md` 第 2 节的关节-路径关系一致）。
3. 输出与 Rust `Path::evaluate_up_to_2nd/3rd` 同结构的 `PathDerivatives`（`q, dq, ddq, dddq`，形状 `(dim, N)`），交给 `robot/robot_core.py` 摄入约束。
4. 记录每段的 `(idx_start, idx_end)` 与其 `local_constraint_overrides`，供 `Robot` 在对应站点区间叠加/覆盖限速限力矩（对应 Rust 环形缓冲支持"滑动窗口分段约束"的设计思路；本版本先做一次性静态叠加，见第 11 节的后续扩展）。

## 5. 轨迹规划模块（`planning/`）——核心设计

### 5.1 `TrajectoryPlanner` 门面类

```python
class TrajectoryPlanner:
    def __init__(self, dim: int, kin: KinematicsModel | None = None, dyn: DynamicsModel | None = None): ...

    def add_command(self, cmd: MotionCommand) -> "TrajectoryPlanner": ...   # 链式调用，累积指令队列

    def set_velocity_limits(self, v_max, v_min) -> "TrajectoryPlanner": ...
    def set_acceleration_limits(self, a_max, a_min) -> "TrajectoryPlanner": ...
    def set_jerk_limits(self, j_max, j_min) -> "TrajectoryPlanner": ...      # 仅 3 阶方法需要
    def set_torque_limits(self, tau_max, tau_min) -> "TrajectoryPlanner": ...  # 需要 dyn

    def plan(
        self,
        method: Literal["topp2_ra", "copp2_socp", "topp3_lp", "topp3_socp", "copp3_socp"],
        objectives: list[CoppObjective] | None = None,   # 仅 COPP 系列需要
        options: PlannerOptions | None = None,
    ) -> TrajectoryResult: ...
```

设计要点：

- **不可变构建期 vs 求解期分离**：`add_command`/`set_*_limits` 只累积描述，真正的路径采样、约束摄入、求解在 `plan()` 内一次性完成，便于同一套指令+约束用不同 `method` 反复求解做算法对比（这正是 `docs/code_reading_guide.md` 里 benchmark 对比表格的生成方式）。
- **方法与阶数的一致性检查**：`plan()` 入口校验 `method` 与已设置约束的阶数是否匹配（例如设置了 `jerk` 限制却调用 `topp2_ra` 应给出明确警告而非静默忽略），错误类型归入 `diag/errors.py`。
- **3 阶方法的隐式两次 SCP 迭代**：`plan(method="topp3_lp"/"topp3_socp"/"copp3_socp", ...)` 内部自动：
  1. 先以 `topp2_ra` 求 `a_ref`（若用户已经跑过 2 阶方法，可通过 `options.linearization_reference` 复用，避免重复计算）；
  2. `Robot.constraints.amax_substitute(a_ref)` 收紧速度上界；
  3. 第一次线性化 + 求解得到 `(a1, b1)`；
  4. 用 `a1` 重新线性化 + 求解得到 `(a2, b2)`（第二次 SCP，通常作为最终解）；
  5. `options.scp_iterations` 允许配置迭代次数（默认 2，与 Rust 版一致），便于测试中验证"迭代次数越多目标值单调不增"这类性质。
- **边界条件默认值**：起止 `a=0, b=0`（静止边界），对应工业场景"从静止开始、到静止结束"的常见需求；`PlannerOptions` 暴露自定义边界的入口（非静止边界衔接多段轨迹，为后续扩展点）。

### 5.2 `TrajectoryResult` 数据结构（`planning/result.py`）

```python
@dataclass
class TrajectoryResult:
    s_grid: np.ndarray          # (N,) 路径参数网格
    a_profile: np.ndarray       # (N,) 或 (N,) 视方法而定；3 阶方法额外有 b_profile
    b_profile: np.ndarray | None
    num_stationary: tuple[int, int] | None   # 仅 3 阶方法：起止静止缓冲站数
    t_final: float
    t_s: np.ndarray              # (N,) 到达每个 s_k 的时间
    method: str
    objective_value: float | None
    solver_status: str
    solve_time_seconds: float

    def sample_uniform_time(self, dt: float) -> "TimeSamples": ...
    def joint_position(self, t: np.ndarray) -> np.ndarray: ...
    def joint_velocity(self, t: np.ndarray) -> np.ndarray: ...
    def joint_acceleration(self, t: np.ndarray) -> np.ndarray: ...
```

`joint_velocity`/`joint_acceleration` 的重建公式与 `code_reading_guide.md` 第 11 节"关节空间重建"完全一致：

```
q̇(t)  = q'(s(t)) · √a(s(t))
q̈(t)  = q''(s(t)) · a(s(t)) + q'(s(t)) · b(t)
```

`b(t)`：2 阶方法由 `a` 的差分估计（`da/dt / (2ṡ)`），3 阶方法直接来自 `b_profile` 插值——这一区分在 `TrajectoryResult` 内部对用户透明。

### 5.3 与 `robot/` 的边界

`TrajectoryPlanner` 不直接操作 `Constraints` 环形缓冲区，而是持有一个 `Robot` 实例，约束设置 API（`set_velocity_limits` 等）本质是对 `Robot.with_axial_*` 的门面封装。这保持了与 Rust 版相同的**职责边界**：`planning/` 只负责编排，约束语义和存储仍归 `robot/` 所有。

## 6. COPP 求解器复刻（`copp/`）

逐一对应 `docs/code_reading_guide.md` 中的算法描述，Python 侧的实现策略：

| Rust 组件 | Python 对应 | 实现策略 |
|-----------|-------------|----------|
| `dp2/reach_set2.rs`（后向可达集） | `copp2/dp2/reach_set2.py` | 逐站点用自实现 2D LP（Seidel 增量，见 `math/lp.py`，对应 2.1 节自实现 AD 一样的教学/性能考量）求 `a_k` 上界；核心循环用 `numpy` 向量化约束系数收集，LP 部分保持逐点（因为存在数据依赖，无法整体向量化） |
| `dp2/topp2_ra.rs`（前向贪心） | `copp2/dp2/topp2_ra.py` | 结构 1:1 复刻：收集约束行 → 代入 `a_{k-1}` 化为 1D LP → 与后向区间取交 → 贪心取最大值 |
| `opt2/copp2_socp.rs` | `copp2/opt2/copp2_socp.py` | 用 `cvxpy` 变量 `a = cp.Variable(n+1, nonneg=True)`，时间目标通过引入 SOC 辅助变量表达 `1/√a` 的上界（`cvxpy` 的 `cp.SOC` 原语或直接用其 `inv_pos`/`quad_over_lin` 原子重述，等价于 Rust 中手工展开的锥形式） |
| `copp3/formulation.rs`（线性化） | `copp3/formulation.py` | 完全复刻线性化公式（第 8.1 节），包括 `a_linearization_floor` 防止除零 |
| `opt3/topp3_lp.rs` / `topp3_socp.rs` / `copp3_socp.rs` | `copp3/opt3/*.py` | 决策变量 `x=[a, b]`，动力学等式约束 `a_{k+1}=a_k+2Δs_k b_k` 用 `cvxpy` 的等式表达；LP 版本目标为线性，SOCP 版本目标含 `1/√a` 型二次锥项 |
| `objectives.rs` | `copp/objectives.py` | `CoppObjective` 用 Python `Enum`/`dataclass` 家族（`Time`, `ThermalEnergy`, `TotalVariationTorque`, `Linear`）表示，各自提供 `to_cvxpy_expr(a, b, aux_vars)` 方法，由具体 solver 组装 |
| `constraints.rs`（环形缓冲） | `copp/constraints.py` | 第一版**不实现环形缓冲**，直接用定长 `numpy` 数组存储（一次性规划场景不需要滑动窗口）；接口预留 `pop_front`/`expand_capacity` 方法签名，行为可后续补齐（见第 11 节），保证未来切换到在线滑动窗口规划时上层 API 不变 |
| `clarabel_backend.rs` | `solver_backend/cvxpy_backend.py` | 封装 `problem.solve(solver=cp.CLARABEL, ...)`，`ClarabelOptionsBuilder` 等价为 `dataclass SolverOptions(allow_almost_solved: bool = True, ...)`，`is_allow(status)` 判定逻辑照搬 |
| `diag/`（Verbosity/错误） | `diag/diagnostics.py`、`diag/errors.py` | `Verbosity` 用 `IntEnum(Silent=0, Summary=1, Debug=2, Trace=3)`；错误层次用异常类继承树 `CoppError -> {ConstraintError, PathError, InvalidInputError, SolverStatusError}` |

**数值一致性目标**：单元测试要求 Python 复刻算法在相同输入（同一组随机样条路径 + 同一组约束）下，与 Rust 版本输出的 `t_final`、`objective_value` 相对误差 `< 1e-4`（与 README 中 TOPP2-RA 对全局最优的误差量级一致）；这是判断"复刻是否正确"的量化标准，具体做法见第 9.3 节。

## 7. 求解器后端抽象（`solver_backend/`）

```python
class SolverBackend(Protocol):
    def solve_socp(self, problem: CvxpyProblemSpec, options: SolverOptions) -> SolverSolution: ...
    def solve_lp(self, problem: CvxpyProblemSpec, options: SolverOptions) -> SolverSolution: ...
```

- `cvxpy_backend.py`：默认后端，`solver=cp.CLARABEL`，失败时按 `options.fallback_solvers=["ECOS","SCS"]` 顺序重试。
- `rust_ffi_backend.py`：**设计占位**，通过 `maturin` 构建的 `copp` PyO3 扩展直接调用 Rust `topp2_ra`/`copp2_socp` 等函数，输入输出格式与 `cvxpy_backend` 保持一致的 `SolverSolution` 结构，使 `planning/` 层可以通过一个 `options.backend: Literal["cvxpy","rust"]` 开关无感切换。这是长期性能路径（见第 10 节），本次不要求实现。

抽象后端的意义：`planning/trajectory_planner.py` 与 `copp/copp2, copp/copp3` 中的求解函数只依赖 `SolverBackend` 协议，不感知具体是 `cvxpy` 还是 Rust FFI，方便测试中用假后端（`FakeSolverBackend`）做纯逻辑单元测试，不必每次都跑真实凸优化。

## 8. 端到端示例（伪代码）

```python
from copp_py.commands import JointMoveCommand, LinearMoveCommand, CircularMoveCommand
from copp_py.planning import TrajectoryPlanner
from copp_py.copp.objectives import Time, ThermalEnergy
from copp_py.robot.kinematics import DHKinematics

kin = DHKinematics.from_dh_table(dh_params)   # 6 自由度机器人

planner = TrajectoryPlanner(dim=6, kin=kin)
planner.add_command(JointMoveCommand(q_start=q_home, q_end=q_approach))
planner.add_command(LinearMoveCommand(pose_start=pose_approach, pose_end=pose_via))
planner.add_command(CircularMoveCommand(pose_start=pose_via, aux=pose_mid, pose_end=pose_depart))

planner.set_velocity_limits(v_max=[1.0]*6, v_min=[-1.0]*6)
planner.set_acceleration_limits(a_max=[2.0]*6, a_min=[-2.0]*6)
planner.set_jerk_limits(j_max=[10.0]*6, j_min=[-10.0]*6)

result = planner.plan(
    method="copp3_socp",
    objectives=[Time(weight=1.0), ThermalEnergy(weight=0.1, axis_weights=motor_weights)],
)

print(result.t_final, result.solver_status)
t = np.arange(0.0, result.t_final, 1e-3)
q_t, qd_t, qdd_t = result.joint_position(t), result.joint_velocity(t), result.joint_acceleration(t)
```

对应 README "Quick Start" 中 Rust 示例的 Python 版意图一致，但输入从"直接给定解析路径闭包"变成了"指令序列"，这是本框架相对 Rust 库的核心增量。

## 9. 测试用例设计（`tests/`）

### 9.1 单元测试（`tests/unit/`）

| 文件 | 覆盖点 |
|------|--------|
| `test_path_autodiff.py` | `Jet3` 等价类对 `sin/cos/exp/ln/sqrt/pow/+-*/` 的三阶导数与 `scipy.misc.derivative`/解析导数数值对照（容差 `1e-8`）；链式法则组合表达式（如 `sin(2πs)·exp(s)`）正确性 |
| `test_path_spline.py` | Hermite 样条 3/5/7 阶在给定边界导数下，样条值与端点约束吻合；随机路点下样条连续性（C²/C⁴/C⁶）数值检验（差分逼近） |
| `test_robot_constraints.py` | 逐条验证第 6 章约束映射公式（速度/加速度/急动度/力矩）与 `code_reading_guide.md` 第 6 节公式的数值一致性；边界情况 `q'_i=0`（对应 `a_max=+∞`） |
| `test_commands_joint_move.py` | 起止点吻合；`boundary="stationary"` 时两端速度为零；插值路径在关节界内单调（无超调） |
| `test_commands_linear_move.py` | IK 往返一致性（`kin.forward(kin.inverse(pose)) ≈ pose`）；路径为直线（正投影残差为零）；奇异位形抛出 `KinematicSingularity` |
| `test_commands_circular_move.py` | 弧上采样点到圆心距离恒为半径（容差 `1e-6`）；角度覆盖与 `direction` 参数（`ccw`/`cw`/`shortest`）一致；共线退化输入抛出 `DegenerateArc` |
| `test_topp2_ra.py` | 与 Rust `examples/topp2_ra.rs` 相同的 3 轴 Lissajous 路径 + 对称 v/a=1 约束，验证 `t_final` 数值（对照 README 已发布数据）；`a_profile` 逐点满足约束不等式 |
| `test_reach_set2.py` | 后向可达集 `a_max[k]` 的单调性质（约束越紧 `a_max` 越小）；退化路径（恒速直线）下解析解对照 |
| `test_copp2_socp.py` | `Time` 目标下与 `topp2_ra` 结果近似一致（凸问题应收敛到相近最优）；`ThermalEnergy` 目标下验证目标值确实低于纯时间目标对应的能耗 |
| `test_topp3_lp.py` / `test_topp3_socp.py` | 两次 SCP 迭代后目标值单调不增；线性化系数公式与 `code_reading_guide.md` 8.1 节公式数值核对；`num_stationary` 在静止边界下非零，在非静止边界下为零 |
| `test_copp3_socp.py` | 力矩约束下逆动力学系数分解（`coeff_a/coeff_b/coeff_g`）与直接 RNE 计算结果核对；多目标线性组合权重生效性 |

### 9.2 集成测试（`tests/integration/`）

| 文件 | 覆盖点 |
|------|--------|
| `test_trajectory_planner_joint_only.py` | 纯 `JointMoveCommand` 序列（不需要 IK），端到端 `plan()` 跑通全部 5 种 `method`，检查返回的 `TrajectoryResult` 字段完整性与约束满足性（逐点检查 `q̇,q̈,q⃛` 在界内） |
| `test_trajectory_planner_mixed_commands.py` | `JointMove + LinearMove + CircularMove` 混合序列，验证段间位置连续（`max|q(s_i^-) - q(s_i^+)| < tol`）与速度连续（若 `boundary="continuous"`） |
| `test_trajectory_planner_objectives.py` | 同一指令序列分别用 `Time` 与 `ThermalEnergy` 目标求解，验证目标切换确实改变了 `a(s)` 分布（而非退化为同一解）；验证 `method` 与约束阶数不匹配时的显式报错路径 |

### 9.3 基准/交叉验证测试（`tests/benchmark/`）

`test_random_spline_benchmark.py`：

- 对应 Rust `tests/test_random_spline.rs`（`--ignored`，慢速基准），标记 `@pytest.mark.slow`，默认 CI 不跑，本地/夜间任务手动触发。
- 生成与 README benchmark 相同规模的数据集（100 条随机 7-DOF 样条路径 × 1000 段离散点），对每种 `method` 记录计算耗时与目标值，输出统计摘要（mean ± std），与仓库 README 中 Rust 基准表格并列比较，用于持续追踪 Python 复刻实现与 Rust 原版之间的性能/质量差距。
- 若安装了 `rust_ffi_backend` 依赖（可选 extras），额外对同一数据集跑 Rust 后端，做 Python vs Rust 的数值一致性断言（`t_final` 相对误差 `<1e-4`），否则跳过该断言并给出 `pytest.skip` 提示。

### 9.4 测试基础设施（`conftest.py`）

- `fixture: lissajous_path_3axis()` —— 复刻 README/`examples/topp2_ra.rs` 中的确定性 3 轴路径，作为跨测试文件复用的"黄金输入"。
- `fixture: symmetric_limits(dim)` —— 生成对称 v/a/jerk 限制。
- `fixture: fake_kinematics_2r()` —— 平面 2R 机械臂的解析正逆运动学（不依赖 `roboticstoolbox-python`），用于不想引入外部运动学库的单元测试。
- `fixture: solver_backend` —— 参数化 `["cvxpy"]`（以及在安装了对应 extras 时的 `"rust"`），使同一组测试可以在两种后端下运行。

## 10. 与 Rust 版本的关系与演进路径

1. **第一阶段（本设计覆盖范围）**：纯 Python 实现，`cvxpy + clarabel-python` 作为求解后端，用于快速原型、教学、与非 Rust 技术栈的机器人系统集成。
2. **第二阶段（可选）**：`solver_backend/rust_ffi_backend.py` 落地，通过 `maturin` 把本仓库 Rust crate 编译为 Python 扩展模块，`planning/` 层通过协议无缝切换，兼顾开发便利性与生产性能。
3. **第三阶段（可选）**：若指令层（`commands/`）证明足够通用，可考虑反向贡献回 Rust 版本（新增 `copp::commands` crate 子模块），使 Rust/Python 两侧共享同一套"指令 → 路径"翻译逻辑，避免长期维护两份实现。

## 11. 待决问题 / 后续扩展点

以下问题在本次设计中给出建议方向，但不阻塞轨迹规划模块的落地，留待后续迭代：

- **在线滑动窗口规划**：Rust 版 `Constraints` 环形缓冲支持 CNC 场景下逐段推送约束、边规划边执行；Python 第一版为一次性静态规划，`copp/constraints.py` 已预留接口签名但不实现窗口滑动语义。
- **段间非静止边界拼接**：多指令序列若要求段边界速度非零连续（例如高速产线连续走多段直线不减速），需要 `PlannerOptions` 暴露自定义边界 `(a_boundary, b_boundary)` 并在 `path_builder.py` 中传递给每个内部子问题，当前设计只描述了默认的静止边界情形。
- **IK 数值稳定性与解分支选择**：`LinearMoveCommand`/`CircularMoveCommand` 依赖数值 IK 的连续性假设，对于强非线性/多解机器人（如 6R 非球形腕），需要更严格的解分支跟踪策略，建议后续引入解析 IK（`roboticstoolbox-python` 对常见机器人族提供）优先于数值 IK。
- **性能剖析**：`copp2/dp2/*.py` 的逐点 LP 循环是否需要 `numba`/`Cython` 加速，取决于第一版基准测试（9.3 节）结果，暂不预判。

<!-- END DESIGN.md -->

---

## 4. 如何在新项目中继续

1. 把本文件（`SESSION_HANDOFF.md`）粘贴或上传到新会话的第一条消息中，作为上下文。
2. 告诉 Claude 你想从哪一步继续，例如：
   - "基于这份设计，先实现 `copp_py/path/autodiff.py`（Jet3 等价类）和对应单元测试"
   - "先落地 `commands/` + `planning/path_builder.py`，只支持关节运动指令，暂不做 IK"
   - "把 `copp2/dp2/topp2_ra.py` 按第 6 节的映射表实现，并对照 Rust `examples/topp2_ra.rs` 的数值做交叉验证"
3. 如果新项目就是当前 `copp` 仓库的另一个 clone/worktree，直接把 `python/DESIGN.md` 和 `python/README.md` 两个文件复制过去即可，无需使用本文档的第 3 节全文（那只是为了脱离原仓库时的自包含性）。
4. 注意：`python/DESIGN.md`、`python/README.md` 在原仓库中**尚未 git commit**，如果需要保留，记得在原仓库里先提交，或者手动确认已经复制到新位置。

---

## 5. 一些可能有用的其它上下文（本会话未展开，但可能相关）

- 原仓库根目录下还有 `docs/code_reading_guide.md`（更详细的 Rust 内部实现讲解，含约束矩阵组装、Clarabel 锥类型等），如果继续深挖 Rust 侧算法细节，应参考该文件而不是仅凭本文档摘录。
- 原仓库 `README.md` 中有 TOPP2-RA/COPP2-SOCP/TOPP3-LP/TOPP3-SOCP/COPP3-SOCP 的基准数据表格（计算耗时、轨迹时间/目标值的 mean±std），可作为 Python 复刻版本做性能/质量对照的基准真值。
- 原仓库当前分支：`feature/wc/add_save_data_plot_script`，工作区中还有一些未跟踪/暂存的文件（`python/DESIGN.md`、`python/README.md` 为本次新增，另有 `docs/paper_notes.md` 等未跟踪文件，与本次任务无关，是会话开始前就存在的工作区状态）。
