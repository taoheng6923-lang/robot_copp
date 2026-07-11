# robot_copp — M2：指令层 + 降维层（可跑）

`trajectory-planning/path/` 的 **commands（三类运动指令）+ lowering（IK 降维）**，
打通"指令序列 → 笛卡尔路径 → UR5 关节路径导数 → copp SPLP 求解"的完整链路。
对应 [`python_framework.md`](./python_framework.md) §9 里程碑 M2 与
[`robot_copp_design.md`](./robot_copp_design.md) §3/§5。M1 内核见
[`README_M1.md`](./README_M1.md)。

## 已实现（M2 范围）

| 模块 | 内容 | 设计对应 |
|------|------|----------|
| `path/types.py` | `CartesianSamples`（位姿+1~3 阶导采样）/ `CartesianPath`、`JointSpacePath` 协议 / `PathDerivatives`（q,q',q'',q''' + cv/cw） | framework §4 / 设计 §5.5 |
| `path/errors.py` | `PathError` 层次（零长/退化弧/不可达/IK 跳变/衔接错位） | framework §5.10 |
| `path/commands/base.py` | `Section` + `MotionCommand` 协议 + SLERP 工具（世界系 ω̂ 恒定） | 设计 §3.4 |
| `path/commands/joint_move.py` | `JointMoveCommand`：**线性关节几何**（与设计 §3.1 的差异见下） | 设计 §3.1 |
| `path/commands/linear_move.py` | `LinearMoveCommand`：直线 + SLERP，全解析导数；纯姿态调整（Δp=0）经 `rot_scale` 退化安全 | 设计 §3.2 |
| `path/commands/circular_move.py` | `CircularMoveCommand`：三点定圆（方向由 via 唯一确定）/ (center,normal)+direction；弧长参数化解析导数；共线/重合/离面显式报错 | 设计 §3.3 |
| `path/commands/assemble.py` | `build_sections` + `lower_sections`：G0 衔接校验 + IK seed 链逐段降维 | 设计 §3.5 |
| `path/lowering/sampling.py` | 曲率驱动自适应采样（弦高 e≈‖r''‖Δs²/8 模型，网格必含 s_breaks） | 设计 §5.1 |
| `path/lowering/ik.py` | 连续解 IK（seed 链）+ FK 回代校验 + 跳变校验 | 设计 §5.2 |
| `path/lowering/derivatives.py` | Jacobian 链式法则求 q',q'',q'''（J'/J'' 方向差分）+ `lower_cartesian`/`lower_joint` 驱动 | 设计 §5.3 |
| `path/lowering/singularity.py` | σ_min/σ_max 奇异检测 + DLS 逆（λ=0.05）+ singular 标记 | 设计 §5.4 |

机器人本体（真实 UR5 DH 运动学 + 解析 IK）在 M1 阶段已提前落地（`robot/ur5.py`），
本里程碑的 lowering 直接以它为 `KinematicsModel` 后端。

## 与设计文档的两处已确认修正

1. **三阶链式法则的系数 2**（设计 §5.3，已改）：对 `J'q' + Jq'' = r₂` 再求导，
   `J q''' = r₃ − 2·J'q'' − J''q'`——设计文档 v0.5 前漏写系数 2。实现按正确式，
   并在自测中以 q 序列有限差分交叉验证（若按原式实现，该断言会以 O(10%) 量级失败）。
2. **JointMove 用线性关节几何而非五次多项式**（设计 §3.1 的偏离，有意）：五次
   多项式参数化 q(u) 是"时间域直接参数化"的经典做法；TOPP 框架下几何与时间律
   解耦——时间平滑（rest-to-rest、jerk 受限）由 copp 时间律保证，几何上的五次
   多项式反而使端点 q'(s)=0（参数化退化、非正则）。段间几何平滑属 M3 blending。

## M2 语义边界（无 blending）

相邻指令段只做 **G0 精确衔接校验**（位置/姿态/关节角必须严丝合缝，错位抛
`JunctionMismatchError`），每段**独立 rest-to-rest 求解**——切向不连续的角点以
停顿衔接，物理严格可行但非全局时间最优。不停顿的 G2 平滑过渡（zone 截短 +
最优五次 Hermite）属 M3（设计 §4）。

IK 解分支连续性：首个笛卡尔段用调用方 `q_seed` 作 IK 种子，其后各段用上一段
末端关节角，跨段解分支一致（设计 §12.2）。

三点圆弧的扫角符号约定：`sweep` 相对 **via 定出的法向 n̂ = (via−p0)×(p1−p0)**
（几何不变量是弧长与"经过 via"，不是世界系下的正负号）。

## 验证方法（每步都有独立交叉验证）

`path/self-test/test_lowering.py`（单一用例 `test_lowering`）：

- **夹具自检**：解析测试路径（位置三角曲线 + 变轴姿态 Rx(α(s))·Ry(β(s))，
  ω̂/ω̂'/ω̂'' 全解析）先用 R(s) 的有限差分验证自身导数推导；
- **IK 回代**：FK(q_k) 复现路径位姿（~1e-9）；细网格下 q 准连续；
- **链式法则**：q 序列有限差分复现 dq/ddq/dddq（尤其锁定三阶式系数 2）；
- **独立恒等式**：网格差分的 J'（不经被测方向差分代码）满足 J'q'+Jq''≈r₂；
- **自适应采样**：网格含段边界、逐区间弦高误差在容差内；
- **纯姿态调整**：p'≡0 时 cv≡0、TCP 线速度为零、全流程不崩；
- **端到端**：自适应网格 → `to_topp3_data`（TCP 速度模约束）→ `solve_splp`，
  迭代单调、速度/加速度/TCP 约束满足。

`path/self-test/test_commands.py`（单一用例 `test_commands`）：

- 三类指令的几何性质（直线残差 ~1e-15、弧上点距圆心恒 r、‖p'‖=1、‖p''‖=1/r、
  ω̂ 与 R 的 FD 一致、端点/via 吻合、劣弧优弧方向、ccw/cw/shortest 语义）；
- 全部退化路径显式报错（零长指令、共线三点、起终点重合、端点离面）；
- 段间衔接错位（笛卡尔 2mm 错位 / 关节 0.01rad 错位）抛 `JunctionMismatchError`；
- 混合序列端到端：JointMove + LinearMove + CircularMove 在 UR5 上逐段降维 +
  SPLP 求解，逐段约束满足、段间关节角衔接 <1e-8。

错误路径另经运行时探针确认：不可达站点 → `UnreachablePoseError`；
接近工作空间边界（肘部趋直、解快速摆动）→ `IkJumpError`。

## 运行

```bash
.\.venv\Scripts\python.exe -m pytest -q        # 全部三个自测（copp 内核 + lowering + commands）
# 或单独跑（额外落图到各自 self-test/output/）：
.\.venv\Scripts\python.exe trajectory-planning\path\self-test\test_lowering.py
.\.venv\Scripts\python.exe trajectory-planning\path\self-test\test_commands.py
```

- `path/self-test/output/lowering_test.png`：q(s)/q'(s)、自适应步长、cv/cw、
  速度剖面与 TCP 折算上界、SPLP 收敛。
- `path/self-test/output/commands_test.png`：三段 TCP 轨迹（3D，含 via）、
  各段 rest-to-rest 速度剖面、拼接后的关节速度时间历程。

## 最小用法

```python
import numpy as np
from robot import UR5Kinematics, Pose
from copp import load_robot_limits, solve_splp, SolveOptions
from path.commands import (JointMoveCommand, LinearMoveCommand,
                           CircularMoveCommand, build_sections, lower_sections)

kin = UR5Kinematics()
sections = build_sections([
    JointMoveCommand(q_home, q_a),
    LinearMoveCommand(pose_a, pose_b),            # pose_a 应等于 FK(q_a)
    CircularMoveCommand(pose_b, pose_c, via=v),
])
pds = lower_sections(sections, kin, q_seed=q_home)   # 每段一个 PathDerivatives

limits = load_robot_limits(v_tcp_max=0.5, w_tcp_max=2.5)
for pd in pds:                                       # M2：逐段 rest-to-rest
    data = limits.to_topp3_data(pd.s_grid, pd.dq, pd.ddq, pd.dddq,
                                tcp_geom=pd.tcp_geom())
    profile, hist = solve_splp(data, SolveOptions(n_iter=3))
```

## 仍待补齐（后续里程碑）

- **M3 blending**：zone 截短 + 最优五次 Hermite G2 过渡（不停顿过角点），
  `path/blending/` 目前仍是占位。
- **planner 门面**（M2+）：上面"最小用法"的手工串联收进
  `trajectory-planning/planner/planner.py`；轨迹合成/验证收进 `planner/synth`。
- **真实动力学**：力矩约束系数仍来自 `UR5RobotModel.torque_coeffs` 对角近似，
  RNE（`DynamicsModel`）未实现，故本里程碑端到端默认不启用力矩约束。
- **HLAW**（M5）：长序列三窗调度。
