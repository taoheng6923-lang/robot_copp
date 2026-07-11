# robot_copp — M1：SPLP 数值内核（可跑）

论文 **TOTP-SPLP**（分段线性目标 PLP + 序列线性化 + LP）的最小可运行内核。
对应 [`python_framework.md`](./python_framework.md) §9 里程碑 M1 与 [`robot_copp_design.md`](./robot_copp_design.md) §7。

> **顶层模块划分**（与 `python_framework.md` §2 描述的单包结构不同，实际实现按职责拆分）：
> - 顶层 `robot/` —— 机器人运动学/动力学本体（FK/IK/Jacobian/逆动力学），独立于其余三者
> - `trajectory-planning/`（纯目录容器，非 Python 包）下：
>   - `copp/` —— 纯 TOTP-SPLP 数值求解核心（types/options/constraints/solve/backend），不依赖 robot/path
>   - `path/` —— 路径构造（M2 已实现 commands/lowering，见 [`README_M2.md`](./README_M2.md)；blending 属 M3 未实现）
>   - `planner/` —— 调度/合成/门面（hlaw/synth/planner.py，M2+/M5，尚未实现）
>
> 三者（`path`/`copp`/`planner`）+ 顶层 `robot` 经 `trajectory-planning/planner/planner.py`
> 编排串联；`copp/` 只通过 `PathDerivatives`/`Topp3Data` 这类数据结构与 path/robot 交互，
> 不做反向导入。`trajectory-planning` 本身不是 Python 包（无 `__init__.py`，含连字符不可作为
> 模块名），只是把 `path`/`copp`/`planner` 三个顶层包归到一处的目录；本地开发经
> `pyproject.toml` 的 `[tool.pytest.ini_options] pythonpath` 把它加入 `sys.path`。

## 已实现（M1 范围）

`copp/{types,solve,backend/cvxpy_backend}` + `robot/{base,ur5}`：

| 模块 | 内容 | 论文/设计对应 |
|------|------|---------------|
| `copp/types.py` | `Topp3Data`（轴向约束输入）/ `Profile(a,b,c)` | framework §4 |
| `copp/solve/state.py` | 梯形动力学、速度上界 ā、时间权重 | Prop.1 / 设计 §7.1 |
| `copp/solve/seed.py` | 种子 a⁽⁰⁾（2 阶 LP，= topp2_ra 角色） | 论文 §5.1 |
| `copp/solve/linearize.py` | jerk 凹约束切线线性化 | 论文 **eq.32** |
| `copp/solve/plp_objective.py` | PLP 割线上包络 + 辅助变量 J_k + 下界 a≥δ0 | 论文 **eq.27 / 29d / Prop.3** |
| `copp/solve/lp_problem.py` | 单次 PLP-LP 组装（cvxpy） | 论文式 29 |
| `copp/solve/splp.py` | Algorithm 2 迭代（eq.30 停止） | 论文 **Algorithm 2** |
| `copp/solve/interp.py` | 解析插值 s↔t（含 Prop.2 静止段） | Prop.1 / **Prop.2** |
| `copp/constraints/model.py` | `RobotLimits`（机器人本体约束配置 → `Topp3Data`） | framework §5.5 |
| `copp/constraints/ingest.py` | TCP 速度模上界、关节力矩约束行（M4） | framework §5.5 / 设计 §6 |
| `copp/options.py` | `ConstraintFlags`（六类约束开关） | framework §5.10 |
| `copp/config.py` | 从 YAML 加载机器人逐关节约束 → `RobotLimits` | framework §5.5 |
| `robot/base.py` | `KinematicsModel` / `DynamicsModel` 协议 | framework §5.1 |
| `robot/ur5.py` | `UR5Kinematics`（真实 UR5 DH 运动学：fk/jacobian 解析解、ik 闭式解析逆解 8 支路 + 退化位姿 DLS 兜底，落地 `KinematicsModel` 协议）+ `UR5RobotModel`（真实 TCP 几何 + 对角近似动力学，供本 self-test 使用） | framework §5.1 |
| `copp/backend/cvxpy_backend.py` | CLARABEL 求解 | framework §5.9 |
| `copp/viz.py` | 结果可视化（2×3 概览图，可选 matplotlib） | framework §5.8 |

> **机器人本体**：M1 起用真实 UR5（标准 DH 参数、官方公开逐关节力矩上限），不再是
> 与关节数无关的"合成 N 轴"占位。`UR5Kinematics.fk`/`jacobian` 是解析解（已用有限差分
> 交叉验证到 ~1e-10）；`ik` 是**闭式解析逆解**（枚举全部 ≤8 支路取离 seed 最近者，
> O(1)，已用 200 组随机位形验证到 ~1e-13），退化位形（腕奇异/不可达）回退阻尼最小
> 二乘（DLS）数值迭代兜底。`UR5RobotModel.torque_coeffs` 仍是**对角近似**（下游集中
> 质量单摆臂估算，量级贴近真实 UR5 但非精确 RNE，无科氏/惯量耦合），真实
> `DynamicsModel`/RNE 待 M2+ 落地。`joint_path` 仍是合成随机轨迹（路径生成层 M2/M3
> 落地前的占位），但 `tcp_geometry`/`tcp_coeffs` 现在真实调用 `UR5Kinematics.jacobian`
> 沿该路径求值，不再是与关节角无关的虚构公式。见 `robot/ur5.py` 模块 docstring 的数据来源。

## 运行

```bash
cd robot_copp
pip install numpy scipy cvxpy clarabel pyyaml matplotlib   # 依赖（matplotlib 供可视化）
python trajectory-planning/copp/self-test/test_splp_kernel.py                  # 或 pytest（含可视化冒烟测试）
```

> 用例统一维护在 `trajectory-planning/copp/self-test/` 下（`test_splp_kernel.py`）。直接运行该脚本会跑全部
> 断言并把三张图落到 `trajectory-planning/copp/self-test/output/`。

**机器人约束配置**：本体逐关节约束（`vmax/amax/jmax/tau`）与两端边界默认从
[`configs/robot_ur5.yaml`](./configs/robot_ur5.yaml) 读取（UR5 六轴）。三档可信度
（详见该文件内注释）：`tau_max` 为 Universal Robots 官网公开的逐关节力矩上限（官方
权威）；`vmax` 取自 ROS-Industrial `ur5.urdf.xacro` 的关节速度限（驱动真实 UR5 硬件
的社区描述文件，逐关节不同，非统一近似值）；`amax/jmax` **无任何官方或社区数据**，
是按 vmax 的经验比例臆造的示例值，不代表真实机器人规格。
[`configs/robot_3axis.yaml`](./configs/robot_3axis.yaml) 保留作通用/教学示例。
TCP 速度模上界（`v_tcp_max/w_tcp_max`）作为“给定”参数在调用处传入（任务/工艺侧设定）：

```python
from copp import load_robot_limits
limits = load_robot_limits("configs/robot_ur5.yaml", v_tcp_max=0.6, w_tcp_max=0.9)
```

改机器人本体约束只改 YAML；测试的唯一约束来源见 `trajectory-planning/copp/self-test/limits_config.py`。

**约束开关（可选约束）**：六类约束（速度 / 加速度 / jerk / 力矩 / TCP 位置速度 /
TCP 姿态角速度）均可单独启用或关闭，开关放在
[`configs/comm_paras.yaml`](./configs/comm_paras.yaml) 的 `constraints` 节：

```python
from copp import load_constraint_flags, SolveOptions
flags = load_constraint_flags()                 # 读 comm_paras.yaml 的 constraints 节
solve_splp(data, SolveOptions(flags=flags))     # 关闭的约束求解时不施加
```

`ConstraintFlags` 默认全开；关闭某约束仅表示求解时不施加它（其上/下界数据仍可存在）。
被关闭约束的对应量便不再被限制（其超限比可 >1），时间相应变化。

**单一测试用例** `test_splp_kernel`（`trajectory-planning/copp/self-test/test_splp_kernel.py`）一次求解后完成全部断言
并输出三张图：SPLP 迭代 `t_final` **单调不增**并收敛、剖面形状/边界、rest-to-rest
（首末速度=0）、s↔t 有限且严格递增、M4（TCP 速度模 + 关节力矩）约束满足且至少一项绑定、
机器人配置加载。典型行为：jerk/acc 约束**恰好贴边**（超限比 ≈1.0）—— jerk-限时最优。

**可视化**（`copp.viz`，三张均由单一用例 `test_splp_kernel` 一次生成，便于对照分析）：

- `self-test/output/splp_test.png`（`plot_splp_result`，2×3）：① SPLP 收敛；② 速度剖面 ṡ(s)
  与上界 √ā；③ 路径加速度 b(s)；④ 约束利用率 vs s；⑤ 时间律 s(t)；⑥ 关节速度。
- `self-test/output/splp_limits_test.png`（`plot_kinematic_limits`，2×3，M4 数据）：① 关节速度
  ② 关节加速度 ③ 关节 jerk ④ **关节力矩** ⑤ TCP 位置速度模 ‖ṗ‖ ⑥ TCP 姿态角速度模
  ‖ω‖，各带约束线。信号用 `reconstruct_time_signals`（`interp.fine_profiles` 的**区间内
  细剖面** Prop.1/2）重构 —— q̈ 与 q⃛ **导数自洽**（网格点+线性插值会在静止段把
  `q̈∝σ^{1/3}` 画成 `∝σ`，使加速度起点斜率与 jerk 值对不上）。
- `self-test/output/fig4_interpolation.png`（`plot_fig4_interpolation`）：**复现论文 Fig.4**
  —— 静止起点（`a_s=b_s=0`）+ 非静止终点的区间解析插值。
  ① `a(u)`：静止段 `(u-u_s)^{4/3}`（Prop.2）过渡到 c-ZOH 二次段（Prop.1）；② `b(u)`：
  静止段 `(u-u_s)^{1/3}`、尾部线性；③ 参数 jerk `⃛u=c√a`：头部 `N_s` 段恒定、尾部
  逐区间且网格点间断。直观印证 interp.py 的 Prop.1 + Prop.2 插值数学。
  示意参数（恒定参数 jerk 宽度 `n_stat`=N_s、`c_tail` 等）在全局参数文件
  [`configs/comm_paras.yaml`](./configs/comm_paras.yaml) 的 `fig4_example` 节，
  经 `load_fig4_example()` 读取（`load_comm_paras()` 取全量）。

## M4 增量（约束扩展，已完成）

在 M1 内核上新增两类约束（`constraints/ingest.py` + `types.TcpConstraint/TorqueConstraint`）：

- **TCP 速度模长**：位置速度模 `‖ṗ‖=cv·√a`、姿态角速度模 `‖ω‖=cw·√a`，均为 a 的
  **线性上界**，折进 `velocity_upper_bound` 的 ā（设计 v0.3 的两项 TCP 约束）。
- **关节力矩**：`τ = n_tor·a + m_tor·b + g_tor`（2 阶、对 (a,b) 精确线性），逐轴
  `τ_min ≤ τ ≤ τ_max`，加进 LP 与种子（论文 eq.44 / robot6dof §5.2.5）。
- 均经 `RobotLimits`（`v_tcp_max/w_tcp_max/tau_max`）+ `to_topp3_data(tcp_geom=…, torque_coeffs=…)`
  一处配置。测试中 TCP 系数由 `UR5RobotModel.tcp_coeffs`（真实 UR5 Jacobian 沿合成关节
  路径求值）给出，力矩系数由 `UR5RobotModel.torque_coeffs`（对角近似，见 `robot/ur5.py`）
  给出；真实逆动力学（M2 的 `DynamicsModel`/RNE）待后续落地。
- 验证：单一测试用例 `test_splp_kernel`（`_assert_solver` 内）断言 TCP 速度模与力矩满足约束；
  `self-test/output/splp_limits_test.png` 面板 ④关节力矩（贴 ±τ_max）、⑤‖ṗ‖（贴住 v_max）
  直观展示约束绑定（各面板限位线均按关节自身实际配置值画，同色虚线，见 `copp/viz.py`）。

## 静止段 / 零进给奇异（Box I，已在优化器落地）

支持真正的 **rest-to-rest**（`a_bnd=(0,0)`，首末速度严格为 0），且**求解/插值/可视化三处自洽**：

- 零进给端的头/尾各 `N_s`/`N_f` 个区间按 **Box I / Proposition 2**（论文式 20）以
  **jerk 零阶保持**离散：段内 `⃛u≡κ` 恒定，`a_k=((u_k-u_s)/(u_1-u_s))^{4/3}a_1`、
  `b_k=2a_k/(3(u_k-u_s))`、`c_k=2a_k/(9(u_k-u_s)²)`——整段由**单一自由变量** `a_1`
  线性决定（`solve/state.py: static_relations`）。
- `N_s` 由 `SolveOptions.n_stationary`（默认 1）与边界共同决定
  （`resolve_num_stationary`；仅 `a_bnd≈0` 端启用，非零端向后兼容走 c-ZOH）。
  jerk 约束在静止点用式 20 点值 `c(u_k)`（`linearize.jerk_constraints` 点式，跳过 rest 端）；
  `interp.py` 把整段静止按 jerk-ZOH 闭式积分（时长 `3(u_{N_s}-u_s)/√a_{N_s}` 有限）；
  `viz.reconstruct_grid_signals` 在 rest 点用 `⃛u=κ` 补正关节 jerk。
- 结果：静止段的 `(a,b,c)` 与 jerk-ZOH 自洽，**关节 jerk 边界值非零 = `q'·κ`**
  （与论文 Fig.4 一致，见 `fig4_interpolation.png` 对照）；零进给奇异根除，`t_final`
  有限且更接近时间最优（c-ZOH 边界的过约束被消除）。

**解析插值（Prop.1 闭式）**：`interp.py` 非静止区间时长用 **Prop.1 闭式 Φ_k**
（`_quad_time`，论文 eq.11，按 c>0/c<0/c=0 分对数/反正弦/根式），已与梯形对拍到
~1e-11。静止段用 Prop.2 闭式。故 `t(s)` 全程解析；`t→s` 反演在解析前向表上插值
（数值等价 Φ_k⁻¹）。**信号重构**（`reconstruct_time_signals`）在区间内闭式细剖面上做，
`q̇/q̈/q⃛` 导数自洽——c-ZOH↔jerk-ZOH 拼接处 `a` 精确连续（`profile.c` 保持纯 c-ZOH
控制、不被静止点值覆写；曾因覆写导致该处加速度突变，已修）。

## 仍待补齐（后续里程碑）

- **备选 SOCP 路线**（`mode="socp"`，精确目标）未实现，默认走 PLP+LP（算力最省）。
- **M2（指令层 + 降维层）已完成**：`trajectory-planning/path/{commands,lowering}` 可跑，
  见 [`README_M2.md`](./README_M2.md)。
- `trajectory-planning/path/blending`（M3）、`trajectory-planning/planner/hlaw`（M5）、
  `trajectory-planning/planner/{synth,planner.py}`（M2+）仍是占位文件，尚无实现，
  见各 `__init__.py` 说明。
- **真实动力学**（`robot.base.DynamicsModel`/RNE）未实现；`UR5RobotModel.torque_coeffs`
  仍是对角近似，见上文"机器人本体"说明。

## 最小用法

```python
import numpy as np
from copp import Topp3Data, solve_splp, SolveOptions

data = Topp3Data(s_grid=..., dq=..., ddq=..., dddq=...,   # q',q'',q''' 于网格
                 vmax=..., amax=..., jmax=...,
                 a_bnd=(0.04, 0.04), b_bnd=(0.0, 0.0))
profile, hist = solve_splp(data, SolveOptions(n_iter=5, verbose=True))
# profile.a / profile.b / profile.c ; hist.t_final（每次迭代的终止时间）
```
