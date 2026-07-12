# Session Handoff — robot_copp 项目会话记录

> 本文件累积记录本仓库（`robot_copp`）关键会话的目标、决策、产出与未完成事项，供跨会话/换环境时续接工作。**按时间倒序排列，最新会话在最前面**。Session 1（下方最后一节）是项目启动前从另一仓库导出的设计讨论存档，与当前代码结构/实现已有很大出入，仅供历史参考——真正"如实反映当前代码"的权威文档是 [`docs/README_M1.md`](../docs/README_M1.md)（内核 M1/M4）与 [`docs/README_M2.md`](../docs/README_M2.md)（指令+降维 M2）。

---

## Session 3（2026-07-13）：依赖修复可跑 + UR5 解析 IK + t–n 约束调通 + Fig.3/Fig.4 出图 + c 平滑惩罚 + 设计文档补全

### 0. 背景与定位

接手时项目**依赖装不上（测试跑不了）**，且 Session 2 搭的 t–n 速度相关力矩约束**未调通**。本会话主线是"把项目真正跑起来 + 补齐 copp 求解层的 t–n 约束、可视化与设计文档"。改动集中在 `trajectory-planning/copp/`、`robot/ur5.py`、`configs/`、`docs/`。**收尾 `pytest -q` = 3/3 通过**（test_copp / test_lowering / test_commands）。

> 与 Session 2 handoff 的关系：Session 2 把"解析 IK / t–n 约束 / test_copp 改名"记为"用户或 linter 修改、未亲见"。**本会话是实际动手实现/调通这些的过程**（尤其 t–n 从不可行→可行、Fig.4 从合成示意→真实用例数据）；下文以本会话实测为准。

### 1. 依赖修复（Windows 段错误）

`cvxpy`/`clarabel` 缺失，直接 `pip install` 会把 numpy 升到 2.x 并**段错误崩溃**（`import cvxpy` → 0xC0000005）。逐个二分定位到 `qdldl 0.1.9`、`osqp 1.1.3` 的最新 Windows 轮子 import 即崩。最终方案：

- 独立虚拟环境 `.venv`（`.\.venv\Scripts\python.exe`，Anaconda Python 3.12 建），**钉稳定栈**：`numpy 1.26.4`、`cvxpy 1.5.4`、`osqp 0.6.7.post3`、`qdldl 0.1.7.post5`、`clarabel 0.11.1`、`scipy 1.13`、`matplotlib`。
- 把这些上/下界写进 `pyproject.toml` `dependencies`（`numpy<2`、`cvxpy<1.6`、`osqp>=0.6.7,<1.0`、`qdldl>=0.1.7,<0.1.8` 等），避免下次 fresh install 再踩崩溃轮子。
- **恢复 Anaconda base**：误升级把 base 的 numpy 弄成 2.x（破坏 base scipy/matplotlib），已还原到 1.26.4 并清掉误装包。
- **CLAUDE.md 已记**：不要用 base 环境装本项目依赖；测试用 `.\.venv\Scripts\python.exe -m pytest -q`。

### 2. CLAUDE.md 中文回复规则

按用户要求在 `CLAUDE.md` 加"始终用中文回复"（项目级持久生效），顺带记了 `.venv` 测试命令与依赖注意事项。

### 3. UR5 IK：数值 DLS → 闭式解析（`robot/ur5.py`）

把 `ik` 从阻尼最小二乘牛顿迭代（≤200 步）换成 **Andersen 标准 DH 闭式解析逆解**：枚举 θ1(2)×θ5(2)×θ3(2)=≤8 支路，逐关节 2π 折叠后取离 seed 最近解，DLS 降为退化位形兜底（`_ik_dls`）。

- **过程中发现并修的真实 bug**：初版平面两连杆（θ2/θ3/θ4）用错了平面——`T14=A2·A3·A4`（α2=α3=0）原点退化的是 **x-y 平面**（z≡d4 常量），我一开始按 x-z 平面解，导致 8 支路都有 ~0.1–0.2m 位置误差（姿态却对，因为 θ4 由 `T34` 反解总能补正姿态）。改用 x-y 平面 + `θ234=atan2(R14[1,0],R14[0,0])` 后往返 `fk∘ik` 误差降到 **~1.2e-13 m / 3e-8 rad**（3000 组随机位形），比 DLS 快 **~5.9×**。

### 4. t–n 速度相关力矩约束：从不可行调到可行（本会话最大修复量）

Session 2 已搭好类型/接线（`types.SpeedTorqueConstraint`、`ingest.speed_torque_constraints`、`ConstraintFlags.speed_torque`、lp_problem/seed/state/model/config 接线），但**种子 LP 不可行、跑不动**。本会话逐个修：

1. **cvxpy `*` 被当矩阵乘**：`sec = a * (aq/w0)`（`Variable(N)*ndarray(N)` → 点积标量，污染约束）→ 改 `cp.multiply`。
2. **静止端不可行（核心）**：原"固定切点"切线在 a=0 处引入伪速度项 `≈Fv·q_j/2`；UR5 肩关节重力 g≈140 已占 τ0=150 的 93%、余量仅 ~5N·m，被伪项吃掉 → 种子不可行。**改为把 √a 在 SCP 迭代点 a_lin 处线性化**（切点随迭代收敛到工作点、静止端无伪项），并把梯形包络写成"平台 + rolloff **两个精确 halfplane 之交**"（无需 facet 逼近梯形本身），跳过静止段（`num_stat`）。seed 在速度上界 a_bar 处线性化。
3. **模型演进**：从"线性 τ0−κ|q̇|（emf_slope）"改为真实电机数据表的**梯形** t–n 特性（低速平台 τ0 到拐点 ω_c=vmax，线性收窄到 0 于空载 ω0）；`SpeedTorqueConstraint` 字段变为 `tau0/rated_speed/noload_speed/viscous/coulomb`（不再有 `emf_slope`）。
4. **viz 用了废弃的 `st.emf_slope`** → 改用 `speed_torque_envelope`（梯形）；`plot_speed_torque` 里 τ0−κ|q̇| 标注也改梯形。
5. **合成 Fv 过大**（肩关节额定转速黏滞 27N·m，物理不合理，导致任何速度都超力矩）→ 降到 `≈0.015·τ0/ω_c`（1.5%·τ0）。
6. **验证**：解回代真实梯形 util ≤ **1.0000**（保守、贴边、绝不越界）；关掉 t–n 时真实 util 达 **~11.5**（必超限，证明约束必要）。启用 t–n 时轴速上界自动取 ω0（见 `state.velocity_upper_bound`）。
7. `SpeedTorqueConstraint.n_facets` 字段现已**无用**（实现固定用两个精确 halfplane），保留未删。

### 5. `tau_scale` 力矩倍率（可配置实验旋钮）

把 `robot_ur5.yaml` 的 `tau_max` 恢复为**官方值**，新增顶层 `tau_scale`（`config.load_robot_limits` 统一乘 tau_max/tau_min、进而 t–n 平台 τ0），实验放开力矩余量看效果：1×→肩重力饱和限速（+226% 时间）、2×→+48%、3×→+26%（瓶颈转到 base 惯量）。**当前用户手动设 `tau_scale: 5.0`**。

### 6. 论文 Fig.3 复现（`viz.plot_tn_convexification`）

逐关节在 `(q̇²,τ)` 平面画：灰=真实梯形可行域（rolloff 段映射后弯曲→下方非凸）、蓝=仿射切线内逼近、红叉纹="切掉的非凸角"、散点=工作点（点色=利用率）。直观显示论文"切一角转凸"的意图与效果。出 `output/tn_convexification_test.png`。

### 7. Fig.4 重构 + `test_splp_kernel.py → test_copp.py` + 唯一用例出图

- **加控制量 c 曲线**：Fig.4 从 3 面板（a/b/jerk）变 4 面板（a/b/**c**/jerk）。c=b'=⃛u/√a：非静止段 c-ZOH 逐区间常值阶梯，静止端 →∞（正是静止段用 jerk-ZOH 的原因）。
- **改为由唯一用例 profile 全程解析重构**：删掉合成示意数据路径（`fig4_interpolation_example`、`plot_fig4_interpolation`、`config.load_fig4_example`、`comm_paras.yaml` 的 `fig4_example` 节全部移除），新增 `viz.plot_interp_profiles(data, profile)` 用 `interp.fine_profiles` 覆盖**全部插补周期**（原来只有 4 个合成周期）。
- **测试改名**：`test_splp_kernel.py → test_copp.py`（`git mv`），函数 `test_splp_kernel → test_copp`；**output 下 5 张图全部来自唯一用例的同一 `(data, profile)`**（splp / limits / speed_torque / tn_convexification / fig4）。

### 8. c 控制平滑惩罚（目标函数扩展）

中段 c 会无谓来回跳变（对时间最优无必要、使 jerk 锯齿）。加 **L1 平滑惩罚**：仅对相邻**非静止**区间对 (k,k+1) 引入 `p_k≥|c_k−c_{k+1}|`（两条 LP 不等式，`c_k=(b_k−b_{k-1})/Δ_k` 仿射于 b），目标 `+= λ·Σp_k`。

- `SolveOptions.smooth_c_weight`（默认 0）+ `config.load_smooth_c_weight` 读 `comm_paras.yaml` 的 `objective.smooth_c_weight`。**用户最终设 0.005**（我建议的轻量值；0.05 太强）。
- 实测：λ=0.05 → 中段 Σ|Δc| 降 ~99%（6053→52）、t_final +20.8%；λ=0.005 属轻量（时间代价数个百分点）。**仅非静止段参与**（静止段 c 发散、用 jerk-ZOH）。保持纯 LP、不影响可行性/收敛证明。

### 9. 设计文档补全（`docs/robot_copp_design.md`）

- **§7.2⑤**：新增 c 平滑惩罚的公式与说明（LP 不等式 + 目标项 + 实测权重表 + 落地位置）。
- **§7.0（新）"本仓库实际求解的完整离散 LP（权威形式）"**：一次性写全**决策变量 / 目标函数 / 全部约束**（边界、非静止动力学、静止 Box I、速度上界、轴向加速度、jerk 线性化、力矩盒式、t–n、PLP 割线/下界、c 平滑），每条标注 `ConstraintFlags` 开关与论文式号，逐条对应代码。
- **§7.0 t–n 描述更正**（用户指出有误）：从含糊的"√a SCP 线性化"改为正确归属**论文（Ardeshiri 2010）凸化方法**——核心是改写到 `(τ,q̇²)` 平面（因 `q̇²=q'²·a` 线性于 a，故约束仿射于 a，论文式 17→18），√a 的 SOCP/切线处理只是 LP 落地细节；引用 `tn_constraint_notes.md §5/§8/§11`。
- 另：`tn_constraint_notes.md` 本身在本会话创建并扩充（§8 粘滞/库仑摩擦如何并入：Fv 进斜率、Fc 进截距；§11 落地实现）；`README_M1.md` 同步更新（七类约束、t–n 节、五张图、test_copp、lp_problem 行）。

### 10. 当前测试与 git 状态

- `pytest -q`（testpaths = copp/self-test + path/self-test）→ **3/3 通过**。`planner/self-test/test_planner.py` 仍是 Session 2 遗留的未通过项（未收进 testpaths），本会话未碰。
- 本会话**未 commit**（遵循"仅用户明确要求才提交"）。`git status` 有大量修改 + `test_copp.py`（由 `git mv` 改名）。
- 用户手动改了两处配置：`robot_ur5.yaml` `tau_scale: 5.0`、`comm_paras.yaml` `objective.smooth_c_weight: 0.005`（均属正常调参，勿回退）。

### 11. 交接提醒

- t–n 约束的 √a 线性化用**切线上界(系数≥0)/过原点割线下界(系数<0)**，跳过静止段；若改走 SOCP（`mode='socp'` 预留、未实现）则 √a 可精确、退化为论文原式。深改前先读 `constraints/ingest.py::speed_torque_constraints` + `docs/tn_constraint_notes.md`。
- Fig.4/所有 output 图**必须来自唯一用例** `(data, profile)`——不要再引入单独示意数据（用户明确要求唯一性）。
- 合成参数（`torque_coeffs` 对角近似、t–n 的 Fv/Fc/ω0、amax/jmax、tau_scale）都不是真实 UR5 规格，代码/YAML/文档均有可信度标注。

---

## Session 2（2026-07-12 ～ 2026-07-13）：git 落地 + UR5 真实机器人 + M2 指令/降维层 + planner 门面（进行中）

### 0. 背景

Session 1 的导出文档来自另一个仓库（Rust `copp` 库项目），当时只产出了设计文档、无可运行代码。`robot_copp` 是**独立的新仓库**：Session 1 的 handoff 被复制进来做背景参考后，项目在这里独立演进，M1 数值内核（TOTP-SPLP）+ M4 约束扩展（TCP 速度、关节力矩）已经先于本次会话被实现并可跑（合成的 3 轴 stand-in 机器人）。本次会话开始时仓库还**没有接入 git**。

### 1. Git 落地

初始化仓库、创建首个 commit、推送到 `https://github.com/taoheng6923-lang/robot_copp`（用户账号 `taoheng6923-lang`）。此环境本机没有装 Git，先用 `winget install Git.Git` 装好，再 `git init` + 设置本地 `user.name/user.email`（`taoheng6923` / `taoheng6923@users.noreply.github.com`）+ 提交 + 推送（首次推送经 Git Credential Manager 弹窗登录）。

### 2. 顶层文件夹层级：三轮重构

**第一轮**（在 `copp/` 包内部按 `python_framework.md` §2 的目标结构对齐）：把原来平铺在 `copp/` 下的 `limits.py`+`constraints.py` 拆成 `constraints/`（`model.py`=`RobotLimits`、`ingest.py`=TCP/力矩摄入函数）子包；`flags.py` 改名 `options.py`（`ConstraintFlags`）；为尚未实现的里程碑模块（`commands/`、`blending/`、`lowering/`、`hlaw/`、`synth/`、`planner.py`）建了空的占位目录/文件（只有 docstring，无实现代码——用户当时明确要求"只设计文件夹层级，不需要具体实现"）。

**第二轮**（用户要求把"路径构造"和"copp 核心"拆成不同顶层模块，新增 controller 模块，robot 独立成模块）：

- 新建 `trajectory-planning/`（**纯目录容器，不是 Python 包**——无 `__init__.py`，且目录名含连字符不能作为合法模块名）
- `copp/`（原顶层）→ 移入 `trajectory-planning/copp/`；`copp/__init__.py` 去掉了对 robot 的 re-export（保持数值核心不反向依赖机器人本体）
- 原 `copp/commands`、`copp/blending`、`copp/lowering` 三个占位包 → 移到新建的顶层 `pathgen/`（用户又要求改名为 `path/`，二者都放进 `trajectory-planning/`）
- 原 `copp/hlaw`、`copp/synth`、`copp/planner.py` → 移到新建的顶层 `controller/`（用户要求改名为 `planner/`，同样放进 `trajectory-planning/`）
- `copp/robot/` → 独立成顶层 `robot/`（不再嵌套在 copp 下）
- 每次移动都同步改了内部 import（尤其 `sys.path` 拼接的相对层数）、`pyproject.toml` 的 `packages`/`testpaths`，并加了 `pythonpath = [".", "trajectory-planning"]` 让 `path`/`copp`/`planner` 三个顶层包和 `robot`（不在 trajectory-planning 下）都能被找到

**第三轮**（测试目录规范化）：`tests/` 挪到 `copp/self-test/`（改名 `self-test`）；生成的分析图目录几经调整，最终定为各自 `self-test/output/`（例如 `trajectory-planning/copp/self-test/output/`），而不是仓库级的 `output/` 或 `test2/`。

三轮重构后的最终结构：

```
robot/                          顶层，独立：机器人本体（ur5.py）
trajectory-planning/            纯目录容器（非 python 包）
├── copp/                       纯 TOTP-SPLP 数值核心（M1+M4，已实现）
│   └── self-test/              test_copp.py（原 test_splp_kernel.py）+ output/
├── path/                       路径构造（M2 已实现 commands+lowering；blending=M3 占位）
│   └── self-test/              test_lowering.py、test_commands.py + output/
└── planner/                    调度/合成/门面（本次会话进行中；hlaw=M5 占位）
    └── self-test/              test_planner.py（未完全通过，见下）+ output/
```

### 3. 机器人模型换成真实 UR5（`robot/ur5.py`）

删除了旧的 `SyntheticRobotModel`（与关节数无关的合成 3 轴 stand-in），新增：

- **`UR5Kinematics`**：标准 DH 参数表（Universal Robots 官网权威数据，已用 WebFetch 核实）。`fk`/`jacobian` 是解析解（用有限差分交叉验证到 ~1e-10）。`ik` 最初实现是阻尼最小二乘（DLS）数值迭代，后来被扩展/替换成**闭式解析逆解**（Andersen 标准 DH 推导，θ1×θ5×θ3 枚举全部 ≤8 支路，逐关节 2π 折叠后取离 seed 最近解，O(1) 求解，比 DLS 快 1~2 个数量级），DLS 降级为退化位形（腕奇异/不可达）的兜底——已用 200 组随机位形数值验证到 ~1e-13。**注意**：这部分解析 IK 代码是在系统提醒里以"用户或 linter 修改"的形式出现的，我没有亲眼看到实现过程，接手前建议自己跑一遍验证（`robot/ur5.py` 里的 `_ik_analytic` 方法）。
- **`UR5RobotModel`**：`joint_path` 仍是合成随机轨迹（占位，因为路径生成层此时还没做出来），但 `tcp_geometry`/`tcp_coeffs` 已改用真实 `UR5Kinematics.jacobian(q)` 沿该合成路径求值（不再是与关节角无关的虚构公式）；`torque_coeffs` 是"下游集中质量单摆臂"近似（基于 ROS-Industrial `ur5.urdf.xacro` 公开的连杆质量/长度），不是精确 RNE。
- **逐关节限值数据来源**（在 `robot/ur5.py` docstring 和 `configs/robot_ur5.yaml` 注释里都做了三档可信度标注，避免被误当真实规格）：
  - `tau_max=[54,150,150,28,28,9]` Nm —— **官方权威**（Universal Robots 官网 "Max. joint torques CB3 and e-Series"）
  - `vmax=[3.15,3.15,3.15,3.2,3.2,3.2]` rad/s —— 来自 ROS-Industrial `ur5.urdf.xacro` 的 `<limit velocity=.../>`（驱动真实 UR5 硬件的社区描述文件，可信但非 UR 官方数据表原文；一开始用的是错误的"统一 180°/s 近似"，后来查证后改成了逐关节的真实值）
  - `amax`/`jmax` —— **无任何官方或社区公开数据**，按 vmax 的经验比例臆造，仅供数值示例

### 4. `copp/viz.py` 限位线修正

原来三处约束线画的是"跨关节 max/min 包络"（比如力矩面板画一条 ±150Nm 的线），对 UR5 这种逐关节力矩上限差异巨大（9~150Nm）的真实机器人是**误导性的**——wrist3 实际上限只有 9Nm，画在 150Nm 包络线下会显得还有巨大裕度。改成逐关节按自己的配置值画、与该关节曲线同色的虚线。已用数值检查确认（关节 1 精确顶到 -150.00Nm 绑定，不是越界）。

### 5. 文档整理

填了此前完全空白的根 `README.md`（项目一句话定位 + 当前状态 + 快速开始 + 目录结构 + **文档索引表**：每份文档是什么、权威性如何——这是之前最缺的东西，6+ 份文档之间没有任何导航）；修了 `CLAUDE.md` 里已经不存在的 `pytest tests/ -q` 命令；`docs/README_M1.md`、`docs/robot_copp_design.md`、`docs/python_framework.md` 都补了"与实际实现的差异"提示框，指向真正如实反映代码的文档。

### 6. M2 层实现：`path/lowering/` + `path/commands/`（本次会话主要工作量）

**`path/types.py`**：`CartesianSamples`（位姿 + 1~3 阶导采样）、`CartesianPath`/`JointSpacePath` protocol、`PathDerivatives`（q/dq/ddq/dddq + cv/cw，可直接喂 `RobotLimits.to_topp3_data`）。

**`path/lowering/`**：
- `sampling.py` —— 曲率驱动自适应采样（弦高误差模型 e≈‖r''‖Δs²/8，网格强制含 `s_breaks`）
- `ik.py` —— 连续解 IK（seed 链）+ FK 回代校验 + 相邻站点跳变校验
- `derivatives.py` —— Jacobian 链式法则求 q'/q''/q'''（J'/J'' 用方向有限差分）。**修正了设计文档 `robot_copp_design.md` §5.3 的一个公式错误**：三阶式漏写了系数 2，正确式是 `Jq''' = r₃ − 2·J'q'' − J''q'`（对 `Jq'=r₁` 逐阶求导两次得到）。这个错误如果照抄设计文档实现，在自测里会体现为 O(10%) 量级的三阶导交叉验证失败——已在代码注释和 `docs/README_M2.md` 里记录。
- `singularity.py` —— σ_min/σ_max 奇异检测 + 阻尼最小二乘逆

**`path/commands/`**：`JointMoveCommand`（**有意偏离设计文档**：用线性关节几何而非设计文档 §3.1 的五次多项式——理由是 TOPP 框架下几何与时间律解耦，时间平滑由 copp 时间律保证，五次多项式反而会让端点 q'(s)=0 使参数化退化）、`LinearMoveCommand`（直线+SLERP，纯姿态调整退化安全）、`CircularMoveCommand`（三点定圆/圆心+法向，弧长参数化）、`assemble.py`（`build_sections`+`lower_sections`：G0 精确衔接校验 + IK seed 链逐段降维）。M2 语义：**段间无 blending**，角点用停顿衔接，每段独立 rest-to-rest 求解（G2 平滑过渡是 M3）。

每层都写了高强度自测（`test_lowering.py`、`test_commands.py`），核心原则是**不相信被测代码自己的输出**，大量用独立手段交叉验证：解析测试路径先用有限差分自检推导对不对、q 序列有限差分反推 dq/ddq/dddq、独立网格差分验证 Jacobian 恒等式 J'q'+Jq''≈r₂、退化情形（纯姿态调整 p'≡0、共线圆弧、零长指令）显式报错。

`docs/README_M2.md` 已写好，如实记录了上述内容。

### 7. t–n（转矩–转速）速度相关力矩约束

这部分是**用户直接推进补充的**（我在系统提醒里看到大量"文件被用户或 linter 修改"的通知，没有亲眼见证实现过程，只读到了最终结果）：

- `docs/tn_constraint_notes.md`：精读 Ardeshiri et al. 2010 论文（*Convex Optimization approach for Time-Optimal Path Tracking of Robots with Speed Dependent Constraints*），梳理"电机转矩上限随转速下降"这类速度相关约束如何在不破坏凸性的前提下并入 TOPP 凸优化——核心是把约束改写成仿射于 `(τ, q̇²)` 而不是仿射于 `q̇` 本身。**这正是用户中途打断时问的"非凸约束切一角转成凸约束的参数设置"**：约束原式含 √a（非凸，a 是路径速度平方），按 SPLP 一贯的序列线性化思路，在当前迭代点 `a_lin` 处对 √a 取**切线上界**（保守内逼近，"切掉"了非凸可行域超出切线的那一角），随迭代收敛到精确解——与 jerk 约束用的是同一套线性化机制（论文 eq.32 那种手法），只是这次用在速度相关力矩上。
- 落地：`types.SpeedTorqueConstraint` + `constraints/ingest.speed_torque_constraints` + `ConstraintFlags.speed_torque`，约束形式 `|τ_dyn + Fv·q̇ + Fc·sgn(q̇)| ≤ τ0 − κ·|q̇|`（κ=反电动势斜率、Fv=粘滞摩擦、Fc=库仑摩擦，τ0 复用 tau_max）。
- 参数（`emf_slope`/`viscous`/`coulomb`）逐关节配在 `configs/robot_ur5.yaml`，同样标注"合成臆造，无官方数据"（κ=0.45·τ_max/vmax、Fv=0.4κ、Fc=0.03·τ_max 这类经验比例）。
- 原来的 `test_splp_kernel.py` 改名成了 `test_copp.py`。

### 8. planner 门面 + synth 层（本次会话进行中，**未完成**）

这是被用户提问打断时正在做的工作：

- **`planner/synth/resample.py`**：`synthesize`（单段：在 copp 的解析细剖面——`fine_profiles`，Prop.1 c-ZOH + Prop.2 静止段 jerk-ZOH 闭式——上重构 q̇/q̈/q⃛ 到等时间栅格）+ `concatenate`（多段按时间顺序拼接，段间 rest 停顿衔接）+ `TrajectoryResult` 数据类。
- **`planner/synth/verify.py`**：`verify_limits`（超限率 R_v / 超限时长比 D_v，论文 §6.1.2 的验收指标）+ `VerifyMetrics`。
- **`planner/planner.py`**：`TrajectoryPlanner` 门面类——`add_command()` 累积指令（链式调用），`plan(q_seed)` 一次性跑完 `build_sections → lower_sections → 逐段 solve_splp → synthesize → concatenate → verify_limits`。多段规划会检查 `limits.a_bnd/b_bnd` 必须是静止边界（M2 段间停顿语义的前提），否则报错。
- **`planner/self-test/test_planner.py`**：JointMove+LinearMove+CircularMove 三段混合指令端到端自测（结构性质、rest-to-rest 连续性、FD 交叉验证、约束校验器灵敏度、防误用断言）。

**过程中发现并修复的一个真实 bug**：`copp/solve/interp.py` 的 `_segment_tail_static`（静止尾段的细分点采样）原来按 `s` 均匀取点，但静止段 `ṡ ∝ ρ^{2/3}`（ρ 是到 rest 端的距离），导致最后一个细分区间在**时间上**占了整段约 35%——这个非均匀性在等时间栅格插值时会造成明显失真。改成按 `ρ_frac³` 立方分布采样后，`FD(s) 中心差分 vs 解析 ṡ` 的交叉验证从误差 5.2e-2 降到通过阈值。这个 bug 之前没被发现是因为旧的可视化/测试代码都是在**网格点**上重构信号，没有在时间栅格上做过独立的 FD 交叉验证。

**当前卡住、尚未解决的问题**：`test_planner.py` 的 `_assert_verify` 断言"真实限值下应该 R_v=D_v=0"失败，实测 `R_v=0.313`（31% 样本"超限"）。中断前已经用一次性诊断脚本细分到每类约束在不同容差下的超限样本占比，结论是：
- `velocity`、`acceleration` 完全不超（0%）
- `jerk`：0 容差下 7.85% 样本超，但全部在 3.4% 以内（`max_util=1.034`）——这大概率是解析细剖面区间内 O(Δ²) 离散化误差的正常范围，不是实现 bug
- `tcp_velocity`：0 容差下 23.73% 样本超，但全部在 0.06% 以内（`max_util=1.0006`）——是"贴线"而非真超限

**下一步建议先查这个方向**：`verify_limits`/`test_planner.py` 的容差设置可能过严（比如 jerk 用了硬 `1.0` 而没考虑区间内 O(Δ²) 误差界、tcp_velocity 的 0.06% 量级可能只是浮点/插值噪声），而不是去怀疑 `resample.py`/`verify.py` 的实现逻辑有错——已经确认没有真正意义上的大幅超限（最坏情况才超 3.4%）。改完容差后记得把 `trajectory-planning/planner/self-test` 加进 `pyproject.toml` 的 `testpaths`（现在还没加）。

### 9. 当前测试状态

```
pytest -q   # testpaths = ["trajectory-planning/copp/self-test", "trajectory-planning/path/self-test"]
```
→ **3/3 通过**（test_copp、test_lowering、test_commands）。

```
python trajectory-planning/planner/self-test/test_planner.py
```
→ **失败**（`_assert_verify` 断言，见上）。这个测试还没被 `pyproject.toml` 的 `testpaths` 收录，所以 `pytest -q` 目前看不出这个失败。

### 10. Git 状态（本次会话结束时）

最近一次 commit 是 `a818c82 feat: 增加lowering和command层；增加t-n约束及摩擦力约束的文档梳理`（不含 planner 门面这部分工作）。截至本次会话结束，`git status` 显示大量文件已修改，且 `trajectory-planning/planner/synth/resample.py`、`verify.py`、整个 `trajectory-planning/planner/self-test/` 目录都是**未跟踪的新文件**，尚未 commit。

### 11. 关键设计决策速查（避免下次重新讨论）

- `trajectory-planning/` 是纯目录容器，不是 Python 包（含连字符不能作模块名，也没有 `__init__.py`）；`path`/`copp`/`planner` 三个顶层包靠 `pyproject.toml` 的 `[tool.pytest.ini_options] pythonpath = [".", "trajectory-planning"]` 找到彼此，`robot` 靠仓库根（`.`）找到。
- `copp/` 不反向依赖 `robot/`/`path/`——只通过 `Topp3Data`/`PathDerivatives` 这类纯数据结构交互。
- M2 无 blending：段间 G0 精确衔接 + 每段独立 rest-to-rest，物理严格可行但非全局时间最优；G2 不停顿平滑过渡是 M3（`path/blending/`，仍是空占位）。
- IK 用 UR5 闭式解析解（8 支路），非数值迭代法（更快更稳；DLS 只做退化位形兜底）。
- `torque_coeffs`（对角近似）和 t–n 参数（emf_slope/viscous/coulomb）都是"物理量级合理但非真实测量"的合成近似，反复在代码注释/YAML 注释/文档里标注来源可信度，不能当作真实 UR5 出厂规格使用。
- 真实动力学（`DynamicsModel`/RNE）、M3 blending、M5 HLAW 均未实现，仍是占位目录。

### 12. 交接提醒

本次会话中，相当一部分文件是以"该文件被用户或 linter 修改"的系统提醒形式出现变化的（UR5 解析 IK、t–n 约束整块、部分文档措辞），我没有亲眼见证实现/推导过程，只是照最终状态做了描述和交叉验证（跑测试、核对公式）。下一次会话如果要深入改这些部分，建议先重新通读一遍相关源码（尤其 `robot/ur5.py::_ik_analytic`、`copp/constraints/ingest.py::speed_torque_constraints`），不要完全依赖本节的文字复述。

---

## Session 1（2026-07-10 导出，另一仓库）：COPP-Python 轨迹规划框架设计（历史存档）

> ⚠️ 以下为项目启动前从另一仓库（Rust `copp` 库项目）导出的会话记录，**当时只产出设计文档、无可运行代码**。`robot_copp` 项目后来独立演进，模块划分、命名、技术选型都已和下文有很大出入（例如：下文设想的单包 `copp_py/` 结构已被 Session 2 的 `robot/`+`trajectory-planning/{copp,path,planner}` 四包结构取代；下文的 `PathSegment`/`PathBuilder` 设计已被 `path/types.py` 的 `CartesianSamples`/`Section` 取代）。仅供了解最初的问题背景与算法参考，**不代表当前设计**——当前设计权威文档见 [`docs/robot_copp_design.md`](../docs/robot_copp_design.md)、[`docs/README_M1.md`](../docs/README_M1.md)、[`docs/README_M2.md`](../docs/README_M2.md)。

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
