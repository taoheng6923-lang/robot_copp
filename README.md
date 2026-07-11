# robot_copp

基于论文 **TOTP-SPLP**（分段线性目标 + 序列线性化 + LP）的机器人时间最优轨迹规划模块：给定几何路径与运动学/动力学约束（速度/加速度/jerk/力矩/TCP 速度），求满足约束、时间最优的时间律 `s(t)`，进而重构关节轨迹 `q(t), q̇(t), q̈(t), q⃛(t)`。

## 当前状态

**M1 数值内核 + M4 约束扩展 + M2 指令/降维层已实现并可跑**：指令序列（关节运动 /
直线 / 圆弧）→ 笛卡尔路径 → 真实 UR5（DH 运动学、解析 IK）连续解降维 → SPLP
时间最优求解，全链路打通（M2 语义：段间 G0 衔接、逐段 rest-to-rest）。路径
blending（M3）/ planner 门面 / HLAW 在线窗口调度（M5）尚未实现（占位目录已建好）。
详细进度见 [`docs/README_M1.md`](docs/README_M1.md)（内核）与
[`docs/README_M2.md`](docs/README_M2.md)（指令+降维）。

## 快速开始

```powershell
.\.venv\Scripts\python.exe -m pytest -q     # 三个自测：copp 内核 / lowering / commands
# 或单独跑（额外生成分析图到各自 self-test/output/）：
.\.venv\Scripts\python.exe trajectory-planning\copp\self-test\test_splp_kernel.py
.\.venv\Scripts\python.exe trajectory-planning\path\self-test\test_lowering.py
.\.venv\Scripts\python.exe trajectory-planning\path\self-test\test_commands.py
```

## 项目结构

```
robot/                    顶层，独立模块：机器人本体（ur5.py 真实 DH 运动学 + 解析 IK）
trajectory-planning/
├── copp/                 TOTP-SPLP 数值求解核心（已实现，M1+M4）
├── path/                 路径构造（M2 已实现 commands/lowering；blending 属 M3 未实现）
│   ├── commands/         三类运动指令 → Section（几何 + 端点元数据）
│   ├── lowering/         自适应采样 + 连续解 IK + 链式法则求 q',q'',q'''
│   └── blending/         G2 Hermite 过渡（占位，未实现）
└── planner/              调度/合成/门面：hlaw/synth/planner.py（占位，未实现）
configs/                  机器人约束 YAML（robot_ur5.yaml 等）+ 通用参数
docs/                     设计文档、算法参考、里程碑状态（见下）
```

## 文档索引

| 文档 | 内容 | 权威性 |
|------|------|--------|
| [`docs/README_M1.md`](docs/README_M1.md) | **M1/M4 内核实际状态**：已实现模块、机器人参数来源可信度、如何运行、已知待办 | 如实反映当前代码，改动后应第一时间更新 |
| [`docs/README_M2.md`](docs/README_M2.md) | **M2 指令+降维层实际状态**：commands/lowering 模块表、与设计文档的两处已确认修正、验证方法 | 如实反映当前代码，改动后应第一时间更新 |
| [`docs/robot_copp_design.md`](docs/robot_copp_design.md) | 架构设计方案：指令层→blending→降维→约束→求解→HLAW 的完整数据流与公式 | 算法/架构层面权威；模块树（§10）与实际目录结构有出入 |
| [`docs/python_framework.md`](docs/python_framework.md) | 配套 `robot_copp_design.md` 的代码骨架/接口签名参考 | 早期单包结构快照，算法原理仍有效，见文档内差异说明 |
| [`docs/paper_notes.md`](docs/paper_notes.md) | TOTP-SPLP 论文笔记（状态变换、无损离散化、PLP、HLAW 等理论） | 论文推导参考 |
| [`docs/code_reading_guide.md`](docs/code_reading_guide.md) | 原 Rust `copp` 库代码阅读指南（公式↔代码映射） | 上游 Rust 库参考，本仓库为其 Python 复刻/延伸 |
| [`claude-record/SESSION_HANDOFF.md`](claude-record/SESSION_HANDOFF.md) | 项目启动前一次会话的存档（早期 Python 框架设计讨论） | 历史存档，不代表当前设计 |

`CLAUDE.md` 是给 Claude Code 的项目级操作说明（环境注意事项、测试命令），非面向人类读者的文档。
