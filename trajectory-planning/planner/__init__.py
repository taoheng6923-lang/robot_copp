"""调度/合成层（M2+/M5，尚未实现）。

编排 `trajectory-planning/path/` → 顶层 `robot/` → `trajectory-planning/copp/`
全流程并产出可执行轨迹：短序列直接调 `copp.solve.solve_splp`，长序列/流式场景
经 `hlaw/` 三窗（种子/可行/最优）分层调度以保证逐窗可行（论文 Theorem 3）；
`synth/` 把求解结果 Profile(a,b,c) 结合路径导数合成时间域轨迹并校验约束满足度；
`planner.py` 是串起以上全部的用户入口门面。见 docs/python_framework.md
§5.7-5.8/§5.10、docs/robot_copp_design.md §8-10。

- hlaw/       长序列分层前瞻窗口调度
- synth/      轨迹合成（Profile → q(t),q̇,q̈,q⃛）与约束校验
- planner.py  TrajectoryPlanner 门面
"""
