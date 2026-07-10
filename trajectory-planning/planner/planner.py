"""TrajectoryPlanner 门面（M2+，尚未实现）。

串起全流程用户入口：path.commands → path.blending → path.lowering →
copp.constraints → (hlaw?) copp.solve → synth。见 docs/python_framework.md
§5.10、docs/robot_copp_design.md §10。
"""
