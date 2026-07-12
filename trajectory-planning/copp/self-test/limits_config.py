"""测试用约束配置（唯一来源）。

机器人本体约束（逐关节 vmax/amax/jmax/力矩 + t–n 反电动势/摩擦 + 边界）全部从
YAML 配置文件读取（configs/robot_ur5.yaml，改本体约束改那里）；TCP 速度模上界
作为“给定”参数在此传入（任务/工艺侧设定，不在机器人配置文件里）。

速度相关力矩（t–n）的物理参数（emf_slope/viscous/coulomb）已随 robot_ur5.yaml
逐关节读入 RobotLimits 的 st_* 字段；约束开关在 comm_paras.yaml 的
`constraints.speed_torque`。
"""

from __future__ import annotations

from copp import load_robot_limits

# TCP 速度模上界（给定，非机器人本体参数）
V_TCP_MAX = 1.0
W_TCP_MAX = 5.0

# 其余参数（逐关节运动学/力矩/摩擦限制 + 边界）直接读配置文件
LIMITS = load_robot_limits(v_tcp_max=V_TCP_MAX, w_tcp_max=W_TCP_MAX)
