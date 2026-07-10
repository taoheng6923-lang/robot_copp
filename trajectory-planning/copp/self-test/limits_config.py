"""测试用约束配置（唯一来源）。

机器人本体约束（逐关节 vmax/amax/jmax/力矩 + 边界）从 YAML 配置文件读取
（configs/robot_3axis.yaml，改本体约束改那里）；TCP 速度模上界作为“给定”
参数在此传入（任务/工艺侧设定，不在机器人配置文件里）。
"""

from __future__ import annotations

from copp import load_robot_limits

# TCP 速度模上界（给定，非机器人本体参数）
V_TCP_MAX = 1.0
W_TCP_MAX = 5.0

# 其余参数（逐关节运动学/力矩限制 + 边界）直接读配置文件
LIMITS = load_robot_limits(v_tcp_max=V_TCP_MAX, w_tcp_max=W_TCP_MAX)
