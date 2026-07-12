"""约束摄入层（framework §5.5 / 设计 §6）。

- model:  RobotLimits（机器人本体约束配置，M1/M4 的 ConstraintSet 前身）
- ingest: 物理约束 → 路径域不等式（TCP 速度模上界、关节力矩行）
"""

from .model import RobotLimits
from .ingest import (
    tcp_a_upper_bound, torque_constraints,
    speed_torque_constraints, speed_torque_utilization, speed_torque_envelope,
)

__all__ = [
    "RobotLimits", "tcp_a_upper_bound", "torque_constraints",
    "speed_torque_constraints", "speed_torque_utilization", "speed_torque_envelope",
]
