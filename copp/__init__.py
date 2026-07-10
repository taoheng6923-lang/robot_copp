"""copp —— 基于 copp 算法（TOTP-SPLP）的机器人时间最优轨迹规划包（robot_copp 项目）。

当前进度：M1 数值内核（论文 TOTP-SPLP：分段线性目标 PLP + 序列线性化 + LP）。
详见 python_framework.md §9 里程碑与 robot_copp_design.md §7。
"""

from .types import Topp3Data, Profile, TcpConstraint, TorqueConstraint
from .limits import RobotLimits
from .flags import ConstraintFlags
from .config import (
    load_robot_limits, load_comm_paras, load_fig4_example, load_constraint_flags,
    DEFAULT_CONFIG, DEFAULT_COMM_CONFIG,
)
from .robot import KinematicsModel, DynamicsModel, SyntheticRobotModel
from .solve.splp import solve_splp, SolveOptions

__all__ = [
    "Topp3Data", "Profile", "TcpConstraint", "TorqueConstraint",
    "RobotLimits", "ConstraintFlags",
    "load_robot_limits", "load_comm_paras", "load_fig4_example", "load_constraint_flags",
    "DEFAULT_CONFIG", "DEFAULT_COMM_CONFIG",
    "KinematicsModel", "DynamicsModel", "SyntheticRobotModel",
    "solve_splp", "SolveOptions",
]
