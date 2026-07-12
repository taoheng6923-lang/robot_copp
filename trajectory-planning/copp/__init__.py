"""copp —— TOTP-SPLP 时间最优求解数值核心（robot_copp 项目）。

只负责"给定路径导数 + 约束 → 时间最优 (a,b,c) 剖面"这一层（论文
TOTP-SPLP：分段线性目标 PLP + 序列线性化 + LP），不依赖顶层 `robot/`
（运动学动力学本体）与同级 `path/`（路径构造）。两者经
`trajectory-planning/planner/planner.py` 与 copp 编排串联。当前进度：
M1 数值内核 + M4 约束扩展（TCP 速度、关节力矩）。详见 python_framework.md
§9 里程碑与 robot_copp_design.md §7。
"""

from .types import Topp3Data, Profile, TcpConstraint, TorqueConstraint, SpeedTorqueConstraint
from .constraints import RobotLimits
from .options import ConstraintFlags
from .config import (
    load_robot_limits, load_comm_paras, load_constraint_flags, load_smooth_c_weight,
    DEFAULT_CONFIG, DEFAULT_COMM_CONFIG,
)
from .solve.splp import solve_splp, SolveOptions

__all__ = [
    "Topp3Data", "Profile", "TcpConstraint", "TorqueConstraint", "SpeedTorqueConstraint",
    "RobotLimits", "ConstraintFlags",
    "load_robot_limits", "load_comm_paras", "load_constraint_flags", "load_smooth_c_weight",
    "DEFAULT_CONFIG", "DEFAULT_COMM_CONFIG",
    "solve_splp", "SolveOptions",
]
