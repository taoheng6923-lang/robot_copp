"""copp 求解层（SPLP 核心，framework §5.6 / 设计 §7）。"""

from .splp import solve_splp, SolveOptions
from .interp import s_to_t, t_to_s

__all__ = ["solve_splp", "SolveOptions", "s_to_t", "t_to_s"]
