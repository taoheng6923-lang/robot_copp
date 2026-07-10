"""求解器后端抽象（framework §5.9）。M1 提供 cvxpy 后端。"""

from .cvxpy_backend import solve_problem, DEFAULT_SOLVER

__all__ = ["solve_problem", "DEFAULT_SOLVER"]
