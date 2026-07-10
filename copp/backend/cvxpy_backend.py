"""cvxpy 求解后端（framework §5.9）。

把"组装好的 cvxpy.Problem → 数值解"与建模解耦。默认用 CLARABEL
（与 Rust copp 同核，便于数值对照），缺失时回退到其它已装 LP/SOCP 求解器。
"""

from __future__ import annotations

import cvxpy as cp

# 优先级：CLARABEL（与 Rust copp 一致）> HIGHS/GLPK（纯 LP）> SCS/ECOS（兜底）
_PREFERENCE = ["CLARABEL", "HIGHS", "GLPK", "ECOS", "SCS"]


def _pick_solver() -> str:
    installed = set(cp.installed_solvers())
    for name in _PREFERENCE:
        if name in installed:
            return name
    raise RuntimeError("未找到可用的 cvxpy 求解器")


DEFAULT_SOLVER = _pick_solver()


def solve_problem(prob: cp.Problem, solver: str | None = None) -> str:
    """就地求解 cvxpy 问题，返回状态字符串（'optimal' / 'infeasible' / ...）。"""
    prob.solve(solver=solver or DEFAULT_SOLVER)
    return prob.status
