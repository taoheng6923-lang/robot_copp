"""单次 LP 组装与求解（framework §5.6 lp_problem.py / 论文式 29）。

决策变量 x=[a, b, J]。约束：
  - 边界：a[0],a[N-1],b[0],b[N-1]
  - 梯形动力学：a[k]-a[k-1]=(b[k]+b[k-1])Δ[k]（state）
  - 速度上界：a ≤ ā（state.velocity_upper_bound）
  - 轴向加速度：±(q''·a + q'·b) ≤ amax（逐轴逐点）
  - 轴向 jerk（在 a_lin 处线性化）：linearize.jerk_constraints
  - PLP 割线 + 下界 a≥δ0：plp_objective.build_plp
目标：min Σ w_k J_k（PLP 分段线性时间代价）。
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp

from ..types import Topp3Data, Profile
from ..backend import solve_problem
from .state import deltas_s, velocity_upper_bound, static_relations
from .linearize import jerk_constraints
from .plp_objective import PlpObjective, build_plp


def build_and_solve(
    data: Topp3Data, a_lin: np.ndarray, plp: PlpObjective,
    solver: str | None = None, num_stat: tuple[int, int] = (0, 0), flags=None,
) -> Profile:
    """组装一次 PLP-LP（论文式 29）并求解，返回 Profile(a,b,c)。

    num_stat 为 Box I 静止段；flags（ConstraintFlags）控制各约束启用。
    """
    from ..flags import ConstraintFlags

    flags = flags or ConstraintFlags()
    N, n = data.n_grid, data.n_axis
    ds = deltas_s(data.s_grid)
    a_bar = velocity_upper_bound(data, flags)

    a = cp.Variable(N, name="a")
    b = cp.Variable(N, name="b")
    jvar = cp.Variable(N, name="J")

    stat_cons, dyn_mask, c_point = static_relations(a, b, data.s_grid, num_stat)
    m = np.where(dyn_mask)[0]      # 仅非静止区间施加 c-ZOH 梯形动力学
    fin = np.where(np.isfinite(a_bar))[0]  # 速度类约束（关闭时为 +inf，跳过）
    cons = [
        a[0] == data.a_bnd[0], a[N - 1] == data.a_bnd[1],
        b[0] == data.b_bnd[0], b[N - 1] == data.b_bnd[1],
        a[m + 1] - a[m] == cp.multiply(b[m + 1] + b[m], ds[m]),  # 动力学（非静止段）
    ]
    if fin.size:
        cons.append(a[fin] <= a_bar[fin])                       # 速度上界
    cons += stat_cons  # Box I 静止段（式 20）
    # 轴向加速度
    if flags.acceleration:
        for i in range(n):
            acc = cp.multiply(data.ddq[i], a) + cp.multiply(data.dq[i], b)
            cons += [acc <= data.amax[i], acc >= -data.amax[i]]
    # 轴向 jerk（线性化；静止点用 c_point 点值，跳过 rest 端）
    if flags.jerk:
        cons += jerk_constraints(data, a_lin, a, b, c_point=c_point, num_stat=num_stat)
    # 关节力矩（M4，2 阶精确线性；速度类 TCP 约束已折进 a_bar）
    if flags.torque and data.torque is not None:
        from ..constraints import torque_constraints

        cons += torque_constraints(data.torque, a, b)
    # PLP 目标 + 下界
    plp_cons, objective = build_plp(plp, a, jvar)
    cons += plp_cons

    prob = cp.Problem(cp.Minimize(objective), cons)
    status = solve_problem(prob, solver=solver)
    if a.value is None:
        raise RuntimeError(f"PLP-LP 求解失败：status={status}")

    a_val = np.asarray(a.value).ravel()
    b_val = np.asarray(b.value).ravel()
    # c 一律为 c-ZOH 区间控制 c[k]=(b[k]-b[k-1])/Δ[k]（interp 据此重构非静止区间，
    # 保证区间末端 a 精确回到 a[k]、拼接处连续）。静止段的参数 jerk ⃛u=κ 不落在 c 上，
    # 由 num_stationary 让 interp/viz 用段几何单独算（见 fine_profiles / reconstruct_*）。
    c_val = np.zeros(N)
    c_val[1:] = (b_val[1:] - b_val[:-1]) / ds
    return Profile(a=a_val, b=b_val, c=c_val, num_stationary=num_stat)
