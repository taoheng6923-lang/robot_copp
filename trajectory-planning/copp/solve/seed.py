"""种子 a⁽⁰⁾（framework §5.6 seed.py / 设计 §7.2④、论文 §5.1）。

用 2 阶问题（忽略 3 阶 jerk 约束）求 a⁽⁰⁾：max Σ w_k a_k，s.t. 梯形动力学、
速度上界、加速度约束、边界。只需保证首次线性化 LP 可行（Theorem 2）。
这里以 LP 直接求解（对应 topp2_ra 的 2 阶最优解角色）。
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp

from ..types import Topp3Data
from ..backend import solve_problem
from .state import deltas_s, trapz_weights, velocity_upper_bound, static_relations

_SEED_FLOOR = 1e-6


def compute_seed(
    data: Topp3Data, num_stat: tuple[int, int] = (0, 0), flags=None
) -> np.ndarray:
    """返回 2 阶最优速度剖面 a⁽⁰⁾，形状 (N,)。num_stat 指定头/尾静止段（Box I）。

    flags（ConstraintFlags）控制各约束启用；种子与主 LP 用同一组开关保证可行性一致。
    """
    from ..options import ConstraintFlags

    flags = flags or ConstraintFlags()
    N, n = data.n_grid, data.n_axis
    ds = deltas_s(data.s_grid)
    a_bar = velocity_upper_bound(data, flags)
    w = trapz_weights(data.s_grid)

    a = cp.Variable(N, name="a_seed")
    b = cp.Variable(N, name="b_seed")

    stat_cons, dyn_mask, _ = static_relations(a, b, data.s_grid, num_stat)
    m = np.where(dyn_mask)[0]  # 仅非静止区间施加 c-ZOH 梯形动力学
    fin = np.where(np.isfinite(a_bar))[0]  # 速度类约束（关闭时为 +inf，跳过）
    cons = [
        a[0] == data.a_bnd[0], a[N - 1] == data.a_bnd[1],
        b[0] == data.b_bnd[0], b[N - 1] == data.b_bnd[1],
        a[m + 1] - a[m] == cp.multiply(b[m + 1] + b[m], ds[m]),   # 梯形动力学（非静止段）
        a[1:-1] >= _SEED_FLOOR,
    ]
    if fin.size:
        cons.append(a[fin] <= a_bar[fin])
    cons += stat_cons  # 静止段 Box I 关系（式 20）
    if flags.acceleration:
        for i in range(n):
            acc = cp.multiply(data.ddq[i], a) + cp.multiply(data.dq[i], b)
            cons += [acc <= data.amax[i], acc >= -data.amax[i]]
    # 力矩约束也纳入种子（M4），保证首次线性化 LP 的可行性一致（Theorem 2）
    if flags.torque and data.torque is not None:
        from ..constraints import torque_constraints

        cons += torque_constraints(data.torque, a, b)
    # 速度相关力矩（t–n 梯形包络）：√a 在速度上界 a_bar 处线性化（种子 max∫a 工作点≈a_bar，
    # 切线在此处紧、且随后 SPLP 迭代在各自 a_lin 处重线性化），保持可行性一致（Theorem 2）
    if flags.speed_torque and data.speed_torque is not None:
        from ..constraints import speed_torque_constraints

        cons += speed_torque_constraints(data.speed_torque, data.dq, a, b, a_bar, num_stat=num_stat)

    prob = cp.Problem(cp.Maximize(w @ a), cons)
    status = solve_problem(prob)
    if a.value is None:
        raise RuntimeError(f"种子 2 阶 LP 求解失败：status={status}")
    return np.maximum(np.asarray(a.value).ravel(), _SEED_FLOOR)
