"""状态与无损离散化辅助（framework §5.6 state.py / paper_notes §4-5）。

系统动力学离散（Prop.1，区间内 c 常值）消去 c 后得到 a、b 的梯形关系：

    a[k] - a[k-1] = (b[k] + b[k-1]) * Δ[k]           （等式约束）
    c[k] = (b[k] - b[k-1]) / Δ[k]                    （区间控制，派生量）

时间代价的梯形权重（对应 ∫ ds/√a 的数值积分）：

    w[k] = 0.5 * (s[k+1] - s[k-1])   （内部点），端点为 0（a 已由边界固定）。

静止边界（a_bnd 端 a=0）的零进给奇异用 **Box I / Proposition 2**（论文式 20）在优化器
里以 jerk 零阶保持离散：头/尾各 N_s / N_f 个静止区间内 ⃛u≡κ 恒定，得

    a_k = ((u_k-u_s)/(u_1-u_s))^{4/3} a_1,  b_k = 2 a_k / (3(u_k-u_s)),  c_k = 2 a_k / (9(u_k-u_s)²)

即整段由**单一自由变量** a_1（尾段 a_{N-2}）线性决定（见 static_relations）。这样
求解出的 (a,b,c) 在静止段与 jerk-ZOH 自洽（关节 jerk 边界值非零 = q'·κ），与 interp/viz 一致。
"""

from __future__ import annotations

import numpy as np

from ..types import Topp3Data


def resolve_num_stationary(data: Topp3Data, n_stat: int) -> tuple[int, int]:
    """由边界是否为零进给 + 用户 N_s，定出实际 (N_s, N_f)。

    a_bnd 端 ≈0 才启用该端静止段；并留至少 1 个非静止中段、两端不重叠。
    """
    N = data.n_grid
    ns = int(n_stat) if abs(data.a_bnd[0]) <= 1e-12 else 0
    nf = int(n_stat) if abs(data.a_bnd[1]) <= 1e-12 else 0
    max_each = max(0, (N - 2) // 2)  # 保证中段 ≥1 个区间
    return max(0, min(ns, max_each)), max(0, min(nf, max_each))


def static_relations(a, b, s_grid: np.ndarray, num_stat: tuple[int, int]):
    """Box I 静止段离散（论文式 20）。返回 (cons, dyn_mask, c_point)。

    cons     : 静止段 a、b 的线性等式约束列表（cvxpy）
    dyn_mask : (N-1,) bool；True 的区间施加 c-ZOH 梯形动力学（静止区间为 False）
    c_point  : {网格点 k: c(u_k) 的 cvxpy 表达式}（式 20 点值，供 jerk 约束）
    """
    ns, nf = num_stat
    N = s_grid.size
    u = s_grid
    cons = []
    dyn_mask = np.ones(N - 1, dtype=bool)
    c_point: dict[int, object] = {}

    if ns > 0:  # 头段：u_s=u[0]，自由变量 a[1]
        us = u[0]
        d1 = u[1] - us
        for k in range(1, ns + 1):
            dk = u[k] - us
            cons.append(b[k] == (2.0 / (3.0 * dk)) * a[k])            # 式20 b_k
            if k >= 2:
                cons.append(a[k] == (dk / d1) ** (4.0 / 3.0) * a[1])  # 式20 a_k
            c_point[k] = (2.0 / (9.0 * dk ** 2)) * a[k]               # 式20 c_k（点值）
            dyn_mask[k - 1] = False                                    # 区间 k = diff 索引 k-1

    if nf > 0:  # 尾段：u_f=u[N-1]，自由变量 a[N-2]；b 取负（减速）
        uf = u[N - 1]
        dlast = uf - u[N - 2]
        for j in range(1, nf + 1):
            k = N - 1 - j
            dk = uf - u[k]
            cons.append(b[k] == -(2.0 / (3.0 * dk)) * a[k])
            if j >= 2:
                cons.append(a[k] == (dk / dlast) ** (4.0 / 3.0) * a[N - 2])
            c_point[k] = (2.0 / (9.0 * dk ** 2)) * a[k]
            dyn_mask[k] = False                                        # 区间 k+1 = diff 索引 k

    return cons, dyn_mask, c_point


def deltas_s(s_grid: np.ndarray) -> np.ndarray:
    """区间宽度 Δ[k]=s[k]-s[k-1]，长度 N-1（对应区间 1..N-1）。"""
    return np.diff(s_grid)


def trapz_weights(s_grid: np.ndarray) -> np.ndarray:
    """时间目标的梯形积分权重，形状 (N,)，端点为 0。"""
    N = s_grid.size
    w = np.zeros(N)
    w[1:-1] = 0.5 * (s_grid[2:] - s_grid[:-2])
    return w


def velocity_upper_bound(data: Topp3Data, flags=None) -> np.ndarray:
    """速度类约束 → a 的逐点上界 ā[k]（framework §6）。

    含（各由 flags 开关控制）：
      - 轴向速度 velocity：min_i vmax_i²/q'_i[k]²；
      - TCP 位置速度模 tcp_velocity、姿态角速度模 tcp_angular_velocity（若 data.tcp 存在）。
    均为 a 的线性上界，逐点取 min；被关闭的项不计入。未受约束处为 +inf。
    """
    from ..options import ConstraintFlags

    flags = flags or ConstraintFlags()
    N = data.n_grid
    a_bar = np.full(N, np.inf)

    if flags.velocity:
        # 启用 t–n 时轴速上界取空载转速 ω0（力矩→0 处，体现在 t–n 曲线里），
        # 比不考虑 t–n 的保守 vmax 更大；否则用 vmax。
        if flags.speed_torque and data.speed_torque is not None:
            vcap = data.speed_torque.noload_speed
        else:
            vcap = data.vmax
        dq2 = data.dq**2  # (n, N)
        with np.errstate(divide="ignore", invalid="ignore"):
            bound = (vcap[:, None] ** 2) / dq2  # (n, N)
        bound[~np.isfinite(bound)] = np.inf
        a_bar = np.minimum(a_bar, bound.min(axis=0))

    if data.tcp is not None and (flags.tcp_velocity or flags.tcp_angular_velocity):
        from ..constraints import tcp_a_upper_bound

        a_bar = np.minimum(a_bar, tcp_a_upper_bound(
            data.tcp, N, position=flags.tcp_velocity, orientation=flags.tcp_angular_velocity,
        ))
    return a_bar
