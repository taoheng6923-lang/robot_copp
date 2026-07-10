"""3 阶 jerk 凹约束的切线线性化（framework §5.6 linearize.py / 论文 eq.32）。

轴向 jerk（paper_notes §4 推导，与 Rust formulation.rs:238 一致）：

    q⃛_i = √a · (q'''_i·a + 3 q''_i·b + q'_i·c)

约束 |q⃛_i| ≤ jmax_i 写成 paper 的 3 阶标准形（f=jmax_i，h=0）：

    ±(q'''_i·a + 3 q''_i·b + q'_i·c) ≤ jmax_i · a^{-1/2}

在参考点 a_lin 处对凹项 a^{-1/2} 取切线（eq.32，非保守）：

    a^{-1/2} ≥ (3·a_lin - a) / (2·a_lin^{3/2})

得仿射约束（把 -a 项移到左侧，两个符号都在 a 系数上 +jmax/(2 a_lin^{3/2})）：

    ( ±q'''_i + jmax_i/(2 a_lin^{3/2}) )·a + ±3 q''_i·b + ±q'_i·c ≤ 3 jmax_i/(2√a_lin)

c 为区间控制：非静止区间 c[k]=(b[k]-b[k-1])/Δ[k]（c 零阶保持）；**静止段**用 Box I
点值 c(u_k)=2a_k/(9(u_k-u_s)²)（式 20，经 c_point 传入）。在每个网格点对**左右极限**
c(u_k^±) 施加（论文离散约束 7c+7d）。零进给 rest 端（a 固定为 0，jerk 上界=∞）跳过。
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp

from ..types import Topp3Data
from .state import deltas_s

_FLOOR = 1e-10  # a_lin 下限，防止 1/√a_lin 溢出（对应 Rust a_linearization_floor）


def jerk_constraints(
    data: Topp3Data,
    a_lin: np.ndarray,
    a: cp.Variable,
    b: cp.Variable,
    c_point: dict | None = None,
    num_stat: tuple[int, int] = (0, 0),
) -> list[cp.Constraint]:
    """返回在 a_lin 处线性化的 jerk 约束列表（cvxpy），点式逐网格点施加。"""
    ds = deltas_s(data.s_grid)  # (N-1,)
    N = data.n_grid
    c_point = c_point or {}
    ns, nf = num_stat
    al = np.maximum(a_lin, _FLOOR)
    inv = 1.0 / (2.0 * al ** 1.5)
    rhs = 3.0 / (2.0 * np.sqrt(al))
    cons: list[cp.Constraint] = []

    for p in range(N):
        # rest 边界点（a≡0，jerk 上界 a^{-1/2}=∞）：约束平凡，跳过
        if (ns > 0 and p == 0) or (nf > 0 and p == N - 1):
            continue
        # 该点左右极限的 c 表达式
        if p in c_point:
            climits = [c_point[p]]                       # 静止点：c(u_k^±) 连续 = 点值
        else:
            climits = []
            if p >= 1:
                climits.append((b[p] - b[p - 1]) / ds[p - 1])   # 左极限 c(u_p^-)：区间 p
            if p <= N - 2:
                climits.append((b[p + 1] - b[p]) / ds[p])       # 右极限 c(u_p^+)：区间 p+1
        for cexpr in climits:
            for i in range(data.n_axis):
                j = data.jmax[i]
                dq_i, ddq_i, dddq_i = data.dq[i, p], data.ddq[i, p], data.dddq[i, p]
                rhs_i = j * rhs[p]
                acoef = dddq_i + j * inv[p]
                cons.append(acoef * a[p] + 3.0 * ddq_i * b[p] + dq_i * cexpr <= rhs_i)
                acoef_m = -dddq_i + j * inv[p]
                cons.append(acoef_m * a[p] - 3.0 * ddq_i * b[p] - dq_i * cexpr <= rhs_i)
    return cons
