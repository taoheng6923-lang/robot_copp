"""轨迹合成：(PathDerivatives, Profile) → 等时间栅格关节轨迹（framework §5.8 / 设计 §9）。

在 copp 的**区间内解析细剖面**（`copp.solve.interp.fine_profiles`，Prop.1 c-ZOH +
Prop.2 静止段 jerk-ZOH 闭式）上重构，与 `copp.viz.reconstruct_time_signals`
同一数学，保证 q̇/q̈/q⃛ 导数自洽（网格点+线性插值会在静止段把 q̈∝σ^{1/3} 画成
∝σ）。重构公式（设计 §9，c=b'=s⃛/ṡ 约定）：

    q̇ = q'·√a,   q̈ = q''·a + q'·b,   q⃛ = q'''·a^{3/2} + 3q''·√a·b + q'·⃛u

其中 ⃛u = c·√a 由细剖面直接给出（静止段 c→∞ 但 ⃛u≡κ 有限）。

插值说明：t↔s 与 (a,b,⃛u) 来自解析细剖面（每区间 24 细分），在其上线性插值到
等时间栅格的误差为细剖面网格的 O(Δ²)，远小于约束容差；q,q',q'',q''' 沿 s_grid
线性插值（几何量光滑、网格已按弦高容差自适应加密）。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from copp.solve.interp import fine_profiles
from copp.types import Profile

from path.types import PathDerivatives


@dataclass
class TrajectoryResult:
    """等时间栅格的关节轨迹（设计 §9 输出）。

    t                : (T,)   时间栅格（含精确终点 t_final；段拼接后严格递增
                              但 seam 处间隔可能小于 dt）
    s                : (T,)   对应路径参数（多段拼接时为**各段局部** s）
    sdot             : (T,)   ṡ=√a
    q/qd/qdd/qddd    : (n,T)  关节位置/速度/加速度/jerk
    v_tcp/w_tcp      : (T,)   TCP 位置速度模 ‖ṗ‖、姿态角速度模 ‖ω‖
    t_final          : float  总时长
    seg_index        : (T,)   int，每个采样点所属的指令段号（单段全 0）
    """

    t: np.ndarray
    s: np.ndarray
    sdot: np.ndarray
    q: np.ndarray
    qd: np.ndarray
    qdd: np.ndarray
    qddd: np.ndarray
    v_tcp: np.ndarray
    w_tcp: np.ndarray
    t_final: float
    seg_index: np.ndarray

    @property
    def n_axis(self) -> int:
        return self.q.shape[0]


def synthesize(pd: PathDerivatives, profile: Profile, dt: float = 1e-3) -> TrajectoryResult:
    """单段合成：解析细剖面 → 等时间栅格（末尾补精确 t_final 点，保证收在 rest）。"""
    fp = fine_profiles(pd.s_grid, profile)
    t_f, s_f = fp["t"], fp["s"]
    t_final = float(t_f[-1])

    t_u = np.arange(0.0, t_final, dt)
    if t_u.size == 0 or t_final - t_u[-1] > 1e-12:
        t_u = np.concatenate([t_u, [t_final]])

    # 时间域插值：s(t) 与 (a,b,⃛u)(t) 都在解析细剖面上取
    s_u = np.interp(t_u, t_f, s_f)
    a_u = np.maximum(np.interp(t_u, t_f, fp["a"]), 0.0)
    b_u = np.interp(t_u, t_f, fp["b"])
    ub_u = np.interp(t_u, t_f, fp["ubar"])
    sdot = np.sqrt(a_u)

    def _interp_rows(M: np.ndarray) -> np.ndarray:
        return np.vstack([np.interp(s_u, pd.s_grid, M[i]) for i in range(M.shape[0])])

    q_u = _interp_rows(pd.q)
    qp = _interp_rows(pd.dq)
    qpp = _interp_rows(pd.ddq)
    qppp = _interp_rows(pd.dddq)

    qd = qp * sdot
    qdd = qpp * a_u + qp * b_u
    qddd = qppp * sdot ** 3 + 3.0 * qpp * sdot * b_u + qp * ub_u

    cv_u = np.interp(s_u, pd.s_grid, pd.cv)
    cw_u = np.interp(s_u, pd.s_grid, pd.cw)

    return TrajectoryResult(
        t=t_u, s=s_u, sdot=sdot,
        q=q_u, qd=qd, qdd=qdd, qddd=qddd,
        v_tcp=cv_u * sdot, w_tcp=cw_u * sdot,
        t_final=t_final,
        seg_index=np.zeros(t_u.size, dtype=int),
    )


def concatenate(results: list[TrajectoryResult]) -> TrajectoryResult:
    """多段轨迹按时间顺序拼接（M2 语义：段间 rest 停顿衔接）。

    后续段丢弃各自 t=0 首样本（与上一段终点同一 rest 状态、同一关节角——
    G0 衔接已在 lower_sections 校验）。seam 处时间间隔 = 上一段末尾残步 + dt，
    不严格等距；s 列保持各段局部值，段归属看 seg_index。
    """
    if len(results) == 1:
        return results[0]
    t_off = 0.0
    parts: dict[str, list[np.ndarray]] = {k: [] for k in
        ("t", "s", "sdot", "q", "qd", "qdd", "qddd", "v_tcp", "w_tcp", "seg_index")}
    for i, r in enumerate(results):
        sl = slice(1, None) if i > 0 else slice(None)
        parts["t"].append(r.t[sl] + t_off)
        for key in ("s", "sdot", "v_tcp", "w_tcp"):
            parts[key].append(getattr(r, key)[sl])
        for key in ("q", "qd", "qdd", "qddd"):
            parts[key].append(getattr(r, key)[:, sl])
        parts["seg_index"].append(np.full(r.t[sl].size, i, dtype=int))
        t_off += r.t_final

    return TrajectoryResult(
        t=np.concatenate(parts["t"]),
        s=np.concatenate(parts["s"]),
        sdot=np.concatenate(parts["sdot"]),
        q=np.hstack(parts["q"]),
        qd=np.hstack(parts["qd"]),
        qdd=np.hstack(parts["qdd"]),
        qddd=np.hstack(parts["qddd"]),
        v_tcp=np.concatenate(parts["v_tcp"]),
        w_tcp=np.concatenate(parts["w_tcp"]),
        t_final=t_off,
        seg_index=np.concatenate(parts["seg_index"]),
    )
