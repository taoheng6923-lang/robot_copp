"""约束摄入辅助（framework §5.5 / 设计 §6，M4）。

- TCP 速度模长 → a 的逐点上界（线性，折进 ā）；
- 关节力矩 → LP 不等式行（2 阶，对 (a,b) 精确线性）。

TCP 与力矩系数在 M1 由合成模型给出；实际管线中 TCP 系数来自 Jacobian、
力矩系数来自逆动力学（M2 的 KinematicsModel / DynamicsModel）。
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp

from ..types import TcpConstraint, TorqueConstraint, SpeedTorqueConstraint


def tcp_a_upper_bound(
    tcp: TcpConstraint, n_grid: int,
    position: bool = True, orientation: bool = True,
) -> np.ndarray:
    """TCP 速度模长 → a 的逐点上界 (N,)：min(v_max²/cv², w_max²/cw²)。

    position/orientation 分别控制是否计入位置速度模、姿态角速度模（约束开关）。
    系数为 0 处（该方向不产生约束）用 +inf 占位。
    """
    bound = np.full(n_grid, np.inf)
    with np.errstate(divide="ignore", invalid="ignore"):
        if position:
            b_pos = np.where(tcp.cv > 0, tcp.v_max**2 / tcp.cv**2, np.inf)
            bound = np.minimum(bound, b_pos)
        if orientation:
            b_ori = np.where(tcp.cw > 0, tcp.w_max**2 / tcp.cw**2, np.inf)
            bound = np.minimum(bound, b_ori)
    return bound


def torque_constraints(
    torque: TorqueConstraint, a: cp.Variable, b: cp.Variable
) -> list[cp.Constraint]:
    """关节力矩约束的 cvxpy 行：τ_min ≤ n_tor·a + m_tor·b + g_tor ≤ τ_max（逐轴）。"""
    cons: list[cp.Constraint] = []
    for i in range(torque.n_tor.shape[0]):
        tau = (
            cp.multiply(torque.n_tor[i], a)
            + cp.multiply(torque.m_tor[i], b)
            + torque.g_tor[i]
        )
        cons += [tau <= torque.tau_max[i], tau >= torque.tau_min[i]]
    return cons


def speed_torque_envelope(st: SpeedTorqueConstraint, i: int, qd_abs: np.ndarray) -> np.ndarray:
    """第 i 轴的梯形可用力矩包络 τ_env(|q̇|)：平台 τ0 到 ω_c，线性收窄到 0 于 ω0。"""
    tau0, wc, w0 = st.tau0[i], st.rated_speed[i], st.noload_speed[i]
    roll = tau0 * np.clip((w0 - qd_abs) / (w0 - wc), 0.0, 1.0)   # ω_c..ω0 线性到 0
    return np.where(qd_abs <= wc, tau0, roll)


_ST_FLOOR = 1e-10  # a_lin 下限，防止 1/√a_lin 溢出（同 linearize._FLOOR）


def speed_torque_constraints(
    st: SpeedTorqueConstraint, dq: np.ndarray, a: cp.Variable, b: cp.Variable,
    a_lin: np.ndarray, num_stat: tuple[int, int] = (0, 0),
) -> list[cp.Constraint]:
    """速度相关力矩（t–n 梯形包络）约束 → 一组仿射-于-(a,b) 的 cvxpy 行（论文式 17→18 + SCP 切线）。

    真实约束 |τ_dyn + Fv·q̇ + Fc·sgn(q̇)| ≤ τ_env(|q̇|)（q̇=q'·√a）含凹项 √a → 非凸。分两步凸化：

    1) **梯形 = 两个 halfplane 的交**（精确，无近似）：
         τ_env(|q̇|) = min( τ0,  E_roll − s_roll·|q̇| ),   E_roll=τ0+s_roll·ω_c, s_roll=τ0/(ω0−ω_c)
       故 |τ_motor| ≤ τ_env 等价于**同时**满足平台约束(s=0,E=τ0)与 rolloff 约束(s=s_roll,E=E_roll)：
         驱动:  τ_dyn + (s·|q'| + Fv·q')·√a ≤ E − Fc·sgn
         制动: −τ_dyn + (s·|q'| − Fv·q')·√a ≤ E + Fc·sgn

    2) **√a 在 SCP 参考点 a_lin 处线性化**（与 jerk 同源，切点随迭代收敛到工作点，
       故静止端无"固定切点截距"伪速度项）：逐点按 √a 系数正负替换（都收紧=保守）——
         系数≥0 用**切线上界** √a ≤ a/(2√a_lin)+√a_lin/2；
         系数<0 用**过原点割线下界** √a ≥ a·|q'|/ω0（在 |q̇|≤ω0 内成立）。
       每条均仿射于 (a,b)。静止段（num_stat 头/尾）低速、力矩由重力主导（|g_tor|≤τ0 已保证），
       跳过以免小 a 处线性化过保守（同 jerk 跳过 rest 端）。见 docs/tn_constraint_notes.md §5、§8。
    """
    n, N = st.n_tor.shape
    ns, nf = num_stat
    al = np.asarray(a_lin, float)
    al = np.where(np.isfinite(al), al, np.nanmax(al[np.isfinite(al)]) if np.any(np.isfinite(al)) else 1.0)
    al = np.maximum(al, _ST_FLOOR)
    sa = np.sqrt(al)                                 # √a_lin，(N,)
    inv2sa = 1.0 / (2.0 * sa)                        # 切线斜率，(N,)
    active = np.ones(N, dtype=bool)                  # 跳过静止段（低速、重力主导）
    if ns > 0:
        active[:ns] = False
    if nf > 0:
        active[N - nf:] = False
    idx = np.where(active)[0]
    if idx.size == 0:
        return []

    cons: list[cp.Constraint] = []
    for i in range(n):
        tau_dyn = cp.multiply(st.n_tor[i], a) + cp.multiply(st.m_tor[i], b) + st.g_tor[i]
        pq = dq[i]                                   # q'_i(s)，(N,)
        aq = np.abs(pq)
        sgn = np.sign(pq)
        Fv, Fc = st.viscous[i], st.coulomb[i]
        wc, w0, tau0 = st.rated_speed[i], st.noload_speed[i], st.tau0[i]
        s_roll = tau0 / (w0 - wc)                     # rolloff 斜率（对 |q̇|）
        tan = cp.multiply(a, inv2sa) + sa / 2.0       # √a 切线上界（a_lin 处相切），(N,)
        sec = cp.multiply(a, aq / w0)                 # √a 过原点割线下界 √a≥a·|q'|/ω0，(N,)
        for s_j, E_j in ((0.0, tau0), (s_roll, tau0 + s_roll * wc)):   # 平台 + rolloff 两个 halfplane
            coef_d = s_j * aq + Fv * pq               # √a 系数（可正可负），(N,)
            coef_b = s_j * aq - Fv * pq
            bd = _select(coef_d, tan, sec)            # 逐点选 √a 收紧替换（≥0→切线上界，<0→割线下界）
            bb = _select(coef_b, tan, sec)
            drive = tau_dyn + cp.multiply(coef_d, bd)
            brake = -tau_dyn + cp.multiply(coef_b, bb)
            rhs_d = E_j - Fc * sgn                     # (N,)
            rhs_b = E_j + Fc * sgn
            cons += [drive[idx] <= rhs_d[idx], brake[idx] <= rhs_b[idx]]
    return cons


def _select(coef: np.ndarray, upper, lower):
    """逐点：coef≥0 取 upper（切线上界），coef<0 取 lower（割线下界）；两者均为 (N,) 仿射表达式。

    upper=pu·a+qu、lower=pl·a（均对 a 仿射）；返回 (mask·pu+~mask·pl)·a + mask·qu 的仿射组合。
    """
    # upper、lower 是 cvxpy 仿射表达式，用掩码线性组合逐元素挑选
    mask = (coef >= 0).astype(float)                  # (N,)
    return cp.multiply(mask, upper) + cp.multiply(1.0 - mask, lower)


def speed_torque_utilization(
    st: SpeedTorqueConstraint, dq: np.ndarray, a: np.ndarray, b: np.ndarray
) -> np.ndarray:
    """真实（未凸化）t–n 约束的利用率 |τ_motor| / τ_env，(n, N)；≤1 即满足。

    τ_motor=τ_dyn+Fv·q̇+Fc·sgn(q̇)，τ_env 为梯形可用力矩包络（q̇=q'·√a）。
    """
    n = st.n_tor.shape[0]
    sa = np.sqrt(np.maximum(a, 0.0))          # ṡ=√a，(N,)
    ratios = np.zeros((n, dq.shape[1]))
    for i in range(n):
        qd = dq[i] * sa                        # q̇_i，(N,)
        tau_dyn = st.n_tor[i] * a + st.m_tor[i] * b + st.g_tor[i]
        tau_motor = tau_dyn + st.viscous[i] * qd + st.coulomb[i] * np.sign(dq[i])
        tau_env = speed_torque_envelope(st, i, np.abs(qd))
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios[i] = np.abs(tau_motor) / np.maximum(tau_env, 1e-9)
    return ratios
