"""解析插值 s↔t（framework §5.6 interp.py / 论文 Prop.1 + Prop.2）。

逐区间/段在细网格上给出闭式剖面 (a, b, ⃛u)，据此累积时间表 (s_fine, t_fine)；
s_to_t / t_to_s 在其上插值，viz 也用同一细剖面重构 q̇/q̈/q⃛（保证导数自洽）。

**非静止区间**（两端 a>0，Proposition 1，c 零阶保持）：

    a(σ) = a[k-1] + 2 b[k-1] σ + c[k] σ²,  b(σ) = b[k-1] + c[k] σ,  ⃛u(σ) = c[k]·√a(σ)
    t 增量 = ∫₀^{Δ} dσ/√a(σ) —— **Prop.1 解析闭式** Φ_k（_quad_time，按 c>0/c<0/c=0 分
             对数/反正弦/根式，论文 eq.11），无梯形离散误差。

**静止段**（rest 端 a=0，Proposition 2，⃛u 零阶保持；由 profile.num_stationary 指定）：

    头段：a(σ)=A_h σ^{4/3}, b(σ)=(2/3)A_h σ^{1/3}, ⃛u≡κ=(2/9)A_h^{3/2}   （σ=u-u_s）
    尾段：a(σ)=A_f ρ^{4/3}, b(σ)=-(2/3)A_f ρ^{1/3}, ⃛u≡(2/9)A_f^{3/2}     （ρ=u_f-u）

其中 ⃛u=c·√a 在静止段虽 c→∞ 但乘积有限恒定（=κ），故用 ⃛u 直接重构 q⃛，避免 c 的奇异。
段总时长 3(u_{N_s}-u_s)/√a_{N_s}（有限）；采样在 a→0 端加密。
"""

from __future__ import annotations

import numpy as np

from ..types import Profile

_A_FLOOR = 1e-12
_SUB = 24              # 每区间细分点数
A_STATIC = 1e-6        # a ≤ 此阈值视为静止端（(0,0) 回退时的逐区间判别）


_C_EPS = 1e-9   # |c| 判零阈值
_B_EPS = 1e-12  # |b| 判零阈值


def _quad_time(a0, b0, c, sig):
    """Prop.1 闭式时长 Φ_k：t(σ)=∫₀^σ dτ/√(a0+2b0τ+cτ²)，按 c 符号分情形（论文 eq.11）。

    c>0 → 对数（asinh 类）；c<0 → 反正弦；c=0 → 根式（b≠0）或线性（b=0）。
    """
    Q = np.maximum(a0 + 2.0 * b0 * sig + c * sig ** 2, _A_FLOOR)
    sqQ = np.sqrt(Q)
    sa0 = np.sqrt(max(a0, _A_FLOOR))

    if abs(c) < _C_EPS:                      # c ≈ 0
        if abs(b0) < _B_EPS:                 #   b ≈ 0：a≈常数
            return sig / sa0
        return (sqQ - sa0) / b0              #   ∫dτ/√(a0+2b0τ)=√Q/b0

    if c > 0.0:                              # c > 0：(1/√c)·ln|√c√Q+cσ+b0|
        rc = np.sqrt(c)
        F = np.log(np.abs(rc * sqQ + c * sig + b0) + _A_FLOOR)
        F0 = np.log(np.abs(rc * sa0 + b0) + _A_FLOOR)
        return (F - F0) / rc

    cc = -c                                  # c < 0：(1/√(-c))·arcsin((σ-b0/(-c))/R)
    rcc = np.sqrt(cc)
    R = np.sqrt(max(b0 ** 2 + cc * a0, _A_FLOOR)) / cc   # 配方半径 R=√(b0²+cc·a0)/cc
    arg = np.clip((sig - b0 / cc) / R, -1.0, 1.0)
    arg0 = np.clip((0.0 - b0 / cc) / R, -1.0, 1.0)
    return (np.arcsin(arg) - np.arcsin(arg0)) / rcc


# ── 区间/段的细剖面：返回 (σ 或 u, a, b, ⃛u, τ 累积时间) ─────────────────
def _interval_regular(a_lo, b_lo, c_hi, d):
    """非静止区间（Prop.1 c-ZOH）。返回 (sig, a, b, ubar, tau)，sig∈[0,d]。

    时长 tau 用 Prop.1 **解析闭式** Φ_k（_quad_time），非梯形近似。
    """
    sig = np.linspace(0.0, d, _SUB + 1)
    a_sig = np.maximum(a_lo + 2.0 * b_lo * sig + c_hi * sig ** 2, _A_FLOOR)
    b_sig = b_lo + c_hi * sig
    ub_sig = c_hi * np.sqrt(a_sig)
    tau = _quad_time(a_lo, b_lo, c_hi, sig)
    return sig, a_sig, b_sig, ub_sig, tau


def _interval_head_static(a_hi, d):
    """头部单区间 jerk-ZOH（(0,0) 回退用）。返回 (sig, a, b, ubar, tau)。"""
    A_h = max(a_hi, _A_FLOOR) / d ** (4.0 / 3.0)
    frac = np.linspace(0.0, 1.0, _SUB + 1) ** 3
    sig = d * frac
    a_sig = np.maximum(A_h * sig ** (4.0 / 3.0), _A_FLOOR)
    b_sig = (2.0 / 3.0) * A_h * sig ** (1.0 / 3.0)
    ub_sig = np.full_like(sig, (2.0 / 9.0) * A_h ** 1.5)
    tau = 3.0 * sig ** (1.0 / 3.0) * d ** (2.0 / 3.0) / np.sqrt(max(a_hi, _A_FLOOR))
    return sig, a_sig, b_sig, ub_sig, tau


def _interval_tail_static(a_lo, d):
    """尾部单区间 jerk-ZOH（(0,0) 回退用）。返回 (sig, a, b, ubar, tau)。"""
    A_f = max(a_lo, _A_FLOOR) / d ** (4.0 / 3.0)
    rho_frac = np.linspace(1.0, 0.0, _SUB + 1) ** 3
    sig = d * (1.0 - rho_frac)
    rho = d * rho_frac
    a_sig = np.maximum(A_f * rho ** (4.0 / 3.0), _A_FLOOR)
    b_sig = -(2.0 / 3.0) * A_f * rho ** (1.0 / 3.0)
    ub_sig = np.full_like(sig, (2.0 / 9.0) * A_f ** 1.5)
    tau = (3.0 * d / np.sqrt(max(a_lo, _A_FLOOR))) * (1.0 - rho_frac ** (1.0 / 3.0))
    return sig, a_sig, b_sig, ub_sig, tau


def _segment_head_static(u_grid, a_end, ns):
    """头部 N_s 个静止区间作一段 jerk-ZOH。返回 (u, a, b, ubar, tau)。"""
    u0 = u_grid[0]
    L = u_grid[ns] - u0
    A_h = max(a_end, _A_FLOOR) / L ** (4.0 / 3.0)
    frac = np.linspace(0.0, 1.0, ns * _SUB + 1) ** 3
    u_seg = u0 + L * frac
    sg = u_seg - u0
    a_sig = np.maximum(A_h * sg ** (4.0 / 3.0), _A_FLOOR)
    b_sig = (2.0 / 3.0) * A_h * sg ** (1.0 / 3.0)
    ub_sig = np.full_like(u_seg, (2.0 / 9.0) * A_h ** 1.5)
    tau = 3.0 * sg ** (1.0 / 3.0) / np.sqrt(A_h)
    return u_seg, a_sig, b_sig, ub_sig, tau


def _segment_tail_static(u_grid, a_start, nf):
    """尾部 N_f 个静止区间作一段 jerk-ZOH。返回 (u, a, b, ubar, tau)。

    细分点按 rho_frac³ 向 rest 端聚点（与 _interval_tail_static/头段一致）：
    ṡ∝ρ^{2/3} 下 t∝L^{1/3}−ρ^{1/3}，ρ 取立方分布使时间近似均匀——否则
    s-均匀采样的最后一个细分在时间上占整段 (1/24)^{1/3}≈35%，时域插值失真。
    """
    N = u_grid.size
    ustart, uf = u_grid[N - 1 - nf], u_grid[N - 1]
    L = uf - ustart
    A_f = max(a_start, _A_FLOOR) / L ** (4.0 / 3.0)
    rho_frac = np.linspace(1.0, 0.0, nf * _SUB + 1) ** 3
    u_seg = ustart + L * (1.0 - rho_frac)
    rho = L * rho_frac
    a_sig = np.maximum(A_f * rho ** (4.0 / 3.0), _A_FLOOR)
    b_sig = -(2.0 / 3.0) * A_f * rho ** (1.0 / 3.0)
    ub_sig = np.full_like(u_seg, (2.0 / 9.0) * A_f ** 1.5)
    tau = 3.0 / np.sqrt(A_f) * (L ** (1.0 / 3.0) - rho ** (1.0 / 3.0))
    return u_seg, a_sig, b_sig, ub_sig, tau


def fine_profiles(s_grid: np.ndarray, profile: Profile) -> dict:
    """逐区间/段拼接的细剖面。返回 dict：s, t, a, b, ubar（⃛u）。

    t 为累积到达时间；静止头/尾按 Box I jerk-ZOH，其余按 c-ZOH。供 s_to_t/t_to_s
    与 viz 的时间域信号重构共用（保证 q̇/q̈/q⃛ 在同一剖面上导数自洽）。
    """
    a, b, c = profile.a, profile.b, profile.c
    N = s_grid.size
    ns, nf = getattr(profile, "num_stationary", (0, 0))
    s_p = [np.array([s_grid[0]])]
    t_p = [np.array([0.0])]
    a_p = [np.array([max(a[0], 0.0)])]
    b_p = [np.array([b[0]])]
    u_p = [np.array([c[1] * np.sqrt(max(a[0], 0.0)) if ns == 0 else 0.0])]
    t_acc = 0.0

    def _push(s_seg, a_seg, b_seg, ub_seg, tau):
        nonlocal t_acc
        s_p.append(s_seg[1:]); a_p.append(a_seg[1:]); b_p.append(b_seg[1:])
        u_p.append(ub_seg[1:]); t_p.append(t_acc + tau[1:])
        t_acc += tau[-1]

    if ns > 0:  # 头段（整段 jerk-ZOH）
        seg = _segment_head_static(s_grid, a[ns], ns)
        u_p[0] = np.array([seg[3][0]])  # 修正起点 ⃛u=κ
        _push(seg[0], seg[1], seg[2], seg[3], seg[4])

    lo = ns + 1 if ns > 0 else 1
    hi = N - 1 - nf if nf > 0 else N - 1
    for k in range(lo, hi + 1):
        d = s_grid[k] - s_grid[k - 1]
        a_lo, a_hi = max(a[k - 1], 0.0), max(a[k], 0.0)
        if ns == 0 and a_lo <= A_STATIC:
            sig, asg, bsg, ubsg, tau = _interval_head_static(a_hi, d)
        elif nf == 0 and a_hi <= A_STATIC:
            sig, asg, bsg, ubsg, tau = _interval_tail_static(a_lo, d)
        else:
            sig, asg, bsg, ubsg, tau = _interval_regular(a[k - 1], b[k - 1], c[k], d)
        _push(s_grid[k - 1] + sig, asg, bsg, ubsg, tau)

    if nf > 0:  # 尾段（整段 jerk-ZOH）
        seg = _segment_tail_static(s_grid, a[N - 1 - nf], nf)
        _push(seg[0], seg[1], seg[2], seg[3], seg[4])

    return {
        "s": np.concatenate(s_p), "t": np.concatenate(t_p),
        "a": np.concatenate(a_p), "b": np.concatenate(b_p),
        "ubar": np.concatenate(u_p),
    }


def s_to_t(s_grid: np.ndarray, profile: Profile) -> tuple[float, np.ndarray]:
    """返回 (t_final, t_s)，t_s[k] 为到达 s[k] 的时间（t_s[0]=0）。"""
    fp = fine_profiles(s_grid, profile)
    t_s = np.interp(s_grid, fp["s"], fp["t"])
    return float(fp["t"][-1]), t_s


def t_to_s(
    s_grid: np.ndarray, profile: Profile, dt: float
) -> tuple[np.ndarray, np.ndarray]:
    """等时间栅格采样：返回 (t_uniform, s_of_t)。"""
    fp = fine_profiles(s_grid, profile)
    t_uniform = np.arange(0.0, fp["t"][-1], dt)
    s_of_t = np.interp(t_uniform, fp["t"], fp["s"])
    return t_uniform, s_of_t
