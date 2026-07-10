"""SPLP 结果可视化（framework §5.8 附带；可选依赖 matplotlib）。

  plot_splp_result        —— 2×3 概览（收敛/速度剖面/加速度/约束利用率/时间律/关节速度）
  plot_kinematic_limits   —— 关节 q̇/q̈/q⃛/力矩（带约束）+ TCP 位置速度模、姿态角速度模
  plot_fig4_interpolation —— 复现论文 Fig.4：可行轨迹的区间解析插值（Prop.1 + Prop.2）

以纯函数暴露，核心求解不依赖本模块。
"""

from __future__ import annotations

import numpy as np

from .types import Topp3Data, Profile
from .solve.state import velocity_upper_bound
from .solve.interp import t_to_s


def _set_cjk_font() -> None:
    """让 matplotlib 正常显示中文（Windows 常见 CJK 字体）。"""
    import matplotlib

    matplotlib.rcParams["font.sans-serif"] = [
        "Microsoft YaHei", "SimHei", "DengXian", "SimSun",
    ]
    matplotlib.rcParams["axes.unicode_minus"] = False


def reconstruct_grid_signals(data: Topp3Data, profile: Profile) -> dict:
    """在网格点重构物理量（时间导数）。返回 dict。"""
    a, b, c = profile.a, profile.b, profile.c
    s = data.s_grid
    N = data.n_grid
    sqrt_a = np.sqrt(np.maximum(a, 0.0))
    qd = data.dq * sqrt_a[None, :]                                    # q'·√a
    qdd = data.ddq * a[None, :] + data.dq * b[None, :]               # q''a+q'b
    # 参数 jerk ⃛u：非静止点 = c·√a（c-ZOH）；静止段 = 段几何常数 κ=(2/9)A^{3/2}
    #   （c 现为纯 c-ZOH 控制，√a·c 在静止段不等于 κ，故静止区间用几何单独算）
    ubar = c * sqrt_a
    ns, nf = getattr(profile, "num_stationary", (0, 0))
    if ns > 0:
        A_h = max(a[ns], 0.0) / (s[ns] - s[0]) ** (4.0 / 3.0)
        ubar[: ns + 1] = (2.0 / 9.0) * A_h ** 1.5          # 头段静止点（含 rest 点 0）
    if nf > 0:
        A_f = max(a[N - 1 - nf], 0.0) / (s[N - 1] - s[N - 1 - nf]) ** (4.0 / 3.0)
        ubar[N - 1 - nf:] = (2.0 / 9.0) * A_f ** 1.5       # 尾段静止点（含 rest 点 N-1）
    # q⃛ = q'''a^{3/2} + 3q''√a·b + q'·⃛u（用 (√a)³ 避免 a 端点微负产生 nan；rest 端 =q'·κ）
    qddd = (
        data.dddq * (sqrt_a ** 3)[None, :]
        + 3.0 * data.ddq * (sqrt_a * b)[None, :]
        + data.dq * ubar[None, :]
    )
    out = {"sdot": sqrt_a, "qd": qd, "qdd": qdd, "qddd": qddd}
    if data.torque is not None:  # M4：关节力矩 τ = n_tor·a + m_tor·b + g_tor
        out["torque"] = (
            data.torque.n_tor * a[None, :]
            + data.torque.m_tor * b[None, :]
            + data.torque.g_tor
        )
    return out


def reconstruct_time_signals(data: Topp3Data, profile: Profile) -> dict:
    """在**区间内细剖面**上重构时间域信号（q̇/q̈/q⃛/力矩），保证导数自洽。

    与 reconstruct_grid_signals（仅网格点）不同：这里用 interp.fine_profiles 给出的
    区间内闭式 a(σ)、b(σ)、⃛u(σ)（Prop.1/2）逐点重构，因此 q̈ 与 q⃛ 在时间上导数一致
    （网格点+线性插值会在静止段把 q̈∝σ^{1/3} 画成 ∝σ，使起点斜率与 jerk 对不上）。

    返回 dict：t（时间栅格，单调）、s、sdot、qd、qdd、qddd、(torque)。
    """
    from .solve.interp import fine_profiles

    fp = fine_profiles(data.s_grid, profile)
    s_f, t_f, a_f, b_f, ub_f = fp["s"], fp["t"], fp["a"], fp["b"], fp["ubar"]
    sqrt_a = np.sqrt(np.maximum(a_f, 0.0))

    def _to_fine(M):  # (n,N) 路径量 → (n, len(s_f))，沿 s 线性插值（q',q'',q''' 光滑）
        return np.vstack([np.interp(s_f, data.s_grid, M[i]) for i in range(M.shape[0])])

    qp, qpp, qppp = _to_fine(data.dq), _to_fine(data.ddq), _to_fine(data.dddq)
    qd = qp * sqrt_a
    qdd = qpp * a_f + qp * b_f
    # q⃛ = q'''a^{3/2} + 3q''√a·b + q'·⃛u（⃛u 直接用，静止段 =κ 有限）
    qddd = qppp * a_f ** 1.5 + 3.0 * qpp * sqrt_a * b_f + qp * ub_f

    out = {"t": t_f, "s": s_f, "sdot": sqrt_a, "qd": qd, "qdd": qdd, "qddd": qddd}
    if data.torque is not None:
        nt, mt, gt = (_to_fine(data.torque.n_tor), _to_fine(data.torque.m_tor),
                      _to_fine(data.torque.g_tor))
        out["torque"] = nt * a_f + mt * b_f + gt
    return out


def _interp_to_time(s_grid, s_of_t, signal_over_s):
    """把网格量（沿 s）插值到时间栅格（沿 s(t)）。"""
    return np.interp(s_of_t, s_grid, signal_over_s)


def plot_splp_result(
    data: Topp3Data,
    profile: Profile,
    hist=None,
    dt: float = 1e-3,
    save_path: str | None = None,
    show: bool = False,
    title: str = "TOTP-SPLP（M1）求解结果",
):
    """SPLP 结果概览图（2×3）。返回 matplotlib Figure。"""
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    _set_cjk_font()
    import matplotlib.pyplot as plt

    s = data.s_grid
    sig = reconstruct_grid_signals(data, profile)
    a_bar = velocity_upper_bound(data)
    t_uniform, s_of_t = t_to_s(s, profile, dt)

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(title, fontsize=14)

    # ① SPLP 收敛
    if hist is not None and getattr(hist, "t_final", None):
        it = np.arange(1, len(hist.t_final) + 1)
        ax[0, 0].plot(it, hist.t_final, "o-", color="C3")
        ax[0, 0].set_xlabel("SPLP 迭代次数 $p$")
        ax[0, 0].set_ylabel(r"终止时间 $t_f$ (s)")
        ax[0, 0].set_title("① SPLP 收敛")
        ax[0, 0].set_xticks(it)
    else:
        ax[0, 0].axis("off")

    # ② 速度剖面 ṡ(s) 与上界
    ax[0, 1].plot(s, sig["sdot"], color="C0", label=r"$\dot s=\sqrt{a}$")
    ax[0, 1].plot(s, np.sqrt(a_bar), "--", color="0.5", label=r"速度上界 $\sqrt{\bar a}$")
    ax[0, 1].set_xlabel("路径参数 $s$"); ax[0, 1].set_ylabel(r"$\dot s$")
    ax[0, 1].set_title("② 速度剖面"); ax[0, 1].legend(fontsize=9)

    # ③ 路径加速度 b(s)=s̈
    ax[0, 2].plot(s, profile.b, color="C1")
    ax[0, 2].axhline(0, color="0.7", lw=0.8)
    ax[0, 2].set_xlabel("路径参数 $s$"); ax[0, 2].set_ylabel(r"$b=\ddot s$")
    ax[0, 2].set_title("③ 路径加速度")

    # ④ 约束利用率 vs s
    r_v = np.max(np.abs(sig["qd"]) / data.vmax[:, None], axis=0)
    r_a = np.max(np.abs(sig["qdd"]) / data.amax[:, None], axis=0)
    r_j = np.max(np.abs(sig["qddd"]) / data.jmax[:, None], axis=0)
    ax[1, 0].plot(s, r_v, label="速度"); ax[1, 0].plot(s, r_a, label="加速度")
    ax[1, 0].plot(s, r_j, label="jerk")
    ax[1, 0].axhline(1.0, ls="--", color="0.5")
    ax[1, 0].set_xlabel("路径参数 $s$"); ax[1, 0].set_ylabel("利用率 (|·|/上限)")
    ax[1, 0].set_title("④ 约束利用率"); ax[1, 0].legend(fontsize=9)
    ax[1, 0].set_ylim(0, 1.15)

    # ⑤ 时间律 s(t)
    ax[1, 1].plot(t_uniform, s_of_t, color="C2")
    ax[1, 1].set_xlabel("时间 $t$ (s)"); ax[1, 1].set_ylabel("路径参数 $s$")
    ax[1, 1].set_title("⑤ 时间律 $s(t)$")

    # ⑥ 关节速度 q̇_i(t)（区间内细剖面重构，与加速度导数自洽；限位线按各关节自身 vmax_i，同色虚线）
    rec = reconstruct_time_signals(data, profile)
    for i in range(data.n_axis):
        color = f"C{i}"
        ax[1, 2].plot(rec["t"], rec["qd"][i], color=color, label=f"关节 {i}")
        ax[1, 2].axhline(data.vmax[i], color=color, ls=":", lw=0.8, alpha=0.6)
        ax[1, 2].axhline(-data.vmax[i], color=color, ls=":", lw=0.8, alpha=0.6)
    ax[1, 2].set_xlabel("时间 $t$ (s)"); ax[1, 2].set_ylabel(r"$\dot q_i$")
    ax[1, 2].set_title(r"⑥ 关节速度 $\dot q_i(t)$（虚线=各关节自身 $v_{max,i}$，同色）")
    ax[1, 2].legend(fontsize=9, ncol=2)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if save_path:
        fig.savefig(save_path, dpi=120)
    if show:
        plt.show()
    return fig


def plot_kinematic_limits(
    data: Topp3Data,
    profile: Profile,
    tcp: dict | None = None,
    dt: float = 1e-3,
    save_path: str | None = None,
    show: bool = False,
    title: str = "运动学/动力学约束（关节速度/加速度/jerk/力矩 与 TCP 位置速度、姿态角速度）随时间",
):
    """逐轴关节速度/加速度/jerk/力矩（带约束）+ 笛卡尔 TCP 位置速度模、姿态角速度模（带约束）。

    力矩面板仅在 data.torque 存在时绘制（M4）。

    tcp: 可选 dict（M1 无 FK，用合成解析 TCP 路径），字段：
        {"dp":  (3,N) 位置对 s 的导数 p'(s)（→ 位置速度模 = ‖p'‖·√a）,
         "wdir":(3,N) 单位路径速度下的角速度 ω/ṡ（→ 姿态角速度模 = ‖wdir‖·√a）,
         "v_max": 位置速度上界, "w_max": 姿态角速度上界}
    """
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    _set_cjk_font()
    import matplotlib.pyplot as plt

    s = data.s_grid
    rec = reconstruct_time_signals(data, profile)   # 区间内细剖面，导数自洽
    t = rec["t"]

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(title, fontsize=14)

    def _joint_panel(axis, sig_key, limit, name, sym):
        """limit：(n,) 各关节自身上界（对称）；限位线按关节自身值画，与曲线同色。"""
        lim = np.asarray(limit)
        for i in range(data.n_axis):
            color = f"C{i}"
            axis.plot(t, rec[sig_key][i], color=color, label=f"关节 {i}", lw=1.2)
            axis.axhline(lim[i], color=color, ls="--", lw=0.8, alpha=0.6)
            axis.axhline(-lim[i], color=color, ls="--", lw=0.8, alpha=0.6)
        axis.set_xlabel("时间 $t$ (s)"); axis.set_ylabel(sym)
        axis.set_title(f"{name}（虚线=各关节自身上限，同色）"); axis.legend(fontsize=8, ncol=2)
        lim_max = float(np.max(lim))
        axis.set_ylim(-1.25 * lim_max, 1.25 * lim_max)

    _joint_panel(ax[0, 0], "qd", data.vmax, r"① 关节速度 $\dot q_i$", r"$\dot q$")
    _joint_panel(ax[0, 1], "qdd", data.amax, r"② 关节加速度 $\ddot q_i$", r"$\ddot q$")
    _joint_panel(ax[0, 2], "qddd", data.jmax, r"③ 关节 jerk $\dddot q_i$", r"$\dddot q$")

    # ④ 关节力矩 τ_i(t)（M4；上下界可非对称、逐关节不同——UR5 各轴力矩上限差异很大，
    #    限位线必须按关节自身 tau_max_i/tau_min_i 画，不能用跨关节 max/min 的单一包络线）
    if data.torque is not None and "torque" in rec:
        tau_max = data.torque.tau_max
        tau_min = data.torque.tau_min
        for i in range(data.n_axis):
            color = f"C{i}"
            ax[1, 0].plot(t, rec["torque"][i], color=color, label=f"关节 {i}", lw=1.2)
            ax[1, 0].axhline(tau_max[i], color=color, ls="--", lw=0.8, alpha=0.6)
            ax[1, 0].axhline(tau_min[i], color=color, ls="--", lw=0.8, alpha=0.6)
        ax[1, 0].set_xlabel("时间 $t$ (s)"); ax[1, 0].set_ylabel(r"$\tau$")
        ax[1, 0].set_title(r"④ 关节力矩 $\tau_i$（虚线=各关节自身上下限，同色）")
        ax[1, 0].legend(fontsize=8, ncol=2)
        tmax, tmin = float(np.max(tau_max)), float(np.min(tau_min))
        pad = 0.25 * max(abs(tmax), abs(tmin), 1e-9)
        ax[1, 0].set_ylim(tmin - pad, tmax + pad)
    else:
        ax[1, 0].axis("off")

    # ⑤ TCP 位置速度模 / ⑥ TCP 姿态角速度模（同一细剖面 sdot）
    if tcp is not None:
        dp = np.asarray(tcp["dp"])            # (3,N)
        wdir = np.asarray(tcp["wdir"])        # (3,N)
        sqrt_a = rec["sdot"]
        cv_f = np.interp(rec["s"], s, np.linalg.norm(dp, axis=0))
        cw_f = np.interp(rec["s"], s, np.linalg.norm(wdir, axis=0))
        v_t = cv_f * sqrt_a
        w_t = cw_f * sqrt_a

        ax[1, 1].plot(t, v_t, color="C0", lw=1.4, label=r"$\|\dot p\|$")
        ax[1, 1].axhline(tcp["v_max"], ls="--", color="0.4", label="约束")
        ax[1, 1].set_xlabel("时间 $t$ (s)"); ax[1, 1].set_ylabel(r"$\|\dot p\|$")
        ax[1, 1].set_title(r"⑤ TCP 位置速度模 $\|\dot p\|$"); ax[1, 1].legend(fontsize=9)
        ax[1, 1].set_ylim(0, 1.25 * max(tcp["v_max"], float(v_t.max())))

        ax[1, 2].plot(t, w_t, color="C4", lw=1.4, label=r"$\|\omega\|$")
        ax[1, 2].axhline(tcp["w_max"], ls="--", color="0.4", label="约束")
        ax[1, 2].set_xlabel("时间 $t$ (s)"); ax[1, 2].set_ylabel(r"$\|\omega\|$")
        ax[1, 2].set_title(r"⑥ TCP 姿态角速度模 $\|\omega\|$"); ax[1, 2].legend(fontsize=9)
        ax[1, 2].set_ylim(0, 1.25 * max(tcp["w_max"], float(w_t.max())))
    else:
        ax[1, 1].axis("off"); ax[1, 2].axis("off")

    fig.text(
        0.5, 0.015,
        "注：TCP 曲线来自真实 UR5 DH 运动学（Jacobian 沿合成关节路径求值）；"
        "力矩系数仍是对角近似（下游集中质量单摆臂估算，非精确 RNE，见 robot/ur5.py）。"
        "各面板限位线均按关节自身实际配置值（同色虚线），非跨关节 max/min 包络。"
        "M4：TCP 速度模与关节力矩均为求解器约束。",
        ha="center", fontsize=9, color="0.35",
    )

    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    if save_path:
        fig.savefig(save_path, dpi=120)
    if show:
        plt.show()
    return fig


def fig4_interpolation_example(
    n_stat: int = 2,
    c_tail: tuple[float, ...] = (0.20, 0.05, -0.18),
    a_head: float = 0.6,
    du: float = 1.0,
    n_sub: int = 80,
):
    """构造论文 Fig.4 的示意轨迹（静止起点 a_s=b_s=0、非静止终点），返回逐区间连续量。

    头部 n_stat 个区间：**恒定参数 jerk** ⃛u≡κ（Prop.2，jerk 零阶保持），
        a(u)=A_h·(u-u_s)^{4/3}, b(u)=B_h·(u-u_s)^{1/3}, A_h=(3/2)B_h；
    其余区间：**c 零阶保持**（Prop.1），a(u)=a_{k-1}+2b_{k-1}σ+c_kσ²，b=b_{k-1}+c_kσ；
        参数 jerk ⃛u=c_k·√a 在区间内随 a 变、可在网格点间断。

    返回 dict：u_fine/a_fine/b_fine/j_fine（连续曲线，j=⃛u）、u_grid/a_grid/b_grid（网格点）、
    u_stat（静止段右端 = n_stat·du）、kappa。
    """
    Ls = n_stat * du
    A_h = a_head / Ls ** (4.0 / 3.0)
    B_h = (2.0 / 3.0) * A_h
    kappa = (B_h / 6.0 ** (1.0 / 3.0)) ** 1.5  # 头部恒定参数 jerk ⃛u_s

    u_grid = [0.0]
    a_grid = [0.0]
    b_grid = [0.0]
    u_parts, a_parts, b_parts, j_parts = [], [], [], []

    # 头部静止段（jerk-ZOH），整段一条 (u-u_s)^{4/3} 曲线
    u_h = np.linspace(0.0, Ls, n_stat * n_sub + 1)
    u_parts.append(u_h)
    a_parts.append(A_h * u_h ** (4.0 / 3.0))
    b_parts.append(B_h * u_h ** (1.0 / 3.0))
    j_parts.append(np.full_like(u_h, kappa))
    for k in range(1, n_stat + 1):
        u_grid.append(k * du)
        a_grid.append(A_h * (k * du) ** (4.0 / 3.0))
        b_grid.append(B_h * (k * du) ** (1.0 / 3.0))

    # 尾部非静止段（c-ZOH），逐区间前向积分
    a_prev, b_prev = a_grid[-1], b_grid[-1]
    for j, c in enumerate(c_tail):
        u_lo = (n_stat + j) * du
        sig = np.linspace(0.0, du, n_sub + 1)
        a_seg = a_prev + 2.0 * b_prev * sig + c * sig ** 2
        b_seg = b_prev + c * sig
        u_parts.append(u_lo + sig)
        a_parts.append(a_seg)
        b_parts.append(b_seg)
        j_parts.append(c * np.sqrt(np.maximum(a_seg, 0.0)))
        a_prev = a_prev + 2.0 * b_prev * du + c * du ** 2
        b_prev = b_prev + c * du
        u_grid.append(u_lo + du)
        a_grid.append(a_prev)
        b_grid.append(b_prev)

    return {
        "u_fine": np.concatenate(u_parts),
        "a_fine": np.concatenate(a_parts),
        "b_fine": np.concatenate(b_parts),
        "j_fine": np.concatenate(j_parts),
        "u_grid": np.asarray(u_grid),
        "a_grid": np.asarray(a_grid),
        "b_grid": np.asarray(b_grid),
        "u_stat": Ls,
        "kappa": kappa,
    }


def plot_fig4_interpolation(
    ex: dict | None = None,
    save_path: str | None = None,
    show: bool = False,
    title: str = "论文 Fig.4 复现：可行轨迹的区间解析插值（静止起点 + 非静止终点）",
):
    """复现论文 Fig.4：由离散网格点 (a_k,b_k,c_k) 解析重构区间内连续 a(u)、b(u) 与参数 jerk。

    ex 为 fig4_interpolation_example(...) 的输出；None 时用论文默认设定（N_s=2）。
    返回 matplotlib Figure。
    """
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    _set_cjk_font()
    import matplotlib.pyplot as plt

    if ex is None:
        ex = fig4_interpolation_example()

    u, us = ex["u_fine"], ex["u_stat"]
    ug = ex["u_grid"]

    fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    fig.suptitle(title, fontsize=13)

    panels = [
        (ax[0], ex["a_fine"], ex["a_grid"], "C0", r"① 路径速度平方 $a=\dot s^2$", r"$a(u)$"),
        (ax[1], ex["b_fine"], ex["b_grid"], "C1", r"② 路径加速度 $b=\ddot s$", r"$b(u)$"),
    ]
    for axis, yf, yg, color, ttl, ylab in panels:
        axis.axvspan(0.0, us, color="0.85", alpha=0.6)  # 静止段
        axis.plot(u, yf, color=color, lw=1.8)
        axis.scatter(ug, yg, color=color, zorder=5, s=28, edgecolor="k", linewidths=0.5)
        axis.set_ylabel(ylab)
        axis.set_title(ttl, fontsize=11)
        axis.grid(alpha=0.25)

    # ③ 参数 jerk ⃛u：头部恒定（jerk-ZOH），尾部 c·√a 随区间变、网格点可间断
    ax[2].axvspan(0.0, us, color="0.85", alpha=0.6)
    ax[2].plot(u, ex["j_fine"], color="C3", lw=1.8)
    for xg in ug:
        ax[2].axvline(xg, color="0.7", lw=0.6, ls=":")
    ax[2].axhline(0.0, color="0.6", lw=0.8)
    ax[2].set_ylabel(r"$\dddot u$")
    ax[2].set_title(r"③ 参数 jerk $\dddot u=c\sqrt{a}$（头部 $N_s$ 段恒定；尾部逐区间、网格点可间断）", fontsize=11)
    ax[2].set_xlabel(r"路径参数 $u$（网格点 $u_k$；灰色 = 静止段 $[u_s, u_{N_s}]$）")
    ax[2].grid(alpha=0.25)

    # 起末点标注
    ax[0].annotate(r"$a_s=b_s=0$", xy=(0, 0), xytext=(0.15, 0.12 * ex["a_grid"].max()),
                   fontsize=10, color="0.25")
    ax[0].annotate(r"$a_f\neq0$", xy=(ug[-1], ex["a_grid"][-1]),
                   xytext=(ug[-1] - 1.1, ex["a_grid"][-1]), fontsize=10, color="0.25")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if save_path:
        fig.savefig(save_path, dpi=120)
    if show:
        plt.show()
    return fig
