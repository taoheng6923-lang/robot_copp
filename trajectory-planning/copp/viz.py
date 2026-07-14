"""SPLP 结果可视化（framework §5.8 附带；可选依赖 matplotlib）。

  plot_splp_result        —— 2×3 概览（收敛/速度剖面/加速度/约束利用率/时间律/关节速度）
  plot_kinematic_limits   —— 关节 q̇/q̈/q⃛/力矩（带约束）+ TCP 位置速度模、姿态角速度模
  plot_speed_torque       —— 速度相关力矩（t–n）分析：转矩–转速包络/利用率/摩擦分量
  plot_tn_convexification —— 论文 Fig.3：(q̇²,τ) 真实非凸域 vs 仿射切角内逼近
  plot_interp_profiles    —— 论文 Fig.4 风格：由本用例 profile 全程解析重构 a/b/c/参数jerk

上列全部以 (data, profile) 等**唯一求解结果**为输入，不另造示意数据。

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


def constraint_utilizations(data: Topp3Data, profile: Profile, flags=None) -> dict:
    """逐约束的利用率曲线（|·|/上限，沿 s 网格），**只含 flags 启用者**。

    返回 {标签: (N,)}，覆盖全部七类约束；=1 即该约束在此处绑定。被开关关闭的约束
    不出现在结果里（它们不参与求解，画出来只会误导——看着"超限"其实是没施加）。

    注意速度一项的分母：启用 t–n 时轴速上界是空载转速 ω0 而非 vmax（与
    solve/state.velocity_upper_bound 一致），否则会把合法的高速误算成超限。
    """
    from .options import ConstraintFlags

    flags = flags or ConstraintFlags()
    a, b = profile.a, profile.b
    sa = np.sqrt(np.maximum(a, 0.0))
    sig = reconstruct_grid_signals(data, profile)
    out: dict[str, np.ndarray] = {}

    if flags.velocity:
        vcap = (data.speed_torque.noload_speed
                if (flags.speed_torque and data.speed_torque is not None) else data.vmax)
        out["速度"] = np.max(np.abs(sig["qd"]) / vcap[:, None], axis=0)
    if flags.acceleration:
        out["加速度"] = np.max(np.abs(sig["qdd"]) / data.amax[:, None], axis=0)
    if flags.jerk:
        out["jerk"] = np.max(np.abs(sig["qddd"]) / data.jmax[:, None], axis=0)
    if flags.torque and data.torque is not None:
        tq = data.torque
        tau = tq.n_tor * a + tq.m_tor * b + tq.g_tor
        with np.errstate(divide="ignore", invalid="ignore"):
            # 双边盒式：τ>0 归一到 tau_max、τ<0 归一到 tau_min（<0，故商为正），取较紧者
            ratio = np.maximum(tau / tq.tau_max[:, None], tau / tq.tau_min[:, None])
        out["力矩（盒式）"] = np.max(np.nan_to_num(ratio, nan=0.0, posinf=0.0), axis=0)
    if flags.speed_torque and data.speed_torque is not None:
        out["力矩（t–n）"] = speed_torque_signals(data, profile)["util"].max(axis=0)
    if flags.tcp_velocity and data.tcp is not None:
        out["TCP 线速度"] = data.tcp.cv * sa / data.tcp.v_max
    if flags.tcp_angular_velocity and data.tcp is not None:
        out["TCP 角速度"] = data.tcp.cw * sa / data.tcp.w_max
    return out


def plot_splp_result(
    data: Topp3Data,
    profile: Profile,
    hist=None,
    dt: float = 1e-3,
    flags=None,
    save_path: str | None = None,
    show: bool = False,
    title: str = "TOTP-SPLP（M1）求解结果",
):
    """SPLP 结果概览图（2×3）。返回 matplotlib Figure。

    flags : ConstraintFlags；决定面板 ④ 画哪些约束的利用率、以及面板 ② 的速度上界
            按哪些约束合成。None 则视为全部启用。
    """
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    _set_cjk_font()
    import matplotlib.pyplot as plt

    from .options import ConstraintFlags

    flags = flags or ConstraintFlags()
    s = data.s_grid
    sig = reconstruct_grid_signals(data, profile)
    a_bar = velocity_upper_bound(data, flags)   # 上界也只由启用的速度类约束合成
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

    # ④ 约束利用率 vs s —— 只画 flags 启用的约束（关闭者不参与求解，画出来会误导）
    utils = constraint_utilizations(data, profile, flags)
    if utils:
        for k, (name, r) in enumerate(utils.items()):
            ax[1, 0].plot(s, r, color=f"C{k}", label=name)
        ax[1, 0].axhline(1.0, ls="--", color="0.5")
        top = max(1.15, 1.05 * max(float(np.max(r)) for r in utils.values()))
        ax[1, 0].set_ylim(0, top)   # 自适应，避免贴边曲线被裁掉看不见
        ax[1, 0].legend(fontsize=8, ncol=2)
    else:
        ax[1, 0].text(0.5, 0.5, "无启用的约束", ha="center", va="center",
                      transform=ax[1, 0].transAxes, color="0.5")
    ax[1, 0].set_xlabel("路径参数 $s$"); ax[1, 0].set_ylabel("利用率 (|·|/上限)")
    ax[1, 0].set_title(f"④ 约束利用率（仅启用者：{len(utils)}/7）")

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


def speed_torque_signals(data: Topp3Data, profile: Profile) -> dict:
    """在网格点重构 t–n 约束相关量（供 plot_speed_torque / 断言核对）。

    返回逐轴 (n, N)：qd（q̇）、tau_dyn、tau_motor、tau_avail、util（|τ_motor|/τ_avail）、
    以及摩擦分量 visc（Fv·q̇）、coul（Fc·sgn(q̇)）；标量 t（网格到达时间）。
    """
    from .solve.interp import s_to_t
    from .constraints import speed_torque_envelope

    st = data.speed_torque
    a, b = profile.a, profile.b
    sa = np.sqrt(np.maximum(a, 0.0))
    qd = data.dq * sa[None, :]                                   # (n,N) q̇=q'·√a
    tau_dyn = st.n_tor * a[None, :] + st.m_tor * b[None, :] + st.g_tor
    visc = st.viscous[:, None] * qd
    coul = st.coulomb[:, None] * np.sign(data.dq)
    tau_motor = tau_dyn + visc + coul
    tau_avail = np.stack([speed_torque_envelope(st, i, np.abs(qd[i]))   # 梯形包络 τ_env(|q̇|)
                          for i in range(st.n_tor.shape[0])])
    with np.errstate(divide="ignore", invalid="ignore"):
        util = np.abs(tau_motor) / np.maximum(tau_avail, 1e-9)
    _, t_grid = s_to_t(data.s_grid, profile)
    return {"t": t_grid, "qd": qd, "tau_dyn": tau_dyn, "tau_motor": tau_motor,
            "tau_avail": tau_avail, "util": util, "visc": visc, "coul": coul}


def plot_speed_torque(
    data: Topp3Data,
    profile: Profile,
    joints: tuple[int, ...] = (0, 1, 2),
    save_path: str | None = None,
    show: bool = False,
    title: str = "速度相关力矩（t–n）约束分析：转矩–转速包络、约束利用率与摩擦分量",
):
    """t–n 约束分析图（论文 Fig.6 风格 + 利用率/摩擦分解）。需 data.speed_torque。

    2×3 面板：
      ①②③ 选定关节的 τ_motor–q̇ 散点 + 可用力矩包络 ±(τ0−κ|q̇|)（点色=利用率）；
      ④ 各关节约束利用率 |τ_motor|/τ_avail 随时间（=1 即绑定）；
      ⑤ 最紧关节的 |τ_motor| 与可用力矩 τ_avail 随时间（可用力矩随速度收窄）；
      ⑥ 最紧关节的力矩构成：动力学 τ_dyn + 粘滞 Fv·q̇ + 库仑 Fc·sgn → τ_motor。
    """
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    _set_cjk_font()
    import matplotlib.pyplot as plt

    st = data.speed_torque
    sig = speed_torque_signals(data, profile)
    t = sig["t"]
    util = sig["util"]
    jstar = int(np.argmax(util.max(axis=1)))          # 最紧（利用率最高）关节

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(title, fontsize=14)

    # ①②③ τ_motor–q̇ 散点 + 可用力矩包络（Fig.6 风格）
    for col, i in enumerate(joints[:3]):
        axis = ax[0, col]
        qd_i = sig["qd"][i]
        qmax = max(np.abs(qd_i).max(), 1e-6) * 1.08
        qg = np.linspace(-qmax, qmax, 200)
        from .constraints import speed_torque_envelope
        env = speed_torque_envelope(st, i, np.abs(qg))           # 梯形包络 τ_env(|q̇|)
        axis.plot(qg, env, "r-", lw=1.3, label=r"$\pm\tau_{env}(|\dot q|)$（梯形）")
        axis.plot(qg, -env, "r-", lw=1.3)
        sc = axis.scatter(qd_i, sig["tau_motor"][i], c=util[i], cmap="viridis",
                          vmin=0.0, vmax=1.0, s=14, zorder=3)
        axis.set_xlabel(r"$\dot q_i$ (rad/s)"); axis.set_ylabel(r"$\tau$ (N·m)")
        axis.set_title(f"① 关节 {i}：$\\tau_{{motor}}$ vs $\\dot q$（点色=利用率）"
                       if col == 0 else f"关节 {i}")
        axis.grid(True, alpha=0.3)
        if col == 0:
            axis.legend(fontsize=8, loc="best")
    fig.colorbar(sc, ax=ax[0, 2], label="利用率 |τ_motor|/τ_avail")

    # ④ 各关节利用率随时间
    for i in range(data.n_axis):
        ax[1, 0].plot(t, util[i], color=f"C{i}", lw=1.1, label=f"关节 {i}")
    ax[1, 0].axhline(1.0, color="0.3", ls="--", lw=1.0, label="绑定线 =1")
    ax[1, 0].set_xlabel("时间 $t$ (s)"); ax[1, 0].set_ylabel("利用率")
    ax[1, 0].set_title("④ 各关节 t–n 约束利用率随时间")
    ax[1, 0].legend(fontsize=8, ncol=2); ax[1, 0].set_ylim(0, 1.15)

    # ⑤ 最紧关节：|τ_motor| 与可用力矩 τ_avail
    ax[1, 1].plot(t, np.abs(sig["tau_motor"][jstar]), color="C3", lw=1.4, label=r"$|\tau_{motor}|$")
    ax[1, 1].plot(t, sig["tau_avail"][jstar], color="0.3", ls="--", lw=1.4,
                  label=r"可用力矩 $\tau_{env}(|\dot q|)$（梯形）")
    ax[1, 1].set_xlabel("时间 $t$ (s)"); ax[1, 1].set_ylabel(r"$\tau$ (N·m)")
    ax[1, 1].set_title(f"⑤ 最紧关节 {jstar}：需求 vs 可用力矩（可用随速度收窄）")
    ax[1, 1].legend(fontsize=9); ax[1, 1].grid(True, alpha=0.3)

    # ⑥ 最紧关节的力矩构成（动力学 + 粘滞 + 库仑）
    ax[1, 2].plot(t, sig["tau_dyn"][jstar], color="C0", lw=1.3, label=r"动力学 $\tau_{dyn}$")
    ax[1, 2].plot(t, sig["visc"][jstar], color="C1", lw=1.3, label=r"粘滞 $F_v\dot q$")
    ax[1, 2].plot(t, sig["coul"][jstar], color="C2", lw=1.3, label=r"库仑 $F_c\,\mathrm{sgn}$")
    ax[1, 2].plot(t, sig["tau_motor"][jstar], color="C3", lw=1.6, label=r"合计 $\tau_{motor}$")
    ax[1, 2].axhline(0.0, color="0.7", lw=0.6)
    ax[1, 2].set_xlabel("时间 $t$ (s)"); ax[1, 2].set_ylabel(r"$\tau$ (N·m)")
    ax[1, 2].set_title(f"⑥ 最紧关节 {jstar} 力矩构成（摩擦并入约束边界）")
    ax[1, 2].legend(fontsize=8); ax[1, 2].grid(True, alpha=0.3)

    fig.text(
        0.5, 0.015,
        r"t–n 约束 $|\tau_{dyn}+F_v\dot q+F_c\,\mathrm{sgn}| \leq \tau_{env}(|\dot q|)$"
        "（梯形可用力矩：低速平台+高速线性收窄；叠加粘滞+库仑摩擦）；求解按 SPLP 在 a_lin 处"
        "对 √a 切线线性化（保守内逼近、收敛处精确）。参数为合成 stand-in，见 docs/tn_constraint_notes.md。",
        ha="center", fontsize=9, color="0.35",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    if save_path:
        fig.savefig(save_path, dpi=120)
    if show:
        plt.show()
    return fig


def plot_tn_convexification(
    data: Topp3Data,
    profile: Profile,
    save_path: str | None = None,
    show: bool = False,
    title: str = "t–n 约束凸化（论文 Fig.3）：$\\dot q^2$ 平面的真实非凸可行域 vs 仿射切角内逼近",
):
    """逐关节复现论文 **Fig.3**：在 $(\\dot q^2,\\tau)$ 平面直观显示"切割意图与效果"。

    **切割意图**：梯形可用力矩 τ_env(|q̇|) 的 rolloff 斜段在 $\\dot q^2=q'^2 a$ 坐标下是
    **凸曲线**（因 τ_env 仿射于 |q̇|=√(q̇²)），其**下方**可行域**非凸**（论文 Fig.3 虚线）。
    求解按论文思路（在 a_lin 处对 √a 切线线性化）用**仿射直线**内逼近，切掉曲线与切线间的角。

    **切割效果**：蓝色凸内逼近 ⊆ 灰色真实域（保守/安全，绝不超真实边界）；工作点（点色=
    利用率）落在凸内逼近内。本图 rolloff 段以**平台角点单切线**示意最保守的一次切割；实际
    求解逐点在各自 a_lin 处相切、凸域更贴近曲线（收敛处精确）。参数为合成 stand-in。
    """
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    _set_cjk_font()
    import matplotlib.pyplot as plt
    from .constraints import speed_torque_envelope

    st = data.speed_torque
    sig = speed_torque_signals(data, profile)
    n = st.n_tor.shape[0]
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4.2 * nrows), squeeze=False)
    fig.suptitle(title, fontsize=14)

    sc = None
    for i in range(n):
        ax = axes[i // ncols, i % ncols]
        wc, w0, tau0 = float(st.rated_speed[i]), float(st.noload_speed[i]), float(st.tau0[i])
        s_roll = tau0 / (w0 - wc)
        xc = wc ** 2                                       # q̇² 拐点（平台末端）

        xg = np.linspace(0.0, w0 ** 2, 400)                # q̇² 全程 [0, ω0²]
        env = speed_torque_envelope(st, i, np.sqrt(xg))    # τ_env(√x)：梯形→q̇²（rolloff 弯曲）

        # ① 真实可行域 |τ|≤τ_env：灰填 + 虚线边界（rolloff 段为凸曲线→下方非凸）
        ax.fill_between(xg, -env, env, color="0.86", zorder=0)
        ax.plot(xg, env, color="0.35", ls="--", lw=1.4, label=r"真实边界 $\pm\tau_{env}$（非凸）")
        ax.plot(xg, -env, color="0.35", ls="--", lw=1.4)

        # ② 仿射切角内逼近：平台角 (ω_c², τ0) 处对 rolloff 凸曲线引切线（切线在曲线下=保守）
        slope = -s_roll / (2.0 * wc)                       # dτ_env/dx |_{x=ω_c²}
        tan = tau0 + slope * (xg - xc)
        approx = np.clip(np.minimum(tau0, tan), 0.0, None)  # 凸上界=min(平台 τ0, 切线)，≥0 不反转
        ax.plot(xg, approx, color="C0", lw=1.7, label=r"凸内逼近（切线，保守）")
        ax.plot(xg, -approx, color="C0", lw=1.7)
        ax.fill_between(xg, -approx, approx, facecolor="C0", alpha=0.10, zorder=1)

        # ③ 被切掉的非凸角：真实域与凸逼近之差（rolloff 段），红斜纹
        cut = env - approx > 1e-9
        ax.fill_between(xg, approx, env, where=cut, facecolor="none",
                        edgecolor="C3", hatch="xxx", lw=0.0, zorder=2,
                        label="切掉的非凸角")
        ax.fill_between(xg, -env, -approx, where=cut, facecolor="none",
                        edgecolor="C3", hatch="xxx", lw=0.0, zorder=2)

        # ④ 平台|rolloff 分界 + 实际工作点 (q̇², τ_motor)，点色=利用率
        ax.axvline(xc, color="0.55", ls=":", lw=1.0)
        ax.text(xc, tau0 * 1.04, r"$\omega_c^2$", fontsize=8, color="0.4", ha="center")
        x_op = sig["qd"][i] ** 2
        sc = ax.scatter(x_op, sig["tau_motor"][i], c=sig["util"][i], cmap="viridis",
                        vmin=0.0, vmax=1.0, s=18, zorder=3, edgecolor="k", linewidth=0.25,
                        label="工作点（点色=利用率）")

        ax.set_xlabel(r"$\dot q_i^2$ (rad²/s²)")
        ax.set_ylabel(r"$\tau_i$ (N·m)")
        ax.set_title(f"关节 {i}：$\\tau_0$={tau0:.0f}, $\\omega_c$={wc:.2f}, $\\omega_0$={w0:.2f}")
        ax.set_ylim(-tau0 * 1.18, tau0 * 1.18)
        ax.grid(True, alpha=0.25)
        if i == 0:
            ax.legend(fontsize=7, loc="lower left")

    for j in range(n, nrows * ncols):        # 关掉多余子图
        axes[j // ncols, j % ncols].axis("off")
    fig.text(0.5, 0.008,
             r"论文 Fig.3：梯形可行域映射到 $\dot q^2=q'^2 a$ 后 rolloff 边界弯曲→下方非凸；"
             r"用仿射切线内逼近切掉非凸角（此处示意单条角点切线，最保守），保证凸且保守。"
             r"实际求解 SPLP 逐点在各自 a_lin 处相切、收敛处贴合真实曲线；工作点全落在真实域内。",
             ha="center", fontsize=9, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 0.9, 0.96))
    if sc is not None:                        # 独立色条轴，避免挤占子图
        cax = fig.add_axes((0.915, 0.12, 0.014, 0.76))
        fig.colorbar(sc, cax=cax, label="利用率 |τ_motor|/τ_env")
    if save_path:
        fig.savefig(save_path, dpi=120)
    if show:
        plt.show()
    return fig


def plot_interp_profiles(
    data: Topp3Data,
    profile: Profile,
    save_path: str | None = None,
    show: bool = False,
    title: str = "区间解析插值（论文 Fig.4 风格）：本用例轨迹全程 a/b/c/参数jerk 重构",
):
    """用**本测试用例**求得的 (data, profile) 解析重构全程 a(s)/b(s)/c(s)/⃛u(s)（论文 Fig.4 复现）。

    数据源自唯一用例（solve/interp.fine_profiles 逐区间闭式插值，**覆盖全部插补周期**）：静止头/尾
    按 Box I / Prop.2 jerk-ZOH、其余按 Prop.1 c-ZOH。控制量 c=⃛u/√a 在非静止段=c_k（逐区间常值、
    网格点可间断），静止端 a→0 时 c→∞（正是静止段用 jerk-ZOH 而非 c-ZOH 的原因）。返回 Figure。
    """
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    _set_cjk_font()
    import matplotlib.pyplot as plt
    from .solve.interp import fine_profiles

    sg = data.s_grid
    N = sg.size
    fp = fine_profiles(sg, profile)
    s, a_f, b_f, j_f = fp["s"], fp["a"], fp["b"], fp["ubar"]          # ⃛u = ubar
    with np.errstate(divide="ignore", invalid="ignore"):
        c_f = np.where(a_f > 1e-9, j_f / np.sqrt(np.maximum(a_f, 1e-12)), np.nan)  # c=⃛u/√a
    ns, nf = getattr(profile, "num_stationary", (0, 0))
    stat_spans = []
    if ns > 0:
        stat_spans.append((sg[0], sg[ns]))
    if nf > 0:
        stat_spans.append((sg[N - 1 - nf], sg[N - 1]))

    fig, ax = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
    fig.suptitle(title, fontsize=13)

    def _shade(axis):
        for lo, hi in stat_spans:
            axis.axvspan(lo, hi, color="0.85", alpha=0.6)            # 静止段

    # ①② 状态 a、b：连续曲线 + 网格点
    for axis, yf, yg, color, ttl, ylab in (
        (ax[0], a_f, profile.a, "C0", r"① 路径速度平方 $a=\dot s^2$", r"$a(s)$"),
        (ax[1], b_f, profile.b, "C1", r"② 路径加速度 $b=\ddot s$", r"$b(s)$"),
    ):
        _shade(axis)
        axis.plot(s, yf, color=color, lw=1.5)
        axis.scatter(sg, yg, color=color, zorder=5, s=14, edgecolor="k", linewidths=0.35)
        axis.set_ylabel(ylab)
        axis.set_title(ttl, fontsize=11)
        axis.grid(alpha=0.25)

    # ③ 控制量 c=b'=⃛u/√a：非静止段逐区间常值（c-ZOH，网格点可间断）；静止端 →∞
    _shade(ax[2])
    ax[2].plot(s, c_f, color="C2", lw=1.3)
    ax[2].axhline(0.0, color="0.6", lw=0.8)
    nonstat = np.ones_like(s, dtype=bool)                            # ylim 只按**非静止段**的 c（=c_k，有界）
    for lo, hi in stat_spans:                                       # 静止端 c→∞ 任其冲出上/下沿
        nonstat &= ~((s >= lo) & (s <= hi))
    cc = c_f[nonstat & np.isfinite(c_f)]
    if cc.size:
        lo, hi = float(cc.min()), float(cc.max())
        pad = 0.1 * (hi - lo) + 1e-6
        ax[2].set_ylim(lo - pad, hi + pad)
    ax[2].set_ylabel(r"$c=b'$")
    ax[2].set_title(r"③ 控制量 $c=b'=\dddot u/\sqrt{a}$（非静止段 c-ZOH 逐区间常值；静止端 →∞ 故用 jerk-ZOH）", fontsize=11)
    ax[2].grid(alpha=0.25)

    # ④ 参数 jerk ⃛u：静止段恒定=κ；非静止段 c·√a 随区间变、可间断
    _shade(ax[3])
    ax[3].plot(s, j_f, color="C3", lw=1.3)
    ax[3].axhline(0.0, color="0.6", lw=0.8)
    ax[3].set_ylabel(r"$\dddot u$")
    ax[3].set_title(r"④ 参数 jerk $\dddot u=c\sqrt{a}$（静止段恒定=κ；非静止段逐区间、可间断）", fontsize=11)
    ax[3].set_xlabel(r"路径参数 $s$（网格点 $s_k$；灰色 = 静止段，两端 rest-to-rest）")
    ax[3].grid(alpha=0.25)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    if save_path:
        fig.savefig(save_path, dpi=120)
    if show:
        plt.show()
    return fig


