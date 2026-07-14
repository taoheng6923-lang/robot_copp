"""copp 数值内核端到端测试（framework §8）。

**唯一测试用例** test_copp：一次求解得 (data, profile)，全部断言 + 由**同一求解结果**输出五张
分析图——output/ 下所有图像均来自该唯一用例，不另造示意数据：
  output/splp_test.png               —— SPLP 概览（收敛 / 速度剖面 / 约束利用率 / 时间律…）
  output/splp_limits_test.png        —— 关节速度/加速度/jerk/力矩 + TCP 速度模（带约束）
  output/speed_torque_test.png       —— 速度相关力矩（t–n）：转矩–转速包络 / 利用率 / 摩擦分量
  output/tn_convexification_test.png —— 论文 Fig.3：逐关节 (q̇²,τ) 真实非凸域 vs 仿射切角内逼近
  output/fig4_interpolation.png      —— 论文 Fig.4 风格：本用例轨迹全程 a/b/c/参数jerk 解析重构

可用 pytest 运行，也可直接 `python trajectory-planning/copp/self-test/test_copp.py`。
"""

from __future__ import annotations

import os
import sys

_TRAJ_PLANNING_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_REPO_ROOT = os.path.dirname(_TRAJ_PLANNING_DIR)
sys.path.insert(0, _TRAJ_PLANNING_DIR)   # 供 `import copp` 使用
sys.path.insert(0, _REPO_ROOT)           # 供 `import robot`（顶层，不在 trajectory-planning 下）使用

import numpy as np

from copp import (
    Topp3Data, ConstraintFlags, solve_splp, SolveOptions,
    load_robot_limits, load_constraint_flags, load_smooth_c_weight,
)
from copp.solve import s_to_t
from robot import UR5RobotModel

# 约束配置唯一来源见 limits_config.py（改约束只改那里）
from limits_config import LIMITS as _LIMITS

# 机器人运动学/动力学计算集中于 UR5RobotModel（真实 UR5 DH 运动学 + 对角近似动力学）
_MODEL = UR5RobotModel(seed=3)


def _make_data(n_grid=81) -> Topp3Data:
    """合成数据（M4：TCP 速度模 + 关节力矩），本体量全部来自 UR5RobotModel。"""
    s = np.linspace(0.0, 1.0, n_grid)
    q0, q1, q2, q3 = _MODEL.joint_path(s)
    return _LIMITS.to_topp3_data(
        s, q1, q2, q3,
        tcp_geom=_MODEL.tcp_coeffs(s),
        torque_coeffs=_MODEL.torque_coeffs(q0, q1, q2),
    )


def _viz_tcp(s):
    """viz.plot_kinematic_limits 用的 TCP dict：几何 {dp,wdir} + 给定上界 v_max/w_max。"""
    return {**_MODEL.tcp_geometry(s), "v_max": _LIMITS.v_tcp_max, "w_max": _LIMITS.w_tcp_max}


def _assert_solver(data, profile, hist, flags=ConstraintFlags()):
    """求解结果的全部性质断言（仅校验 flags 启用的约束；单调收敛 / rest-to-rest / s↔t）。"""
    a, b, c = profile.a, profile.b, profile.c
    N = data.n_grid
    sa = np.sqrt(np.maximum(a, 0.0))

    # 形状
    assert a.shape == (N,) and b.shape == (N,) and c.shape == (N,)

    # 1) t_final 单调不增并收敛
    tf = hist.t_final
    assert all(x >= y - 1e-6 for x, y in zip(tf, tf[1:])), tf
    assert hist.converged

    # 2) 轴向 vel/acc/jerk 约束满足（仅校验启用者；网格点，含 2% 数值裕度）
    #    jerk 用 reconstruct_grid_signals 的正确 ⃛u（静止段用段几何 κ，非 c-ZOH 的 √a·c）
    from copp.viz import reconstruct_grid_signals
    sig = reconstruct_grid_signals(data, profile)
    qd, qdd, qddd = sig["qd"], sig["qdd"], sig["qddd"]
    binders = []
    if flags.velocity:
        # 轴速上界：启用 t–n 时取空载转速 ω0（=noload_speed），否则取 vmax——
        # 与求解端 solve/state.velocity_upper_bound 一致。此前死用 vmax，仅因 jerk
        # 约束恰好把速度压在 vmax 内而未暴露；关掉 jerk 后速度顶到 ω0 即误报越界。
        vcap = (data.speed_torque.noload_speed
                if (flags.speed_torque and data.speed_torque is not None) else data.vmax)
        r_v = np.max(np.abs(qd) / vcap[:, None]); assert r_v <= 1.02, r_v
    if flags.acceleration:
        r_a = np.max(np.abs(qdd) / data.amax[:, None]); assert r_a <= 1.02, r_a; binders.append(r_a)
    if flags.jerk:
        r_j = np.max(np.abs(qddd) / data.jmax[:, None]); assert r_j <= 1.02, r_j; binders.append(r_j)

    # 3) rest-to-rest（a_bnd=(0,0)）：边界 a≈0、首末关节速度严格为 0（Prop.2 静止段）
    assert data.a_bnd == (0.0, 0.0)
    assert abs(a[0]) < 1e-6 and abs(a[-1]) < 1e-6
    assert np.max(np.abs(qd[:, 0])) < 1e-8 and np.max(np.abs(qd[:, -1])) < 1e-8

    # 4) s↔t：时长有限、到达时间严格递增（未修复前零进给奇异会发散到数百秒）
    t_final, t_s = s_to_t(data.s_grid, profile)
    assert np.isfinite(t_final) and 0.0 < t_final < 10.0, t_final
    assert np.all(np.diff(t_s) > 0)

    # 5) M4：TCP 速度模、关节力矩满足约束（仅校验启用者）
    if flags.tcp_velocity:
        r = np.max(data.tcp.cv * sa) / data.tcp.v_max; assert r <= 1.02, r; binders.append(r)
    if flags.tcp_angular_velocity:
        r = np.max(data.tcp.cw * sa) / data.tcp.w_max; assert r <= 1.02, r; binders.append(r)
    if flags.torque:
        tau = data.torque.n_tor * a + data.torque.m_tor * b + data.torque.g_tor
        assert np.all(tau <= data.torque.tau_max[:, None] + 1e-3)
        assert np.all(tau >= data.torque.tau_min[:, None] - 1e-3)
        binders.append(float(np.max(np.abs(tau) / data.torque.tau_max[:, None])))
    # 速度相关力矩（t–n）：回代**真实（未凸化）**约束，利用率≤1、可用力矩恒正、且绑定
    if flags.speed_torque:
        from copp.viz import speed_torque_signals
        sig_st = speed_torque_signals(data, profile)
        assert np.all(sig_st["tau_avail"] > 0), "可用力矩 τ_env(|q̇|) 应恒为正"
        assert sig_st["util"].max() <= 1.0 + 2e-3, sig_st["util"].max()
        binders.append(float(sig_st["util"].max()))

    # 至少一个启用的约束真正绑定（否则用例无意义）
    assert binders and max(binders) >= 0.5, binders


def _assert_config(flags):
    """机器人配置加载：逐关节 (n,) 约束 + 给定 TCP 上界；约束开关来自 comm_paras.yaml。"""
    lim = load_robot_limits(v_tcp_max=0.6, w_tcp_max=0.9)
    for arr in (lim.vmax, lim.amax, lim.jmax, lim.tau_max, lim.tau_min):
        assert np.asarray(arr).shape == (6,)
    assert np.all(lim.tau_min == -np.asarray(lim.tau_max))
    assert lim.v_tcp_max == 0.6 and lim.w_tcp_max == 0.9
    assert lim.a_bnd == (0.0, 0.0)
    # 约束开关：类型正确、未知键会报错、默认缺省为启用
    assert isinstance(flags, ConstraintFlags)
    assert ConstraintFlags.from_dict({}).jerk is True
    try:
        ConstraintFlags.from_dict({"bogus": True}); raised = False
    except ValueError:
        raised = True
    assert raised, "未知约束开关应报错"


def test_copp():
    """唯一端到端用例：求解 → 断言 → 由同一 (data, profile) 输出五张分析图（无 matplotlib 则跳过出图）。"""
    # ── 求解（约束开关取自 configs/comm_paras.yaml）────────────────────
    flags = load_constraint_flags()
    data = _make_data()
    profile, hist = solve_splp(data, SolveOptions(
        n_iter=6, flags=flags, smooth_c_weight=load_smooth_c_weight()))

    # ── 断言 ──────────────────────────────────────────────────────────
    _assert_config(flags)
    _assert_solver(data, profile, hist, flags)

    # ── 出图（全部来自唯一用例数据，落 self-test/output/，每次覆盖）──────
    try:
        from copp.viz import (
            plot_splp_result, plot_kinematic_limits, plot_speed_torque,
            plot_tn_convexification, plot_interp_profiles,
        )
    except ImportError:
        print("SKIP 出图部分（matplotlib 未安装）")
        return

    import matplotlib.pyplot as plt

    out_dir = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out_dir, exist_ok=True)

    # 图 1：SPLP 概览
    out1 = os.path.join(out_dir, "splp_test.png")
    fig1 = plot_splp_result(data, profile, hist, save_path=out1, show=False)
    assert os.path.exists(out1) and os.path.getsize(out1) > 0
    plt.close(fig1)

    # 图 2：关节 q̇/q̈/q⃛/力矩 + TCP 速度模（带约束）
    out2 = os.path.join(out_dir, "splp_limits_test.png")
    fig2 = plot_kinematic_limits(data, profile, tcp=_viz_tcp(data.s_grid),
                                 save_path=out2, show=False)
    assert os.path.exists(out2) and os.path.getsize(out2) > 0
    plt.close(fig2)

    # 图 3：速度相关力矩（t–n）分析——转矩–转速包络 / 约束利用率 / 摩擦分量
    if flags.speed_torque and data.speed_torque is not None:
        out3 = os.path.join(out_dir, "speed_torque_test.png")
        fig3 = plot_speed_torque(data, profile, save_path=out3, show=False)
        assert os.path.exists(out3) and os.path.getsize(out3) > 0
        plt.close(fig3)

        # 图 3b：论文 Fig.3 复现——逐关节 (q̇², τ) 平面的真实非凸可行域 vs 仿射切角内逼近
        out3b = os.path.join(out_dir, "tn_convexification_test.png")
        fig3b = plot_tn_convexification(data, profile, save_path=out3b, show=False)
        assert os.path.exists(out3b) and os.path.getsize(out3b) > 0
        plt.close(fig3b)

    # 图 4：论文 Fig.4 风格——由本用例 profile 全程解析重构 a/b/c/参数jerk（覆盖全部插补周期）
    out4 = os.path.join(out_dir, "fig4_interpolation.png")
    fig4 = plot_interp_profiles(data, profile, save_path=out4, show=False)
    assert os.path.exists(out4) and os.path.getsize(out4) > 0
    plt.close(fig4)


if __name__ == "__main__":
    test_copp()
    print("PASS test_copp")
