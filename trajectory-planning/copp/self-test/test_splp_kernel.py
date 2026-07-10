"""M1 SPLP 数值内核端到端测试（framework §8）。

**单一测试用例**：一次求解 + 全部断言 + 同时输出三张分析图，便于对照分析问题：
  output/splp_test.png          —— SPLP 概览（收敛 / 速度剖面 / 约束利用率 / 时间律…）
  output/splp_limits_test.png   —— 关节速度/加速度/jerk/力矩 + TCP 速度模（带约束）
  output/fig4_interpolation.png —— 论文 Fig.4 复现（区间解析插值 Prop.1 + Prop.2）

可用 pytest 运行，也可直接 `python trajectory-planning/copp/self-test/test_splp_kernel.py`。
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
    load_robot_limits, load_fig4_example, load_constraint_flags,
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
        r_v = np.max(np.abs(qd) / data.vmax[:, None]); assert r_v <= 1.02, r_v
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


def _assert_fig4(ex, cfg):
    """Fig.4 示意：静止起点 + 非静止终点 + 头部恒定参数 jerk。"""
    assert ex["a_grid"][0] == 0.0 and ex["b_grid"][0] == 0.0
    assert ex["a_grid"][-1] > 0.0 and ex["b_grid"][-1] != 0.0
    head = ex["j_fine"][ex["u_fine"] < ex["u_stat"] - 1e-9]  # 头尾共享 u_stat，严格小于只取头部
    assert np.allclose(head, ex["kappa"], rtol=1e-6)
    assert abs(ex["a_grid"][cfg["n_stat"]] - cfg["a_head"]) < 1e-9


def test_splp_kernel():
    """端到端单一用例：求解 + 断言 + 输出三张分析图（无 matplotlib 则仅跳过出图）。"""
    # ── 求解（约束开关取自 configs/comm_paras.yaml）────────────────────
    flags = load_constraint_flags()
    data = _make_data()
    profile, hist = solve_splp(data, SolveOptions(n_iter=6, flags=flags))

    # ── 断言 ──────────────────────────────────────────────────────────
    _assert_config(flags)
    _assert_solver(data, profile, hist, flags)
    cfg = load_fig4_example()

    # ── 出图（三张，落 self-test/output/，每次覆盖）────────────────────
    try:
        from copp.viz import (
            plot_splp_result, plot_kinematic_limits,
            fig4_interpolation_example, plot_fig4_interpolation,
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

    # 图 3：论文 Fig.4 复现（区间解析插值）——参数取自 comm_paras.yaml
    ex = fig4_interpolation_example(**cfg)
    _assert_fig4(ex, cfg)
    out3 = os.path.join(out_dir, "fig4_interpolation.png")
    fig3 = plot_fig4_interpolation(ex, save_path=out3, show=False)
    assert os.path.exists(out3) and os.path.getsize(out3) > 0
    plt.close(fig3)


if __name__ == "__main__":
    test_splp_kernel()
    print("PASS test_splp_kernel")
