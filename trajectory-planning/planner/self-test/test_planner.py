"""planner 门面 + synth 合成层端到端测试（framework §5.8/§5.10 / 设计 §9-10，M2+）。

**单一测试用例** `test_planner`：JointMove + LinearMove + CircularMove 三段混合
指令经 `TrajectoryPlanner.plan()` 一次跑通（降维 → 逐段 rest-to-rest SPLP →
等时间栅格合成 → 拼接 → R_v/D_v 校验），断言：

  1. 结构：三段、总时长 = 各段之和、时间严格递增、段内等距 dt；
  2. 物理连续性：起点/终点/段 seam 处 rest（q̇≈0）、q 无跳变、终点 FK 命中目标；
  3. 独立一致性：s(t) 的时间差分 ≈ ṡ（解析时间律交叉验证）；q 的时间差分 ≈ q̇；
  4. 约束校验器：真实限值下 R_v=D_v=0（论文 <0.1% 目标）；人为收紧限值后
     verify_limits 必须报警（校验器本身的灵敏度测试）；
  5. 防误用：空指令队列报错；非静止边界的多段规划报错。

输出 output/planner_test.png。
可 pytest 运行，也可直接 `python trajectory-planning/planner/self-test/test_planner.py`。
"""

from __future__ import annotations

import os
import sys

_PLANNER_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TRAJ_DIR = os.path.dirname(_PLANNER_DIR)
_REPO_ROOT = os.path.dirname(_TRAJ_DIR)
sys.path.insert(0, _TRAJ_DIR)     # 供 import copp / path / planner
sys.path.insert(0, _REPO_ROOT)    # 供 import robot

import numpy as np

from robot import UR5Kinematics, Pose
from copp import load_robot_limits

from path.commands import JointMoveCommand, LinearMoveCommand, CircularMoveCommand
from planner import TrajectoryPlanner, PlanOptions, verify_limits

_KIN = UR5Kinematics()
_Q_HOME = np.array([0.0, -np.pi / 2, np.pi / 2, -np.pi / 2, -np.pi / 2, 0.0])
_DT = 1e-3


def _rot(rotvec) -> np.ndarray:
    from scipy.spatial.transform import Rotation
    return Rotation.from_rotvec(np.asarray(rotvec, dtype=float)).as_matrix()


def _make_commands():
    """与 commands 自测同款三段几何：Joint→Linear→Circular（G0 精确衔接）。"""
    qA = _Q_HOME + 0.3 * np.array([1.0, -0.6, 0.5, 0.6, -0.4, 1.0])
    poseA = _KIN.fk(qA)
    pB = poseA.position + np.array([0.06, 0.08, -0.04])
    poseB = Pose(pB, poseA.rotation @ _rot([0.1, -0.05, 0.08]))
    pC = pB + np.array([0.08, -0.05, 0.03])
    chord = pC - pB
    perp = np.cross(chord, np.array([0.0, 0.0, 1.0]))
    perp = perp / np.linalg.norm(perp)
    via = 0.5 * (pB + pC) + 0.03 * perp
    poseC = Pose(pC, poseB.rotation @ _rot([-0.08, 0.12, 0.05]))
    return [
        JointMoveCommand(_Q_HOME, qA),
        LinearMoveCommand(poseA, poseB),
        CircularMoveCommand(poseB, poseC, via=via),
    ], poseC


def _assert_structure(res, n_seg: int):
    traj = res.trajectory
    assert len(res.segments) == n_seg
    assert traj.t[0] == 0.0 and np.all(np.diff(traj.t) > 0.0), "时间应严格递增"
    tf_sum = sum(sp.t_final for sp in res.segments)
    assert abs(traj.t_final - tf_sum) < 1e-9, "总时长应为各段之和"
    assert abs(traj.t[-1] - traj.t_final) < 1e-12, "末样本应落在精确 t_final"
    # 段内等距 dt（排除各段最后一步的残步与 seam）
    for i in range(n_seg):
        idx = np.where(traj.seg_index == i)[0]
        dts = np.diff(traj.t[idx])
        if dts.size > 2:
            assert np.allclose(dts[:-1], _DT, atol=1e-12), f"段 {i} 采样间距非 dt"
    # SPLP 单调
    for sp in res.segments:
        tfs = sp.splp_t_final
        assert all(x >= y - 1e-6 for x, y in zip(tfs, tfs[1:])), "SPLP 非单调"


def _assert_continuity(res, pose_end: Pose):
    traj = res.trajectory
    # 起点/终点 rest（rest 端 a 有 1e-12 数值下限 → ṡ~1e-6，阈值取 1e-4 rad/s）
    assert np.max(np.abs(traj.qd[:, 0])) < 1e-4, "起点应静止"
    assert np.max(np.abs(traj.qd[:, -1])) < 1e-4, "终点应静止"
    assert np.max(np.abs(traj.q[:, 0] - _Q_HOME)) < 1e-9
    fk_end = _KIN.fk(traj.q[:, -1])
    assert np.linalg.norm(fk_end.position - pose_end.position) < 1e-6, "终点位置未命中"
    assert np.linalg.norm(fk_end.rotation - pose_end.rotation) < 1e-6, "终点姿态未命中"
    # seam 处 rest + q 连续（最后样本 vs 下一段首样本已在拼接中去重，
    # 这里检查跨段相邻样本）
    seam = np.where(np.diff(res.trajectory.seg_index) != 0)[0]
    for j in seam:
        assert np.max(np.abs(traj.qd[:, j])) < 1e-3, "seam 前样本应近静止"
        assert np.max(np.abs(traj.qd[:, j + 1])) < 1e-3, "seam 后样本应近静止"
        assert np.max(np.abs(traj.q[:, j + 1] - traj.q[:, j])) < 1e-3, "seam 处 q 跳变"
    # 全程无关节跳变
    assert np.max(np.abs(np.diff(traj.q, axis=1))) < 0.05, "相邻样本 q 跳变过大"


def _assert_consistency(res):
    """独立交叉验证：FD(s)≈ṡ、FD(q)≈q̇（逐段内部样本，避开 seam）。"""
    for sp in res.segments:
        r = sp.result
        if r.t.size < 8:
            continue
        inner = slice(2, -2)
        ds_dt = np.gradient(r.s, r.t)
        err_s = np.max(np.abs(ds_dt[inner] - r.sdot[inner]))
        assert err_s < 2e-2 * (1.0 + np.max(r.sdot)), f"FD(s) vs ṡ 偏差 {err_s:.3e}"
        dq_dt = np.gradient(r.q, r.t, axis=1)
        err_q = np.max(np.abs(dq_dt[:, inner] - r.qd[:, inner]))
        assert err_q < 5e-2 * (1.0 + np.max(np.abs(r.qd))), f"FD(q) vs q̇ 偏差 {err_q:.3e}"


def _assert_verify(res, limits):
    m = res.metrics
    assert m is not None
    assert m.ok, f"真实限值下应无超限：{m.summary()}"
    assert m.r_v == 0.0 and m.d_v == 0.0, m.summary()
    for key in ("velocity", "acceleration", "tcp_velocity", "tcp_angular_velocity"):
        assert m.max_util[key] <= 1.0 + 1e-6, f"{key} 利用率 {m.max_util[key]:.4f} 超 1"
    assert m.max_util["jerk"] <= 1.05, f"jerk 利用率 {m.max_util['jerk']:.4f}（区间内 O(Δ²) 界）"

    # 校验器灵敏度：人为把速度上限压到 1/4 → 必须报警
    tight = load_robot_limits(v_tcp_max=limits.v_tcp_max, w_tcp_max=limits.w_tcp_max)
    tight.vmax = np.asarray(limits.axis_arrays(6)[0]) / 4.0
    m_bad = verify_limits(res.trajectory, tight)
    assert not m_bad.ok and m_bad.r_v > 0.0, "收紧限值后 verify 应报警"


def _assert_guards(limits):
    try:
        TrajectoryPlanner(_KIN, limits).plan(q_seed=_Q_HOME)
        raise AssertionError("空指令队列未报错")
    except ValueError:
        pass
    bad_limits = load_robot_limits(v_tcp_max=0.5, w_tcp_max=2.5)
    bad_limits.a_bnd = (0.04, 0.0)
    cmds, _ = _make_commands()
    planner = TrajectoryPlanner(_KIN, bad_limits)
    for c in cmds:
        planner.add_command(c)
    try:
        planner.plan(q_seed=_Q_HOME)
        raise AssertionError("非静止边界的多段规划未报错")
    except ValueError:
        pass


def test_planner():
    limits = load_robot_limits(v_tcp_max=0.5, w_tcp_max=2.5)
    commands, pose_end = _make_commands()

    planner = TrajectoryPlanner(_KIN, limits)
    for c in commands:
        planner.add_command(c)
    res = planner.plan(q_seed=_Q_HOME, opts=PlanOptions(dt=_DT))

    _assert_structure(res, n_seg=3)
    _assert_continuity(res, pose_end)
    _assert_consistency(res)
    _assert_verify(res, limits)
    _assert_guards(limits)

    print(f"t_final = {res.t_final:.4f} s; {res.metrics.summary()}")
    _plot(res, limits)
    print("PASS test_planner")


def _plot(res, limits):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        matplotlib.rcParams["font.sans-serif"] = [
            "Microsoft YaHei", "SimHei", "DengXian", "SimSun"]
        matplotlib.rcParams["axes.unicode_minus"] = False
    except ImportError:
        print("SKIP 出图（matplotlib 未安装）")
        return

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(out_dir, exist_ok=True)
    traj = res.trajectory
    t = traj.t
    vmax, amax, jmax = limits.axis_arrays(traj.n_axis)
    seams = t[np.where(np.diff(traj.seg_index) != 0)[0] + 1]

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(
        f"planner 门面：三段混合指令 → 等时间栅格轨迹（M2+ 段间停顿；t_f={traj.t_final:.3f}s）",
        fontsize=13)

    def _mark_seams(axis):
        for ts in seams:
            axis.axvline(ts, color="0.75", lw=0.8, ls=":")

    for i in range(traj.n_axis):
        ax[0, 0].plot(t, traj.q[i], lw=1.1, label=f"关节 {i}")
    ax[0, 0].set_title("① 关节角 q(t)"); ax[0, 0].legend(fontsize=7, ncol=2)

    for i in range(traj.n_axis):
        color = f"C{i}"
        ax[0, 1].plot(t, traj.qd[i], lw=1.1, color=color)
        ax[0, 1].axhline(vmax[i], color=color, ls="--", lw=0.7, alpha=0.5)
        ax[0, 1].axhline(-vmax[i], color=color, ls="--", lw=0.7, alpha=0.5)
    ax[0, 1].set_title("② 关节速度 q̇(t)（虚线=各关节 vmax）")

    for i in range(traj.n_axis):
        color = f"C{i}"
        ax[0, 2].plot(t, traj.qdd[i], lw=1.1, color=color)
        ax[0, 2].axhline(amax[i], color=color, ls="--", lw=0.7, alpha=0.5)
        ax[0, 2].axhline(-amax[i], color=color, ls="--", lw=0.7, alpha=0.5)
    ax[0, 2].set_title("③ 关节加速度 q̈(t)（虚线=各关节 amax）")

    ax[1, 0].plot(t, traj.sdot, color="C0")
    ax[1, 0].set_title(r"④ 路径速度 $\dot s(t)$（逐段 rest-to-rest）")

    ax[1, 1].plot(t, traj.v_tcp, color="C0", label=r"$\|\dot p\|$")
    ax[1, 1].axhline(limits.v_tcp_max, color="C0", ls="--", lw=0.9, label="v_tcp_max")
    ax[1, 1].plot(t, traj.w_tcp, color="C4", label=r"$\|\omega\|$")
    ax[1, 1].axhline(limits.w_tcp_max, color="C4", ls="--", lw=0.9, label="w_tcp_max")
    ax[1, 1].set_title("⑤ TCP 速度模（带约束）"); ax[1, 1].legend(fontsize=8)

    for i in range(traj.n_axis):
        ax[1, 2].plot(t, traj.qddd[i], lw=0.9, color=f"C{i}")
    ax[1, 2].axhline(float(np.max(jmax)), color="0.4", ls="--", lw=0.9)
    ax[1, 2].axhline(-float(np.max(jmax)), color="0.4", ls="--", lw=0.9)
    ax[1, 2].set_title("⑥ 关节 jerk q⃛(t)")

    for a_ in ax.flat:
        _mark_seams(a_)
        a_.grid(alpha=0.25)
        a_.set_xlabel("t (s)")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(out_dir, "planner_test.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    assert os.path.exists(out) and os.path.getsize(out) > 0


if __name__ == "__main__":
    test_planner()
