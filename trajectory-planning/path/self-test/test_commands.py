"""commands 指令层端到端测试（framework §8 / 设计 §3，M2）。

**单一测试用例** `test_commands`：
  1. 纯几何性质（不依赖机器人）：
     - JointMove：端点吻合、q' 恒定、零长报错；
     - LinearMove：直线性（残差 ~1e-15）、SLERP 端点吻合、ω̂ 恒定且与 R 的
       有限差分一致、纯姿态调整（Δp=0）安全、零长报错；
     - CircularMove：弧上各点到圆心距离恒为 r、端点/途经点吻合、三点模式方向
       由 via 唯一确定（劣弧/优弧都验）、(center,normal) 模式 ccw/cw/shortest
       语义、p'/p''/p''' 与有限差分一致、共线/重合/离面退化显式报错；
  2. 段间衔接：不衔接的序列抛 JunctionMismatchError（笛卡尔/关节两种）；
  3. 混合序列端到端：JointMove + LinearMove + CircularMove 在 UR5 上逐段
     降维 + SPLP 求解（M2 语义：段间停顿、每段 rest-to-rest），约束满足。

输出 output/commands_test.png。
可 pytest 运行，也可直接 `python trajectory-planning/path/self-test/test_commands.py`。
"""

from __future__ import annotations

import os
import sys

_PATH_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TRAJ_DIR = os.path.dirname(_PATH_DIR)
_REPO_ROOT = os.path.dirname(_TRAJ_DIR)
sys.path.insert(0, _TRAJ_DIR)
sys.path.insert(0, _REPO_ROOT)

import numpy as np

from robot import UR5Kinematics, Pose
from copp import load_robot_limits, solve_splp, SolveOptions
from copp.solve import s_to_t

from path.errors import (
    ZeroLengthCommandError, DegenerateArcError, JunctionMismatchError,
)
from path.commands import (
    JointMoveCommand, LinearMoveCommand, CircularMoveCommand,
    build_sections, lower_sections,
)
from path.commands.base import rotvec_between
from path.lowering import SampleOptions

_KIN = UR5Kinematics()
_Q_HOME = np.array([0.0, -np.pi / 2, np.pi / 2, -np.pi / 2, -np.pi / 2, 0.0])
_EYE = np.eye(3)


def _rot(rotvec) -> np.ndarray:
    from scipy.spatial.transform import Rotation
    return Rotation.from_rotvec(np.asarray(rotvec, dtype=float)).as_matrix()


def _fd_path(path, s: np.ndarray, h: float):
    """路径位置/姿态的中心差分（校验解析导数用）。"""
    sp, sm, s0 = path.eval(s + h), path.eval(s - h), path.eval(s)
    dp_fd = (sp.p - sm.p) / (2 * h)
    ddp_fd = (sp.dp - sm.dp) / (2 * h)
    dddp_fd = (sp.ddp - sm.ddp) / (2 * h)
    w_fd = np.zeros((3, s.size))
    for k in range(s.size):
        M = ((sp.R[k] - sm.R[k]) / (2 * h)) @ s0.R[k].T
        w_fd[:, k] = 0.5 * np.array([M[2, 1] - M[1, 2], M[0, 2] - M[2, 0], M[1, 0] - M[0, 1]])
    return dp_fd, ddp_fd, dddp_fd, w_fd


# ── 1a. JointMove 几何 ─────────────────────────────────────────────────────
def _assert_joint_move():
    q0 = _Q_HOME.copy()
    q1 = _Q_HOME + np.array([0.4, -0.3, 0.2, 0.3, -0.2, 0.5])
    sec = JointMoveCommand(q0, q1).to_section()
    assert sec.native_space == "joint"
    L = sec.path.s_total
    assert abs(L - np.linalg.norm(q1 - q0)) < 1e-12
    s = np.linspace(0.0, L, 11)
    q, dq, ddq, dddq = sec.path.eval_joint(s)
    assert np.allclose(q[:, 0], q0) and np.allclose(q[:, -1], q1)
    assert np.allclose(dq - dq[:, :1], 0.0), "线性关节几何 q' 应恒定"
    assert np.allclose(ddq, 0.0) and np.allclose(dddq, 0.0)
    try:
        JointMoveCommand(q0, q0).to_section()
        raise AssertionError("零长 JointMove 未报错")
    except ZeroLengthCommandError:
        pass


# ── 1b. LinearMove 几何 ────────────────────────────────────────────────────
def _assert_linear_move():
    p0 = np.array([0.3, -0.2, 0.4])
    p1 = np.array([0.42, -0.05, 0.33])
    R0 = _rot([0.3, -0.2, 0.1])
    R1 = R0 @ _rot([0.15, 0.25, -0.2])
    sec = LinearMoveCommand(Pose(p0, R0), Pose(p1, R1)).to_section()
    path = sec.path
    L = path.s_total
    assert abs(L - np.linalg.norm(p1 - p0)) < 1e-12

    s = np.linspace(0.0, L, 21)
    smp = path.eval(s)
    assert np.allclose(smp.p[:, 0], p0) and np.allclose(smp.p[:, -1], p1)
    assert np.allclose(smp.R[0], R0, atol=1e-12) and np.allclose(smp.R[-1], R1, atol=1e-12)
    chord = (p1 - p0) / L
    dev = smp.p - p0[:, None] - chord[:, None] * s[None, :]
    assert np.max(np.abs(dev)) < 1e-12, "直线性失败"
    rotvec = rotvec_between(R0, R1)
    assert np.allclose(smp.w, (R0 @ rotvec)[:, None] / L), "ω̂ 应恒为 R0Θ/L"

    s_in = np.linspace(0.05 * L, 0.95 * L, 9)
    dp_fd, _, _, w_fd = _fd_path(path, s_in, 1e-7 * L)
    smp_in = path.eval(s_in)
    assert np.allclose(dp_fd, smp_in.dp, rtol=1e-6, atol=1e-9)
    assert np.allclose(w_fd, smp_in.w, rtol=1e-5, atol=1e-8), "SLERP ω̂ 与 R 的 FD 不符"

    # 纯姿态调整：Δp=0 → L 由转角×rot_scale 决定
    sec_rot = LinearMoveCommand(Pose(p0, R0), Pose(p0.copy(), R1), rot_scale=0.1).to_section()
    theta = float(np.linalg.norm(rotvec_between(R0, R1)))
    assert abs(sec_rot.path.s_total - 0.1 * theta) < 1e-12
    assert np.allclose(sec_rot.path.eval(np.array([0.0])).dp, 0.0)
    try:
        LinearMoveCommand(Pose(p0, R0), Pose(p0.copy(), R0.copy())).to_section()
        raise AssertionError("零长 LinearMove 未报错")
    except ZeroLengthCommandError:
        pass


# ── 1c. CircularMove 几何 ──────────────────────────────────────────────────
def _assert_circular_move():
    r = 0.12
    c0 = np.array([0.1, 0.2, 0.3])
    p0 = c0 + np.array([r, 0.0, 0.0])
    p1 = c0 + np.array([0.0, r, 0.0])

    # 三点定圆（劣弧：via 在 45°）
    via = c0 + r * np.array([np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0])
    sec = CircularMoveCommand(Pose(p0, _EYE), Pose(p1, _EYE.copy()), via=via).to_section()
    path = sec.path
    assert abs(path.sweep - np.pi / 2) < 1e-9, f"劣弧扫角应为 π/2，得 {path.sweep}"
    assert abs(path.s_total - r * np.pi / 2) < 1e-9
    s = np.linspace(0.0, path.s_total, 33)
    smp = path.eval(s)
    assert np.allclose(np.linalg.norm(smp.p - path.center[:, None], axis=0), r, atol=1e-9), \
        "弧上点到圆心距离应恒为 r"
    assert np.allclose(smp.p[:, 0], p0) and np.allclose(smp.p[:, -1], p1, atol=1e-9)
    assert np.min(np.linalg.norm(smp.p - via[:, None], axis=0)) < r * 0.08, "弧未经过 via"

    # 优弧：via 在 225° → 走长边，|扫角|=3π/2（sweep 符号相对 via 定出的法向
    # n̂=a×b，via 在下半平面时 n̂=−ẑ，故 (e1,e2) 系内 sweep 为 +3π/2——几何不变量
    # 是弧长与"经过 via"，不是世界系符号）
    via_far = c0 + r * np.array([np.cos(1.25 * np.pi), np.sin(1.25 * np.pi), 0.0])
    sec2 = CircularMoveCommand(Pose(p0, _EYE), Pose(p1, _EYE.copy()), via=via_far).to_section()
    assert abs(abs(sec2.path.sweep) - 1.5 * np.pi) < 1e-9, f"优弧|扫角|应为 3π/2，得 {sec2.path.sweep}"
    assert abs(sec2.path.s_total - r * 1.5 * np.pi) < 1e-9
    s2 = np.linspace(0.0, sec2.path.s_total, 65)
    smp2 = sec2.path.eval(s2)
    assert np.min(np.linalg.norm(smp2.p - via_far[:, None], axis=0)) < r * 0.05, "优弧未经过 via"
    assert np.allclose(smp2.p[:, -1], p1, atol=1e-9)

    # (center, normal) 模式方向语义
    for direction, sweep_ref in (("ccw", np.pi / 2), ("cw", np.pi / 2 - 2 * np.pi),
                                 ("shortest", np.pi / 2)):
        sec3 = CircularMoveCommand(Pose(p0, _EYE), Pose(p1, _EYE.copy()),
                                   center=c0, normal=np.array([0.0, 0.0, 1.0]),
                                   direction=direction).to_section()
        assert abs(sec3.path.sweep - sweep_ref) < 1e-9, f"{direction} 扫角错误"

    # 解析导数 vs 有限差分
    s_in = np.linspace(0.05 * path.s_total, 0.95 * path.s_total, 9)
    dp_fd, ddp_fd, dddp_fd, w_fd = _fd_path(path, s_in, 1e-7 * path.s_total)
    smp_in = path.eval(s_in)
    assert np.allclose(dp_fd, smp_in.dp, rtol=1e-6, atol=1e-8)
    assert np.allclose(ddp_fd, smp_in.ddp, rtol=1e-5, atol=1e-7)
    assert np.allclose(dddp_fd, smp_in.dddp, rtol=1e-4, atol=1e-5)
    assert np.allclose(w_fd, smp_in.w, rtol=1e-5, atol=1e-8)
    assert np.allclose(np.linalg.norm(smp_in.dp, axis=0), 1.0, atol=1e-12), "弧长参数化 ‖p'‖=1"
    assert np.allclose(np.linalg.norm(smp_in.ddp, axis=0), 1.0 / r, atol=1e-12), "曲率 ‖p''‖=1/r"

    # 退化情形
    for bad in (
        lambda: CircularMoveCommand(Pose(p0, _EYE), Pose(p1, _EYE.copy()),
                                    via=0.5 * (p0 + p1)).to_section(),          # 共线
        lambda: CircularMoveCommand(Pose(p0, _EYE), Pose(Pose(p0, _EYE).position.copy(), _EYE.copy()),
                                    via=via).to_section(),                       # 起终点重合
        lambda: CircularMoveCommand(Pose(p0, _EYE), Pose(p1 + np.array([0, 0, 0.05]), _EYE.copy()),
                                    center=c0, normal=np.array([0.0, 0.0, 1.0])).to_section(),  # 离面
    ):
        try:
            bad()
            raise AssertionError("退化圆弧未报错")
        except DegenerateArcError:
            pass


# ── 2. 段间衔接校验 ────────────────────────────────────────────────────────
def _assert_junction_checks():
    qA = _Q_HOME + 0.3 * np.array([1.0, -0.6, 0.5, 0.6, -0.4, 1.0])
    poseA = _KIN.fk(qA)
    poseB = Pose(poseA.position + np.array([0.06, 0.08, -0.04]),
                 poseA.rotation @ _rot([0.1, -0.05, 0.08]))

    # 笛卡尔段起点位姿与上游不衔接（偏 2mm）
    bad_pose = Pose(poseA.position + np.array([0.002, 0.0, 0.0]), poseA.rotation.copy())
    secs = build_sections([
        JointMoveCommand(_Q_HOME, qA),
        LinearMoveCommand(bad_pose, poseB),
    ])
    try:
        lower_sections(secs, _KIN, q_seed=_Q_HOME)
        raise AssertionError("笛卡尔段衔接错位未报错")
    except JunctionMismatchError:
        pass

    # 关节段起点与上游不衔接
    secs2 = build_sections([
        JointMoveCommand(_Q_HOME, qA),
        JointMoveCommand(qA + 0.01, _Q_HOME),
    ])
    try:
        lower_sections(secs2, _KIN, q_seed=_Q_HOME)
        raise AssertionError("关节段衔接错位未报错")
    except JunctionMismatchError:
        pass


# ── 3. 混合序列端到端 ──────────────────────────────────────────────────────
def _run_mixed_pipeline():
    qA = _Q_HOME + 0.3 * np.array([1.0, -0.6, 0.5, 0.6, -0.4, 1.0])
    poseA = _KIN.fk(qA)
    pB = poseA.position + np.array([0.06, 0.08, -0.04])
    RB = poseA.rotation @ _rot([0.1, -0.05, 0.08])
    poseB = Pose(pB, RB)
    pC = pB + np.array([0.08, -0.05, 0.03])
    chord = pC - pB
    perp = np.cross(chord, np.array([0.0, 0.0, 1.0]))
    perp = perp / np.linalg.norm(perp)
    via = 0.5 * (pB + pC) + 0.03 * perp
    poseC = Pose(pC, RB @ _rot([-0.08, 0.12, 0.05]))

    commands = [
        JointMoveCommand(_Q_HOME, qA),
        LinearMoveCommand(poseA, poseB),
        CircularMoveCommand(poseB, poseC, via=via),
    ]
    sections = build_sections(commands)
    pds = lower_sections(sections, _KIN, q_seed=_Q_HOME,
                         sample_opts=SampleOptions(eps_pos=1e-4, eps_ori=1e-3))
    assert len(pds) == 3

    limits = load_robot_limits(v_tcp_max=0.5, w_tcp_max=2.5)
    results = []
    for pd in pds:
        data = limits.to_topp3_data(pd.s_grid, pd.dq, pd.ddq, pd.dddq,
                                    tcp_geom=pd.tcp_geom())
        profile, hist = solve_splp(data, SolveOptions(n_iter=3))
        tf = hist.t_final
        assert all(x >= y - 1e-6 for x, y in zip(tf, tf[1:])), "SPLP 非单调"
        assert np.isfinite(tf[-1]) and tf[-1] > 0
        a, b = profile.a, profile.b
        sa = np.sqrt(np.maximum(a, 0.0))
        vmax, amax, _ = limits.axis_arrays(pd.n_axis)
        assert np.all(np.abs(pd.dq * sa[None, :]) <= vmax[:, None] * 1.02 + 1e-9)
        assert np.all(np.abs(pd.ddq * a[None, :] + pd.dq * b[None, :])
                      <= amax[:, None] * 1.05 + 1e-9)
        assert np.all(pd.cv * sa <= limits.v_tcp_max * 1.02 + 1e-9)
        results.append((pd, profile, hist))

    # 段间关节角衔接（IK 精确 → 应严丝合缝）
    assert np.max(np.abs(pds[0].q[:, -1] - pds[1].q[:, 0])) < 1e-8
    assert np.max(np.abs(pds[1].q[:, -1] - pds[2].q[:, 0])) < 1e-8
    return results, via


def test_commands():
    _assert_joint_move()
    _assert_linear_move()
    _assert_circular_move()
    _assert_junction_checks()
    results, via = _run_mixed_pipeline()
    _plot(results, via)
    print("PASS test_commands")


def _plot(results, via):
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
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle("commands 指令层：JointMove + LinearMove + CircularMove → 逐段 SPLP（M2，段间停顿）",
                 fontsize=12)
    names = ["JointMove", "LinearMove", "CircularMove"]

    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    for i, (pd, _, _) in enumerate(results):
        tcp = np.stack([_KIN.fk(pd.q[:, k]).position for k in range(pd.n_grid)], axis=1)
        ax1.plot(tcp[0], tcp[1], tcp[2], lw=1.6, label=names[i])
    ax1.scatter(*via, color="k", s=25, label="via")
    ax1.set_title("① TCP 轨迹（三段）"); ax1.legend(fontsize=7)

    ax2 = fig.add_subplot(1, 3, 2)
    for i, (pd, profile, _) in enumerate(results):
        ax2.plot(pd.s_grid / pd.s_grid[-1] + i, np.sqrt(np.maximum(profile.a, 0.0)),
                 lw=1.4, label=names[i])
    ax2.set_xlabel("归一化 s + 段号"); ax2.set_title(r"② 各段速度剖面 $\dot s(s)$（rest-to-rest）")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.25)

    ax3 = fig.add_subplot(1, 3, 3)
    t_off = 0.0
    for pd, profile, _ in results:
        t_final, t_s = s_to_t(pd.s_grid, profile)
        sa = np.sqrt(np.maximum(profile.a, 0.0))
        for i in range(6):
            ax3.plot(t_off + t_s, pd.dq[i] * sa, lw=1.0,
                     color=f"C{i}", label=f"关节 {i}" if t_off == 0.0 else None)
        t_off += t_final
    ax3.set_xlabel("时间 t (s)"); ax3.set_title(r"③ 关节速度 $\dot q(t)$（三段拼接）")
    ax3.legend(fontsize=7, ncol=2); ax3.grid(alpha=0.25)

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(out_dir, "commands_test.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    assert os.path.exists(out) and os.path.getsize(out) > 0


if __name__ == "__main__":
    test_commands()
