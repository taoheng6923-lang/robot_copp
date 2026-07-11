"""lowering 降维层端到端测试（framework §8 / 设计 §5，M2）。

**单一测试用例** `test_lowering`：解析笛卡尔位姿路径（位置三角曲线 + 变轴姿态
Rx(α(s))·Ry(β(s))，ω̂/ω̂'/ω̂'' 全解析）→ 自适应采样 → UR5 连续解 IK →
Jacobian 链式法则求 q',q'',q''' → 有限差分交叉验证 → 喂给 copp SPLP 求解。

验证链条（每步都有独立的数值交叉验证，不依赖被测代码自身）：
  1. 夹具自检：路径的解析 w/dw/ddw 与 R(s) 的有限差分一致（防夹具推导错误）；
  2. IK 回代：FK(q_k) 复现路径位姿（~1e-9）；
  3. 链式法则：q 序列的有限差分复现 dq/ddq/dddq（尤其验证三阶式 J'q'' 的系数 2——
     设计文档 §5.3 原式漏了该系数，若按原式实现此处会差 O(10%) 量级）；
  4. 独立 Jacobian 恒等式：网格差分的 J' 满足 J'q' + Jq'' ≈ r₂；
  5. 自适应采样：网格含段边界、弦高误差在容差内；
  6. 纯姿态调整（p'≡0）退化安全：cv≡0、IK/求导不崩；
  7. 端到端：PathDerivatives → Topp3Data → solve_splp，约束满足、迭代单调。

输出 output/lowering_test.png（六面板概览）。
可 pytest 运行，也可直接 `python trajectory-planning/path/self-test/test_lowering.py`。
"""

from __future__ import annotations

import os
import sys

_PATH_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TRAJ_DIR = os.path.dirname(_PATH_DIR)
_REPO_ROOT = os.path.dirname(_TRAJ_DIR)
sys.path.insert(0, _TRAJ_DIR)     # 供 import copp / path
sys.path.insert(0, _REPO_ROOT)    # 供 import robot

from dataclasses import dataclass, field

import numpy as np

from robot import UR5Kinematics
from copp import load_robot_limits, solve_splp, SolveOptions

from path.types import CartesianSamples
from path.lowering import (
    SampleOptions, adaptive_sample,
    lower_cartesian, lower_joint,
    min_singular_ratio, damped_inverse_solve,
)

_KIN = UR5Kinematics()
_Q_HOME = np.array([0.0, -np.pi / 2, np.pi / 2, -np.pi / 2, -np.pi / 2, 0.0])
_X = np.array([1.0, 0.0, 0.0])
_Y = np.array([0.0, 1.0, 0.0])


def _rx(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _ry(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


@dataclass
class _RichPath:
    """解析测试路径：位置三角曲线 + 变轴姿态 R0·Rx(α(u))·Ry(β(u))，u=s/L。

    α=ra·sin(πu)、β=rb·(3u²−2u³)；ω̂ = R0·(α_s x̂ + β_s Rx(α)ŷ)（世界系），
    ω̂'/ω̂'' 按乘积法则解析展开（见各行注释），供链式法则的 r₂/r₃ 输入。
    """

    p0: np.ndarray
    R0: np.ndarray
    L: float = 0.3
    amp: float = 0.04
    om: float = 1.5 * np.pi
    ra: float = 0.30
    rb: float = 0.25
    s_breaks: np.ndarray = field(default=None)  # type: ignore[assignment]

    def __post_init__(self):
        self.s_total = float(self.L)
        if self.s_breaks is None:
            self.s_breaks = np.array([0.0, self.s_total])

    def eval(self, s: np.ndarray) -> CartesianSamples:
        s = np.atleast_1d(np.asarray(s, dtype=float))
        N = s.size
        u = s / self.L
        A, om, L = self.amp, self.om, self.L

        # ── 位置（对 u 解析求导后除 L^k）─────────────────────────────────
        t1, t2 = om * u, 2.0 * om * u + 0.7
        p = self.p0[:, None] + A * np.vstack([np.sin(t1), np.sin(t2), np.cos(t1) - 1.0])
        dp = A * np.vstack([om * np.cos(t1), 2 * om * np.cos(t2), -om * np.sin(t1)]) / L
        ddp = A * np.vstack([-om**2 * np.sin(t1), -4 * om**2 * np.sin(t2), -om**2 * np.cos(t1)]) / L**2
        dddp = A * np.vstack([-om**3 * np.cos(t1), -8 * om**3 * np.cos(t2), om**3 * np.sin(t1)]) / L**3

        # ── 姿态角及 s 导 ────────────────────────────────────────────────
        al = self.ra * np.sin(np.pi * u)
        al_s = self.ra * np.pi * np.cos(np.pi * u) / L
        al_ss = -self.ra * np.pi**2 * np.sin(np.pi * u) / L**2
        al_sss = -self.ra * np.pi**3 * np.cos(np.pi * u) / L**3
        be_s = self.rb * (6 * u - 6 * u**2) / L
        be_ss = self.rb * (6 - 12 * u) / L**2
        be_sss = np.full(N, -12.0 * self.rb / L**3)
        be = self.rb * (3 * u**2 - 2 * u**3)

        R = np.zeros((N, 3, 3))
        w = np.zeros((3, N))
        dw = np.zeros((3, N))
        ddw = np.zeros((3, N))
        for k in range(N):
            Rx = _rx(al[k])
            y = Rx @ _Y                      # Rx(α)ŷ
            xy = np.cross(_X, y)             # x̂×(Rxŷ)  ← d(Rxŷ)/ds = α_s·x̂×(Rxŷ)
            xxy = np.cross(_X, xy)           # x̂×(x̂×(Rxŷ))
            m = al_s[k] * _X + be_s[k] * y
            m1 = al_ss[k] * _X + be_ss[k] * y + be_s[k] * al_s[k] * xy
            m2 = (al_sss[k] * _X + be_sss[k] * y
                  + (2 * be_ss[k] * al_s[k] + be_s[k] * al_ss[k]) * xy
                  + be_s[k] * al_s[k] ** 2 * xxy)
            R[k] = self.R0 @ Rx @ _ry(be[k])
            w[:, k] = self.R0 @ m
            dw[:, k] = self.R0 @ m1
            ddw[:, k] = self.R0 @ m2
        return CartesianSamples(p=p, dp=dp, ddp=ddp, dddp=dddp, R=R, w=w, dw=dw, ddw=ddw)


@dataclass
class _PureRotPath:
    """纯姿态调整：p 恒定、绕固定世界轴匀速转 θ。p'≡0 的退化安全用例。"""

    p0: np.ndarray
    R0: np.ndarray
    theta: float = 0.8
    L: float = 0.1              # 特征长度（rot_scale·θ 的角色）
    s_breaks: np.ndarray = field(default=None)  # type: ignore[assignment]

    def __post_init__(self):
        self.s_total = float(self.L)
        if self.s_breaks is None:
            self.s_breaks = np.array([0.0, self.s_total])

    def eval(self, s: np.ndarray) -> CartesianSamples:
        s = np.atleast_1d(np.asarray(s, dtype=float))
        N = s.size
        zero3 = np.zeros((3, N))
        R = np.zeros((N, 3, 3))
        w = np.zeros((3, N))
        for k in range(N):
            R[k] = self.R0 @ _rx(self.theta * s[k] / self.L)
            w[:, k] = self.R0 @ (_X * self.theta / self.L)   # Ṙ Rᵀ = [R0·(θ/L)x̂]×
        return CartesianSamples(
            p=np.repeat(self.p0[:, None], N, axis=1),
            dp=zero3, ddp=zero3.copy(), dddp=zero3.copy(),
            R=R, w=w, dw=zero3.copy(), ddw=zero3.copy(),
        )


def _vee(M: np.ndarray) -> np.ndarray:
    return 0.5 * np.array([M[2, 1] - M[1, 2], M[0, 2] - M[2, 0], M[1, 0] - M[0, 1]])


def _fd(arr: np.ndarray, s: np.ndarray) -> np.ndarray:
    """(n,N) 数组沿均匀网格的中心差分（端点回退为一侧，校验时只取内部）。"""
    return np.gradient(arr, s, axis=1)


def _assert_fixture_selfcheck(path) -> None:
    """夹具自检：解析 w/dw/ddw 与 R(s)/w(s) 的有限差分一致；R 正交。"""
    L = path.s_total
    s = np.linspace(0.02 * L, 0.98 * L, 41)
    h = 1e-6 * L
    smp = path.eval(s)
    for k in range(s.size):
        assert np.allclose(smp.R[k] @ smp.R[k].T, np.eye(3), atol=1e-12)
    sp, sm = path.eval(s + h), path.eval(s - h)
    for k in range(s.size):
        w_fd = _vee(((sp.R[k] - sm.R[k]) / (2 * h)) @ smp.R[k].T)
        assert np.allclose(w_fd, smp.w[:, k], rtol=1e-5, atol=1e-7), f"w 自检失败 k={k}"
    dw_fd = (sp.w - sm.w) / (2 * h)
    ddw_fd = (sp.dw - sm.dw) / (2 * h)
    dp_fd = (sp.p - sm.p) / (2 * h)
    ddp_fd = (sp.dp - sm.dp) / (2 * h)
    dddp_fd = (sp.ddp - sm.ddp) / (2 * h)
    assert np.allclose(dw_fd, smp.dw, rtol=1e-5, atol=1e-6), "dw 自检失败"
    assert np.allclose(ddw_fd, smp.ddw, rtol=1e-5, atol=1e-4), "ddw 自检失败"
    assert np.allclose(dp_fd, smp.dp, rtol=1e-5, atol=1e-7), "dp 自检失败"
    assert np.allclose(ddp_fd, smp.ddp, rtol=1e-5, atol=1e-5), "ddp 自检失败"
    assert np.allclose(dddp_fd, smp.dddp, rtol=1e-5, atol=1e-3), "dddp 自检失败"


def _assert_chain_rule(pd_fine) -> None:
    """q 序列有限差分交叉验证链式法则输出（内部站点）。"""
    s, q = pd_fine.s_grid, pd_fine.q
    inner = slice(3, -3)
    dq_fd = _fd(q, s)
    err1 = np.max(np.abs(dq_fd[:, inner] - pd_fine.dq[:, inner]))
    scale1 = 1.0 + np.max(np.abs(pd_fine.dq))
    assert err1 < 1e-4 * scale1, f"dq 交叉验证失败：err={err1:.3e}, scale={scale1:.3g}"

    ddq_fd = _fd(pd_fine.dq, s)
    err2 = np.max(np.abs(ddq_fd[:, inner] - pd_fine.ddq[:, inner]))
    scale2 = 1.0 + np.max(np.abs(pd_fine.ddq))
    assert err2 < 1e-4 * scale2, f"ddq 交叉验证失败：err={err2:.3e}, scale={scale2:.3g}"

    dddq_fd = _fd(pd_fine.ddq, s)
    err3 = np.max(np.abs(dddq_fd[:, inner] - pd_fine.dddq[:, inner]))
    scale3 = 1.0 + np.max(np.abs(pd_fine.dddq))
    # 若三阶式漏掉 J'q'' 的系数 2，此处误差为 O(10%)·scale，远超阈值
    assert err3 < 1e-3 * scale3, f"dddq 交叉验证失败：err={err3:.3e}, scale={scale3:.3g}"


def _assert_jacobian_identity(path, pd_fine) -> None:
    """独立恒等式：网格差分 J'（不经方向差分代码）满足 J'q' + Jq'' ≈ r₂。"""
    s = pd_fine.s_grid
    smp = path.eval(s)
    r2 = np.vstack([smp.ddp, smp.dw])
    Js = np.stack([_KIN.jacobian(pd_fine.q[:, k]) for k in range(s.size)])  # (N,6,6)
    dJ_grid = np.gradient(Js, s, axis=0)
    res = 0.0
    for k in range(3, s.size - 3):
        res = max(res, float(np.max(np.abs(
            dJ_grid[k] @ pd_fine.dq[:, k] + Js[k] @ pd_fine.ddq[:, k] - r2[:, k]
        ))))
    assert res < 1e-3, f"J'q'+Jq''=r₂ 恒等式残差 {res:.3e} 过大"


def _assert_adaptive_sampling(path) -> None:
    """网格含段边界；相邻站点间弦高误差 ≤ ~2×eps_pos。"""
    opts = SampleOptions(eps_pos=1e-4, eps_ori=1e-3)
    g = adaptive_sample(path, opts)
    assert g[0] == 0.0 and g[-1] == path.s_total
    for b in np.asarray(path.s_breaks):
        assert np.min(np.abs(g - b)) < 1e-8 * path.s_total, f"网格缺少边界 {b}"
    assert np.all(np.diff(g) > 0)
    mids = 0.5 * (g[:-1] + g[1:])
    p_mid = path.eval(mids).p
    p_lo = path.eval(g[:-1]).p
    p_hi = path.eval(g[1:]).p
    sag = np.linalg.norm(p_mid - 0.5 * (p_lo + p_hi), axis=0)
    assert np.max(sag) < 2.0 * opts.eps_pos, f"弦高误差 {np.max(sag):.3e} 超容差"


def _assert_solver_constraints(pd, limits, profile, hist) -> None:
    """端到端：SPLP 收敛单调 + 网格点约束满足（速度/加速度/TCP 速度模）。"""
    tf = hist.t_final
    assert all(x >= y - 1e-6 for x, y in zip(tf, tf[1:])), f"t_final 非单调: {tf}"
    assert np.isfinite(tf[-1]) and tf[-1] > 0

    a, b = profile.a, profile.b
    sa = np.sqrt(np.maximum(a, 0.0))
    vmax, amax, _ = limits.axis_arrays(pd.n_axis)
    qd = pd.dq * sa[None, :]
    qdd = pd.ddq * a[None, :] + pd.dq * b[None, :]
    assert np.all(np.abs(qd) <= vmax[:, None] * 1.02 + 1e-9), "关节速度越界"
    assert np.all(np.abs(qdd) <= amax[:, None] * 1.05 + 1e-9), "关节加速度越界"
    assert np.all(pd.cv * sa <= limits.v_tcp_max * 1.02 + 1e-9), "TCP 位置速度模越界"
    assert np.all(pd.cw * sa <= limits.w_tcp_max * 1.02 + 1e-9), "TCP 姿态角速度模越界"


def test_lowering():
    home = _KIN.fk(_Q_HOME)
    assert min_singular_ratio(_KIN.jacobian(_Q_HOME)) > 1e-4, "家位姿不应奇异"

    # ── 1. 夹具自检 ─────────────────────────────────────────────────────
    rich = _RichPath(p0=home.position.copy(), R0=home.rotation.copy())
    _assert_fixture_selfcheck(rich)
    pure = _PureRotPath(p0=home.position.copy(), R0=home.rotation.copy())
    _assert_fixture_selfcheck(pure)

    # ── 2/3/4. 细网格降维 + FD 交叉验证 + 独立恒等式 ────────────────────
    s_fine = np.linspace(0.0, rich.s_total, 1601)
    pd_fine = lower_cartesian(rich, _KIN, q_seed=_Q_HOME, s_grid=s_fine)
    assert not np.any(pd_fine.singular), "测试路径不应过奇异点"
    smp = rich.eval(s_fine)
    for k in (0, 400, 800, 1200, 1600):       # IK 回代抽查
        fk = _KIN.fk(pd_fine.q[:, k])
        assert np.linalg.norm(fk.position - smp.p[:, k]) < 1e-8
        assert np.linalg.norm(fk.rotation - smp.R[k]) < 1e-8
    assert np.max(np.abs(np.diff(pd_fine.q, axis=1))) < 0.05, "细网格下 q 应准连续"
    _assert_chain_rule(pd_fine)
    _assert_jacobian_identity(rich, pd_fine)

    # ── 5. 自适应采样性质 ───────────────────────────────────────────────
    _assert_adaptive_sampling(rich)

    # ── 6. 纯姿态调整退化安全 ───────────────────────────────────────────
    pd_rot = lower_cartesian(pure, _KIN, q_seed=_Q_HOME)
    assert np.max(pd_rot.cv) < 1e-9, "纯姿态路径 cv 应为 0"
    assert np.min(pd_rot.cw) > 0.0
    res = np.stack([_KIN.jacobian(pd_rot.q[:, k]) @ pd_rot.dq[:, k]
                    for k in range(pd_rot.n_grid)], axis=1)
    assert np.allclose(res[:3], 0.0, atol=1e-9), "纯姿态路径 TCP 线速度应为 0"

    # ── 奇异处理单元检查 ────────────────────────────────────────────────
    J_sing = np.diag([1.0, 1.0, 1.0, 1.0, 1.0, 1e-9])
    x = damped_inverse_solve(J_sing, np.ones(6), lam=0.05)
    assert np.all(np.isfinite(x)) and np.max(np.abs(x)) < 1e3, "DLS 应有界"
    assert min_singular_ratio(J_sing) < 1e-6

    # ── 7. 端到端：自适应网格 → Topp3Data → SPLP ───────────────────────
    limits = load_robot_limits(v_tcp_max=0.5, w_tcp_max=2.5)
    pd = lower_cartesian(rich, _KIN, q_seed=_Q_HOME,
                         sample_opts=SampleOptions(eps_pos=1e-4, eps_ori=1e-3))
    data = limits.to_topp3_data(pd.s_grid, pd.dq, pd.ddq, pd.dddq,
                                tcp_geom=pd.tcp_geom())
    profile, hist = solve_splp(data, SolveOptions(n_iter=4))
    _assert_solver_constraints(pd, limits, profile, hist)

    # ── lower_joint 快路径（解析关节正弦路径）────────────────────────────
    @dataclass
    class _JPath:
        s_total: float = 1.0

        def eval_joint(self, s):
            th = 0.4 * np.sin(np.pi * s[None, :] + np.linspace(0, 1, 6)[:, None])
            q0 = _Q_HOME[:, None] + th
            q1 = 0.4 * np.pi * np.cos(np.pi * s[None, :] + np.linspace(0, 1, 6)[:, None])
            q2 = -np.pi ** 2 * th
            q3 = -np.pi ** 2 * q1
            return q0, q1, q2, q3

    pd_j = lower_joint(_JPath(), _KIN)
    assert pd_j.q.shape[0] == 6 and pd_j.n_grid >= 2
    assert np.all(np.isfinite(pd_j.cv)) and np.max(pd_j.cv) > 0.0

    # ── 出图（可选）────────────────────────────────────────────────────
    _plot(pd, pd_fine, limits, profile, hist)
    print("PASS test_lowering")


def _plot(pd, pd_fine, limits, profile, hist):
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
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle("lowering 降维层：解析笛卡尔路径 → UR5 关节路径 → SPLP（M2）", fontsize=13)

    s = pd.s_grid
    for i in range(6):
        ax[0, 0].plot(pd_fine.s_grid, pd_fine.q[i], label=f"关节 {i}", lw=1.1)
    ax[0, 0].set_title("① 关节角 q(s)（细网格）"); ax[0, 0].legend(fontsize=7, ncol=2)
    for i in range(6):
        ax[0, 1].plot(pd_fine.s_grid, pd_fine.dq[i], lw=1.1)
    ax[0, 1].set_title("② 关节路径导数 q'(s)")
    ax[0, 2].plot(s[:-1], np.diff(s), ".-", ms=3)
    ax[0, 2].set_title("③ 自适应步长 Δs(s)"); ax[0, 2].set_ylim(bottom=0)
    ax[1, 0].plot(s, pd.cv, label=r"$c_v=\|p'\|$"); ax[1, 0].plot(s, pd.cw, label=r"$c_w=\|\hat{\omega}\|$")
    ax[1, 0].set_title("④ TCP 速度模系数"); ax[1, 0].legend(fontsize=8)
    sa = np.sqrt(np.maximum(profile.a, 0.0))
    ax[1, 1].plot(s, sa, label=r"$\dot s=\sqrt{a}$")
    ax[1, 1].plot(s, np.where(pd.cv > 1e-12, limits.v_tcp_max / np.maximum(pd.cv, 1e-12), np.nan),
                  "--", color="0.5", label="TCP 速度模上界折算")
    ax[1, 1].set_title("⑤ 速度剖面与 TCP 折算上界"); ax[1, 1].legend(fontsize=8)
    ax[1, 2].plot(np.arange(1, len(hist.t_final) + 1), hist.t_final, "o-", color="C3")
    ax[1, 2].set_title(f"⑥ SPLP 收敛（t_f={hist.t_final[-1]:.3f}s）")
    for a_ in ax.flat:
        a_.grid(alpha=0.25)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(out_dir, "lowering_test.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    assert os.path.exists(out) and os.path.getsize(out) > 0


if __name__ == "__main__":
    test_lowering()
