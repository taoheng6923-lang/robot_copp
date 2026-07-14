"""三维模型仿真环境：UR5 运动链的 3D 可视化与轨迹动画回放。

定位与依赖方向
--------------
只依赖本包的运动学（`UR5Kinematics`，标准 DH）与 matplotlib，**不依赖**
`trajectory-planning/`——保持 `robot/`（机器人本体）不反向依赖规划层。规划侧
只需把求得的关节轨迹 `q(t)` 递进来即可回放（见 `animate_joint_motion`）。

渲染方式
--------
"骨架"式（非网格模型，仓库内无 URDF/STL 资产，也不引入 VTK 等重依赖）：

  · 连杆     —— 相邻 DH 坐标系原点连成的线段（粗线 + 关节圆点）
  · 地面投影 —— 运动链在 z=地面 上的灰色投影（"影子"，显著提升 3D 纵深可读性）
  · TCP 轨迹 —— 全程路径（浅色）+ 已走过部分的拖尾（高亮）
  · 工具坐标系 —— 末端 T6 的 RGB 三轴（x/y/z），体现姿态而不只是位置
  · 读数     —— 当前时刻 t 与六个关节角（度）

几何来源：`UR5Kinematics._frames(q)` 返回 T[0..6]（T[0]=基座单位阵），其原点
即各关节中心，末端 T[6] 即 TCP 位姿——与 `fk`/`jacobian` 同一套 DH，故动画所
显示的位形与求解器所用的运动学**严格一致**（同源，不是另画一套近似模型）。
"""

from __future__ import annotations

import numpy as np

from .ur5 import N_AXIS, UR5Kinematics

# 配色（与 copp/viz.py 的 matplotlib 风格一致，浅底深线）
_C_LINK = "#2c3e50"    # 连杆
_C_JOINT = "#e67e22"   # 关节
_C_BASE = "#7f8c8d"    # 基座
_C_TRACE = "#c0c0c0"   # TCP 全程路径（未走过）
_C_TRAIL = "#e74c3c"   # TCP 已走过的拖尾
_C_SHADOW = "#bdc3c7"  # 地面投影
_C_AXES = ("#e74c3c", "#27ae60", "#2980b9")  # 工具坐标系 x/y/z


def link_origins(kin: UR5Kinematics, q: np.ndarray) -> np.ndarray:
    """单个位形 q (6,) → 运动链各坐标系原点 (7, 3)：基座 + 6 个关节（末点即 TCP）。"""
    # 同包内部访问：_frames 与 fk/jacobian 同源，保证动画位形与求解器运动学一致
    return np.array([T[:3, 3] for T in kin._frames(np.asarray(q, dtype=float))])


def chain_positions(kin: UR5Kinematics, q_t: np.ndarray) -> np.ndarray:
    """关节轨迹 q_t (6, N) → 逐帧运动链原点 (N, 7, 3)。"""
    q_t = np.asarray(q_t, dtype=float)
    assert q_t.ndim == 2 and q_t.shape[0] == N_AXIS, f"q_t 应为 ({N_AXIS}, N)，得到 {q_t.shape}"
    return np.stack([link_origins(kin, q_t[:, k]) for k in range(q_t.shape[1])])


def _tool_frame(kin: UR5Kinematics, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """末端工具坐标系：返回 (原点 (3,), 旋转 (3,3)——列为 x/y/z 轴方向)。"""
    T6 = kin._frames(np.asarray(q, dtype=float))[-1]
    return T6[:3, 3].copy(), T6[:3, :3].copy()


def _set_cjk_font() -> None:
    """让 matplotlib 正常显示中文（Windows 常见 CJK 字体）——与 copp/viz.py 同一套。"""
    import matplotlib

    matplotlib.rcParams["font.sans-serif"] = [
        "Microsoft YaHei", "SimHei", "DengXian", "SimSun",
    ]
    matplotlib.rcParams["axes.unicode_minus"] = False


def _setup_axes(ax, pts: np.ndarray, margin: float = 0.12):
    """按整段轨迹扫过的空间设**等比例**立方视界。返回地面高度 z_floor（供投影用）。

    三轴用同一半边长 + set_box_aspect((1,1,1))，否则机器人会被各轴独立缩放而扭曲变形。
    地面参照直接用 matplotlib 3D 自带的窗格网格（z 下界处），不再自绘。
    """
    lo, hi = pts.reshape(-1, 3).min(axis=0), pts.reshape(-1, 3).max(axis=0)
    center = 0.5 * (lo + hi)
    half = 0.5 * float(np.max(hi - lo)) + margin  # 立方体半边长 → 等比例，不失真
    z_floor = float(min(0.0, lo[2]) - 0.02)
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(z_floor, center[2] + half)
    ax.set_box_aspect((1.0, 1.0, 1.0))

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    return z_floor


def animate_joint_motion(
    q_t: np.ndarray,
    t: np.ndarray | None = None,
    *,
    kin: UR5Kinematics | None = None,
    save_path: str | None = None,
    show: bool = False,
    fps: int = 25,
    dpi: int = 90,
    title: str = "UR5 三维仿真回放",
    elev: float = 22.0,
    azim: float = 40.0,
):
    """回放关节轨迹 q(t)：3D 骨架动画。

    q_t       : (6, N) 关节角轨迹（**按等时间间隔采样**，逐帧播放）
    t         : (N,) 各帧时刻 [s]；None 则只显示帧序号
    kin       : UR5Kinematics；None 则新建
    save_path : 存 GIF 路径（.gif，用 PillowWriter）；None 不存
    show      : True 则弹交互窗口播放（可拖拽旋转视角）
    fps       : 播放/导出帧率
    dpi       : 导出 GIF 的分辨率（3D 逐帧重绘较慢，导出耗时 ∝ 帧数×dpi²）

    返回 (fig, anim)。**注意**：`anim` 必须被持有引用，否则会被 GC 回收导致动画停摆。
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    _set_cjk_font()
    kin = kin or UR5Kinematics()
    q_t = np.asarray(q_t, dtype=float)
    pts = chain_positions(kin, q_t)          # (N, 7, 3)
    tcp = pts[:, -1, :]                      # (N, 3)
    n_frames = pts.shape[0]

    fig = plt.figure(figsize=(8.0, 7.0))
    ax = fig.add_subplot(111, projection="3d")
    z_floor = _setup_axes(ax, pts)
    ax.view_init(elev=elev, azim=azim)
    ax.set_title(title)

    # 静态层：TCP 全程路径（浅色，作为"将要走的路"的参照）
    ax.plot(tcp[:, 0], tcp[:, 1], tcp[:, 2], color=_C_TRACE, lw=1.0, ls="--", zorder=1)

    # 动态层
    (shadow,) = ax.plot([], [], [], color=_C_SHADOW, lw=3.0, solid_capstyle="round", zorder=1)
    (link,) = ax.plot([], [], [], color=_C_LINK, lw=5.0, solid_capstyle="round", zorder=3)
    (joints,) = ax.plot([], [], [], "o", color=_C_JOINT, ms=7.0, zorder=4)
    (trail,) = ax.plot([], [], [], color=_C_TRAIL, lw=1.8, zorder=2)
    (base,) = ax.plot([], [], [], "s", color=_C_BASE, ms=11.0, zorder=3)
    triad = [ax.plot([], [], [], color=c, lw=2.0, zorder=5)[0] for c in _C_AXES]
    readout = ax.text2D(0.02, 0.96, "", transform=ax.transAxes, family="monospace",
                        fontsize=9, va="top")

    axis_len = 0.09  # 工具坐标系三轴长度 [m]

    def update(k: int):
        P = pts[k]                              # (7,3)
        link.set_data(P[:, 0], P[:, 1]); link.set_3d_properties(P[:, 2])
        joints.set_data(P[1:, 0], P[1:, 1]); joints.set_3d_properties(P[1:, 2])
        base.set_data(P[:1, 0], P[:1, 1]); base.set_3d_properties(P[:1, 2])
        # 地面投影（把 z 全压到地面）
        shadow.set_data(P[:, 0], P[:, 1])
        shadow.set_3d_properties(np.full(P.shape[0], z_floor))
        # TCP 拖尾（已走过的部分）
        trail.set_data(tcp[: k + 1, 0], tcp[: k + 1, 1])
        trail.set_3d_properties(tcp[: k + 1, 2])
        # 末端工具坐标系三轴
        o, R = _tool_frame(kin, q_t[:, k])
        for j, ln in enumerate(triad):
            e = o + axis_len * R[:, j]
            ln.set_data([o[0], e[0]], [o[1], e[1]])
            ln.set_3d_properties([o[2], e[2]])
        # 读数
        head = f"frame {k + 1}/{n_frames}" if t is None else f"t = {t[k]:6.3f} s / {t[-1]:.3f} s"
        deg = np.rad2deg(q_t[:, k])
        readout.set_text(head + "\nq = [" + " ".join(f"{v:7.1f}" for v in deg) + "]  deg")
        return (link, joints, base, shadow, trail, readout, *triad)

    # blit=False：3D 轴的 blit 在多数后端不可靠（artist 的 z 排序需重算）
    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=1000.0 / fps, blit=False, repeat=True)

    if save_path:
        from matplotlib.animation import PillowWriter  # pillow 随 matplotlib 装，无需额外依赖
        anim.save(save_path, writer=PillowWriter(fps=fps), dpi=dpi)
    if show:
        plt.show()
    return fig, anim


def plot_pose(
    q: np.ndarray,
    *,
    kin: UR5Kinematics | None = None,
    save_path: str | None = None,
    show: bool = False,
    title: str = "UR5 位形",
):
    """静态快照：画单个位形 q (6,)。返回 fig。"""
    import matplotlib.pyplot as plt

    _set_cjk_font()
    kin = kin or UR5Kinematics()
    P = link_origins(kin, q)
    fig = plt.figure(figsize=(7.0, 6.5))
    ax = fig.add_subplot(111, projection="3d")
    z_floor = _setup_axes(ax, P[None, ...])
    ax.view_init(elev=22.0, azim=40.0)
    ax.set_title(title)

    ax.plot(P[:, 0], P[:, 1], np.full(P.shape[0], z_floor), color=_C_SHADOW, lw=3.0, zorder=1)
    ax.plot(P[:, 0], P[:, 1], P[:, 2], color=_C_LINK, lw=5.0, solid_capstyle="round", zorder=3)
    ax.plot(P[1:, 0], P[1:, 1], P[1:, 2], "o", color=_C_JOINT, ms=7.0, zorder=4)
    ax.plot(P[:1, 0], P[:1, 1], P[:1, 2], "s", color=_C_BASE, ms=11.0, zorder=3)
    o, R = _tool_frame(kin, q)
    for j, c in enumerate(_C_AXES):
        e = o + 0.09 * R[:, j]
        ax.plot([o[0], e[0]], [o[1], e[1]], [o[2], e[2]], color=c, lw=2.0, zorder=5)

    if save_path:
        fig.savefig(save_path, dpi=130, bbox_inches="tight")
    if show:
        plt.show()
    return fig
