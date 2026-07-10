"""M1 核心数据类型（对齐设计文档 §7 / framework §4）。

约定：路径参数 s，网格 N 个点（索引 0..N-1），区间 k=1..N-1。
copp 状态：a=ṡ²，b=s̈，c=b'=s⃛/ṡ（每区间常值，见 paper_notes §4）。

M1 仅承载 3 阶 TOPP（时间最优）所需的**轴向**约束数据（速度/加速度/jerk），
路径导数 q',q'',q''' 由上层（合成路径或 lowering 层）预先算好后填入。
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass
class TcpConstraint:
    """笛卡尔 TCP 速度模长约束（M4，设计 §6 / v0.3）。

    位置速度模 ‖ṗ‖=cv·√a、姿态角速度模 ‖ω‖=cw·√a，均为 a=ṡ² 的线性上界，
    最终折进 a 的逐点上界 ā（见 solve/state.velocity_upper_bound）。

    cv    : (N,)  位置速度系数 ‖p'(s)‖
    cw    : (N,)  姿态角速度系数 ‖J_ω(q)·q'(s)‖
    v_max : 位置速度模上界
    w_max : 姿态角速度模上界
    """

    cv: np.ndarray
    cw: np.ndarray
    v_max: float
    w_max: float


@dataclass
class TorqueConstraint:
    """关节力矩约束（M4，2 阶、对 (a,b) 精确线性；论文 eq.44 / robot6dof §5.2.5）。

    τ = n_tor·a + m_tor·b + g_tor，其中
      n_tor = M·q'' + C(q,q')·q'   （a=ṡ² 系数：惯性曲率 + 科氏）
      m_tor = M·q'                 （b=s̈ 系数）
      g_tor = g(q) + 摩擦           （常数项）
    约束 τ_min ≤ τ ≤ τ_max（逐轴、上下各一行）。系数由逆动力学预计算（M2 的
    DynamicsModel），M1 演示中用合成模型给出。

    n_tor,m_tor,g_tor : (n, N)
    tau_max,tau_min   : (n,)
    """

    n_tor: np.ndarray
    m_tor: np.ndarray
    g_tor: np.ndarray
    tau_max: np.ndarray
    tau_min: np.ndarray


@dataclass
class Topp3Data:
    """3 阶 TOPP 的离散输入（轴向约束版本）。

    Attributes
    ----------
    s_grid : (N,)      路径参数网格（单调递增）
    dq     : (n, N)    q'(s)   —— 关节路径一阶导
    ddq    : (n, N)    q''(s)  —— 二阶导
    dddq   : (n, N)    q'''(s) —— 三阶导
    vmax   : (n,)      轴向速度上界（对称）
    amax   : (n,)      轴向加速度上界（对称）
    jmax   : (n,)      轴向 jerk 上界（对称）
    a_bnd  : (2,)      两端 a=ṡ² 边界 (a_start, a_final)
    b_bnd  : (2,)      两端 b=s̈ 边界 (b_start, b_final)
    """

    s_grid: np.ndarray
    dq: np.ndarray
    ddq: np.ndarray
    dddq: np.ndarray
    vmax: np.ndarray
    amax: np.ndarray
    jmax: np.ndarray
    a_bnd: tuple[float, float] = (0.0, 0.0)
    b_bnd: tuple[float, float] = (0.0, 0.0)
    tcp: "TcpConstraint | None" = None       # M4：TCP 速度模长约束（可选）
    torque: "TorqueConstraint | None" = None  # M4：关节力矩约束（可选）

    @property
    def n_axis(self) -> int:
        return self.dq.shape[0]

    @property
    def n_grid(self) -> int:
        return self.s_grid.size

    def validate(self) -> None:
        N, n = self.n_grid, self.n_axis
        assert self.s_grid.ndim == 1 and N >= 4, "s_grid 至少 4 点"
        assert np.all(np.diff(self.s_grid) > 0), "s_grid 必须严格递增"
        for name, arr in (("dq", self.dq), ("ddq", self.ddq), ("dddq", self.dddq)):
            assert arr.shape == (n, N), f"{name} 形状应为 ({n},{N})"
        for name, arr in (("vmax", self.vmax), ("amax", self.amax), ("jmax", self.jmax)):
            assert arr.shape == (n,) and np.all(arr > 0), f"{name} 应为正的 ({n},)"
        if self.tcp is not None:
            assert self.tcp.cv.shape == (N,) and self.tcp.cw.shape == (N,), "TCP cv/cw 形状应为 (N,)"
            assert self.tcp.v_max > 0 and self.tcp.w_max > 0, "TCP v_max/w_max 应为正"
        if self.torque is not None:
            for nm, arr in (("n_tor", self.torque.n_tor), ("m_tor", self.torque.m_tor), ("g_tor", self.torque.g_tor)):
                assert arr.shape == (n, N), f"torque.{nm} 形状应为 ({n},{N})"
            assert self.torque.tau_max.shape == (n,) and self.torque.tau_min.shape == (n,), "tau_max/min 形状应为 (n,)"
            assert np.all(self.torque.tau_max >= self.torque.tau_min), "需 tau_max ≥ tau_min"


@dataclass
class Profile:
    """求解结果剖面 (a, b, c)（设计 §7 3 阶方法）。

    a : (N,)  ṡ²
    b : (N,)  s̈
    c : (N,)  每区间控制 c=b'（c[0] 置 0；c[k] 对应区间 (k-1,k)）
    """

    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    num_stationary: tuple[int, int] = (0, 0)
