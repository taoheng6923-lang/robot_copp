"""合成机器人模型（M1 stand-in）。

M1 无真实 FK / 逆动力学，用解析式集中提供 SPLP 求解与可视化所需的本体量：
  - joint_path(s)     合成关节路径 q(s) 及导数（stand-in：实际由 lowering 的 IK 链给出）
  - tcp_geometry(s)   TCP 位置导数 p'(s)、角速度方向 ω/ṡ（stand-in：FK + Jacobian 链）
  - tcp_coeffs(s)     TCP 速度模系数 {cv=‖p'‖, cw=‖ω/ṡ‖}（供 to_topp3_data）
  - torque_coeffs(...) 力矩系数 {n_tor, m_tor, g_tor}（stand-in：逆动力学 RNE）

把这些“机器人运动学 / 动力学计算”单独管理在此，测试 / 上层只调用、不内联。
M2 起由 kinematics.base 的 KinematicsModel + DynamicsModel 真实链路替换本模块。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SyntheticRobotModel:
    """解析 TCP 螺旋路径 + 对角常惯性动力学的合成本体。

    n_axis        关节数
    seed          合成关节路径的随机种子（决定幅值/频率/相位）
    inertia       对角常惯性 M_i（力矩 a、b 系数）
    gravity       重力幅值 G_i（力矩常数项）
    tcp_radius    TCP 位置螺旋半径 R
    tcp_pitch     TCP 位置螺旋螺距分量 h（z 向 p'）
    """

    n_axis: int = 3
    seed: int = 3
    inertia: float = 0.4
    gravity: float = 1.5
    tcp_radius: float = 0.06
    tcp_pitch: float = 0.05

    # ── 运动学：关节路径（lowering/IK 的 stand-in） ──────────────────────
    def joint_path(self, s: np.ndarray):
        """返回 (q0, q1, q2, q3)，各 (n, N)：q(s) 及其 1/2/3 阶导。"""
        rng = np.random.default_rng(self.seed)
        A = rng.uniform(0.3, 0.5, self.n_axis)
        w = rng.uniform(1.0, 2.0, self.n_axis)
        phi = rng.uniform(0.0, 2.0 * np.pi, self.n_axis)
        th = w[:, None] * s[None, :] + phi[:, None]
        q0 = A[:, None] * np.sin(th)
        q1 = A[:, None] * w[:, None] * np.cos(th)
        q2 = -A[:, None] * w[:, None] ** 2 * np.sin(th)
        q3 = -A[:, None] * w[:, None] ** 3 * np.cos(th)
        return q0, q1, q2, q3

    # ── 运动学：TCP 几何（FK + Jacobian 的 stand-in） ────────────────────
    def tcp_geometry(self, s: np.ndarray) -> dict:
        """返回 {dp, wdir}，各 (3, N)：位置速度模 = ‖dp‖·√a、姿态角速度模 = ‖wdir‖·√a。"""
        w = 2.0 * np.pi
        dp = np.vstack([
            -w * self.tcp_radius * np.sin(w * s),
            w * self.tcp_radius * np.cos(w * s),
            self.tcp_pitch * np.ones_like(s),
        ])
        k = 3.0 * np.pi
        wdir = np.vstack([
            0.25 * np.cos(k * s),
            0.20 + 0.10 * s,
            0.20 * np.sin(k * s),
        ])
        return {"dp": dp, "wdir": wdir}

    def tcp_coeffs(self, s: np.ndarray) -> dict:
        """返回 {cv, cw}，各 (N,)：TCP 速度模系数（供 to_topp3_data(tcp_geom=...)）。"""
        g = self.tcp_geometry(s)
        return {
            "cv": np.linalg.norm(g["dp"], axis=0),
            "cw": np.linalg.norm(g["wdir"], axis=0),
        }

    # ── 动力学：力矩系数（逆动力学 RNE 的 stand-in） ─────────────────────
    def torque_coeffs(self, q0: np.ndarray, q1: np.ndarray, q2: np.ndarray) -> dict:
        """τ = n_tor·a + m_tor·b + g_tor 的系数（对角惯性 + 重力，略科氏/摩擦）。

        n_tor = M·q''（a 系数），m_tor = M·q'（b 系数），g_tor = G·sin(q)（常数项）。
        返回 {n_tor, m_tor, g_tor}，各 (n, N)（供 to_topp3_data(torque_coeffs=...)）。
        """
        return {
            "n_tor": self.inertia * q2,
            "m_tor": self.inertia * q1,
            "g_tor": self.gravity * np.sin(q0),
        }
