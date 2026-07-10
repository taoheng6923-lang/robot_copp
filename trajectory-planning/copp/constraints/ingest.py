"""约束摄入辅助（framework §5.5 / 设计 §6，M4）。

- TCP 速度模长 → a 的逐点上界（线性，折进 ā）；
- 关节力矩 → LP 不等式行（2 阶，对 (a,b) 精确线性）。

TCP 与力矩系数在 M1 由合成模型给出；实际管线中 TCP 系数来自 Jacobian、
力矩系数来自逆动力学（M2 的 KinematicsModel / DynamicsModel）。
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp

from ..types import TcpConstraint, TorqueConstraint


def tcp_a_upper_bound(
    tcp: TcpConstraint, n_grid: int,
    position: bool = True, orientation: bool = True,
) -> np.ndarray:
    """TCP 速度模长 → a 的逐点上界 (N,)：min(v_max²/cv², w_max²/cw²)。

    position/orientation 分别控制是否计入位置速度模、姿态角速度模（约束开关）。
    系数为 0 处（该方向不产生约束）用 +inf 占位。
    """
    bound = np.full(n_grid, np.inf)
    with np.errstate(divide="ignore", invalid="ignore"):
        if position:
            b_pos = np.where(tcp.cv > 0, tcp.v_max**2 / tcp.cv**2, np.inf)
            bound = np.minimum(bound, b_pos)
        if orientation:
            b_ori = np.where(tcp.cw > 0, tcp.w_max**2 / tcp.cw**2, np.inf)
            bound = np.minimum(bound, b_ori)
    return bound


def torque_constraints(
    torque: TorqueConstraint, a: cp.Variable, b: cp.Variable
) -> list[cp.Constraint]:
    """关节力矩约束的 cvxpy 行：τ_min ≤ n_tor·a + m_tor·b + g_tor ≤ τ_max（逐轴）。"""
    cons: list[cp.Constraint] = []
    for i in range(torque.n_tor.shape[0]):
        tau = (
            cp.multiply(torque.n_tor[i], a)
            + cp.multiply(torque.m_tor[i], b)
            + torque.g_tor[i]
        )
        cons += [tau <= torque.tau_max[i], tau >= torque.tau_min[i]]
    return cons
