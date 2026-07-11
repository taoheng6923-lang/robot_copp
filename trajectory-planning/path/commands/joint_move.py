"""JointMove 关节运动指令（framework §5.2 / 设计 §3.1）。

几何采用**线性关节插值** q(s) = q_start + Δq·(s/L)，L = ‖Δq‖₂：
q' = Δq/L 恒定、q'' = q''' = 0。

与设计文档 §3.1 的差异（有意为之）：设计稿建议五次多项式参数化 q(u)（两端
q̇=q̈=0），那是"时间域直接参数化"的经典做法；TOPP 框架下几何与时间律解耦——
时间平滑（rest-to-rest、jerk 受限）由 copp 时间律 s(t) 保证，几何上的五次
多项式反而使端点 q'(s)=0（参数化退化、路径非正则，端点约束系数消失）。
故取线性几何；段间的几何平滑过渡属 M3 blending。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..errors import ZeroLengthCommandError
from .base import Section


@dataclass
class _LinearJointPath:
    """q(s) = q_start + dq_unit·s，s ∈ [0, L]（实现 JointSpacePath 协议）。"""

    q_start: np.ndarray
    dq_unit: np.ndarray        # Δq/L，(n,)
    s_total: float

    def eval_joint(self, s: np.ndarray):
        s = np.atleast_1d(np.asarray(s, dtype=float))
        q = self.q_start[:, None] + self.dq_unit[:, None] * s[None, :]
        dq = np.repeat(self.dq_unit[:, None], s.size, axis=1)
        zeros = np.zeros_like(q)
        return q, dq, zeros, zeros.copy()


@dataclass
class JointMoveCommand:
    """关节运动指令：q_start → q_end（线性关节几何，见模块 docstring）。"""

    q_start: np.ndarray
    q_end: np.ndarray

    def to_section(self) -> Section:
        q0 = np.asarray(self.q_start, dtype=float)
        q1 = np.asarray(self.q_end, dtype=float)
        L = float(np.linalg.norm(q1 - q0))
        if L < 1e-12:
            raise ZeroLengthCommandError("JointMove 起止关节角重合")
        return Section(
            path=_LinearJointPath(q_start=q0, dq_unit=(q1 - q0) / L, s_total=L),
            native_space="joint",
            q_start=q0, q_end=q1,
        )
