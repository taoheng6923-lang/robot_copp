"""LinearMove 直线指令（framework §5.2 / 设计 §3.2）。

位置线性插值 + 姿态 SLERP（世界系角速度密度恒定），对 s 全解析：
    p(s) = p0 + Δp·(s/L)，p' = Δp/L，p'' = p''' = 0
    R(s) = R0·exp((s/L)·[Θ]×)，ω̂ = R0·Θ/L 恒定，ω̂' = ω̂'' = 0
其中 Θ = log(R0ᵀR1)（最短路径）。

段长 L = max(‖Δp‖, rot_scale·θ)：纯姿态调整（Δp=0）时用转角 θ 乘特征长度
rot_scale 参数化，避免 s 区间退化（设计 §6.1 的"纯姿态退化安全"）。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from robot import Pose

from ..types import CartesianSamples
from ..errors import ZeroLengthCommandError
from .base import Section, rotvec_between, slerp_frames


@dataclass
class _LinePath:
    """直线 + SLERP 的 CartesianPath 实现。"""

    p0: np.ndarray
    dp_unit: np.ndarray          # Δp/L，(3,)
    R0: np.ndarray
    rotvec: np.ndarray           # Θ = log(R0ᵀR1)，(3,)
    s_total: float

    def __post_init__(self):
        self.s_breaks = np.array([0.0, self.s_total])
        self._w = (self.R0 @ self.rotvec) / self.s_total   # ω̂ 恒定（世界系）

    def eval(self, s: np.ndarray) -> CartesianSamples:
        s = np.atleast_1d(np.asarray(s, dtype=float))
        N = s.size
        zero3 = np.zeros((3, N))
        return CartesianSamples(
            p=self.p0[:, None] + self.dp_unit[:, None] * s[None, :],
            dp=np.repeat(self.dp_unit[:, None], N, axis=1),
            ddp=zero3, dddp=zero3.copy(),
            R=slerp_frames(self.R0, self.rotvec, s / self.s_total),
            w=np.repeat(self._w[:, None], N, axis=1),
            dw=zero3.copy(), ddw=zero3.copy(),
        )


@dataclass
class LinearMoveCommand:
    """直线指令：pose_start → pose_end（笛卡尔直线 + SLERP）。

    rot_scale : 纯姿态/姿态主导时的角度→弧长折算 [m/rad]（决定 s 参数尺度，
                不影响几何形状，只影响该段网格密度的量纲基准）。
    """

    pose_start: Pose
    pose_end: Pose
    rot_scale: float = 0.1

    def to_section(self) -> Section:
        p0 = np.asarray(self.pose_start.position, dtype=float)
        p1 = np.asarray(self.pose_end.position, dtype=float)
        R0 = np.asarray(self.pose_start.rotation, dtype=float)
        R1 = np.asarray(self.pose_end.rotation, dtype=float)
        dp = p1 - p0
        rotvec = rotvec_between(R0, R1)
        L = max(float(np.linalg.norm(dp)), self.rot_scale * float(np.linalg.norm(rotvec)))
        if L < 1e-12:
            raise ZeroLengthCommandError("LinearMove 起止位姿重合（位移与转角均为零）")
        return Section(
            path=_LinePath(p0=p0, dp_unit=dp / L, R0=R0, rotvec=rotvec, s_total=L),
            native_space="cartesian",
            pose_start=self.pose_start, pose_end=self.pose_end,
        )
