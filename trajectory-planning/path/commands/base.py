"""指令层公共类型（framework §5.2 base / 设计 §3.4）。

一条运动指令翻译为一个 `Section`：其 `path` 是几何路径（笛卡尔或关节原生），
加上端点元数据（供段间衔接校验与 IK seed 传递）。M2 无 blending：相邻段只做
G0 精确衔接校验，每段独立 rest-to-rest 求解（设计 §4 的 G2 过渡属 M3）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Protocol, runtime_checkable

import numpy as np

from robot import Pose

from ..types import CartesianPath, JointSpacePath


@dataclass
class Section:
    """一条指令的翻译结果：几何路径 + 端点元数据。

    path         : CartesianPath（Linear/Circular）或 JointSpacePath（JointMove）
    native_space : "cartesian" / "joint"（决定降维走 IK 链还是快路径）
    q_start/q_end       : 关节段端点关节角；笛卡尔段为 None
    pose_start/pose_end : 笛卡尔段端点位姿；关节段为 None（需要时经 FK 求）
    """

    path: "CartesianPath | JointSpacePath"
    native_space: Literal["cartesian", "joint"]
    q_start: np.ndarray | None = None
    q_end: np.ndarray | None = None
    pose_start: Pose | None = None
    pose_end: Pose | None = None


@runtime_checkable
class MotionCommand(Protocol):
    """运动指令协议：翻译为 Section（纯几何，不做 IK）。"""

    def to_section(self) -> Section: ...


def rotvec_between(R0: np.ndarray, R1: np.ndarray) -> np.ndarray:
    """R0→R1 的旋转向量 log(R0ᵀR1)（角度 ≤ π，最短路径；θ=π 时轴向有二义性）。"""
    from scipy.spatial.transform import Rotation

    return Rotation.from_matrix(R0.T @ R1).as_rotvec()


def slerp_frames(R0: np.ndarray, rotvec: np.ndarray, u: np.ndarray) -> np.ndarray:
    """R(u) = R0·exp(u·[rotvec]×)，u (N,) → (N,3,3)。

    世界系角速度密度恒定：Ṙ Rᵀ = [R0·rotvec]×·u̇（rotvec 是 exp 的旋转轴，
    被自身旋转固定，故 R(u)·rotvec = R0·rotvec 与 u 无关）。
    """
    from scipy.spatial.transform import Rotation

    Rrel = Rotation.from_rotvec(np.outer(np.atleast_1d(u), rotvec)).as_matrix()
    return R0[None, :, :] @ Rrel
