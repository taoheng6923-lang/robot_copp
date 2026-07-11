"""运动学 / 动力学后端抽象（framework §5.1）。

把机器人本体（FK/IK/Jacobian/逆动力学）经协议解耦，供上层注入：
  - blending  用 FK 桥接混合指令；
  - lowering  用 IK + Jacobian 链把笛卡尔路径降维为关节路径导数；
  - constraints 用逆动力学生成力矩行。

`ur5.UR5Kinematics` 用真实 UR5 DH 参数落地 KinematicsModel 协议（fk/jacobian
解析解，ik 解析逆解 + DLS 兜底）；DynamicsModel 仍是 `ur5.UR5RobotModel.torque_coeffs`
的对角近似 stand-in，真实 RNE 待 M2+ 落地。其余机型可仿 `ur5.py` 新增
DhPoeKinematics / RtbKinematics 适配器。

本模块只定义接口与 `Pose` 数据类型，不含规划算法（纯适配层）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np


@dataclass
class Pose:
    """位姿：位置 + 旋转矩阵（相对基座）。fk/ik 协议的交换格式。"""

    position: np.ndarray   # (3,)
    rotation: np.ndarray   # (3,3)


@runtime_checkable
class KinematicsModel(Protocol):
    """FK / IK / Jacobian（M2：DhPoeKinematics / RtbKinematics）。"""

    def fk(self, q: np.ndarray):
        """正运动学：关节角 q → TCP 位姿。"""
        ...

    def jacobian(self, q: np.ndarray) -> np.ndarray:
        """几何 Jacobian，(6, n)。"""
        ...

    def jacobian_derivative(self, q: np.ndarray, dq: np.ndarray) -> np.ndarray:
        """沿路径方向的 dJ/ds（有限差分），(6, n)。"""
        ...

    def ik(self, pose, seed: np.ndarray) -> np.ndarray:
        """逆运动学：位姿 → 离 seed 最近的关节解（连续解选择）。"""
        ...


@runtime_checkable
class DynamicsModel(Protocol):
    """逆动力学 RNE（M4：力矩约束 / 热能目标）。"""

    def inverse_dynamics(
        self, q: np.ndarray, dq: np.ndarray, ddq: np.ndarray
    ) -> np.ndarray:
        """RNE：(q, q̇, q̈) → 关节力矩 (n,)。"""
        ...
