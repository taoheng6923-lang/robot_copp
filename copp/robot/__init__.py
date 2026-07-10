"""机器人本体：运动学 / 动力学后端（framework §5.1）。

集中管理机器人本体的计算（FK/IK/Jacobian/逆动力学），与规划算法（solve/）解耦。

- base:      KinematicsModel / DynamicsModel 协议（M2 目标接口）
- synthetic: SyntheticRobotModel（M1 stand-in：解析 TCP 路径 + 对角惯性动力学）
"""

from .base import KinematicsModel, DynamicsModel
from .synthetic import SyntheticRobotModel

__all__ = ["KinematicsModel", "DynamicsModel", "SyntheticRobotModel"]
