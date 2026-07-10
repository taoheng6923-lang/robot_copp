"""机器人本体：运动学 / 动力学后端（framework §5.1）。

集中管理机器人本体的计算（FK/IK/Jacobian/逆动力学），与规划算法（solve/）解耦。

- base: KinematicsModel / DynamicsModel 协议
- ur5:  UR5Kinematics（真实 DH 运动学，实现 KinematicsModel 协议）+
        UR5RobotModel（真实 TCP 几何 + 对角近似动力学，供 SPLP 测试/可视化）
"""

from .base import KinematicsModel, DynamicsModel
from .ur5 import UR5Kinematics, UR5RobotModel, Pose

__all__ = ["KinematicsModel", "DynamicsModel", "UR5Kinematics", "UR5RobotModel", "Pose"]
