"""机器人本体：运动学 / 动力学后端（framework §5.1）。

集中管理机器人本体的计算（FK/IK/Jacobian/逆动力学），与规划算法（solve/）解耦。

- base:  KinematicsModel / DynamicsModel 协议
- ur5:   UR5Kinematics（真实 DH 运动学，实现 KinematicsModel 协议）+
         UR5RobotModel（真实 TCP 几何 + 对角近似动力学，供 SPLP 测试/可视化）
- sim3d: 三维仿真环境——运动链的 3D 骨架渲染与 q(t) 轨迹动画回放（matplotlib，
         函数内惰性导入，故未装 matplotlib 时本包仍可正常 import）
"""

from .base import KinematicsModel, DynamicsModel, Pose
from .ur5 import UR5Kinematics, UR5RobotModel
from .sim3d import animate_joint_motion, plot_pose, link_origins, chain_positions

__all__ = [
    "KinematicsModel", "DynamicsModel", "Pose", "UR5Kinematics", "UR5RobotModel",
    "animate_joint_motion", "plot_pose", "link_origins", "chain_positions",
]
