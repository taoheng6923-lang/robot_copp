"""降维层：笛卡尔/关节路径 → PathDerivatives（framework §5.4 / 设计 §5，M2 已实现）。

- sampling:    曲率驱动自适应离散化（网格含 s_breaks）
- ik:          连续解 IK（seed 链 + FK 回代/跳变校验）
- derivatives: Jacobian 链式法则求 q',q'',q'''（J'/J'' 方向差分）+ lower_* 驱动
- singularity: σ_min/σ_max 奇异检测 + 阻尼最小二乘逆
"""

from .sampling import SampleOptions, adaptive_sample, uniform_sample
from .ik import IkOptions, solve_ik_sequence
from .derivatives import (
    DerivativeOptions, joint_derivatives, lower_cartesian, lower_joint,
)
from .singularity import min_singular_ratio, damped_inverse_solve

__all__ = [
    "SampleOptions", "adaptive_sample", "uniform_sample",
    "IkOptions", "solve_ik_sequence",
    "DerivativeOptions", "joint_derivatives", "lower_cartesian", "lower_joint",
    "min_singular_ratio", "damped_inverse_solve",
]
