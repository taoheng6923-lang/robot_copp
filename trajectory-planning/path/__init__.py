"""路径构造模块（指令 → 笛卡尔路径 → 关节路径，framework §5.2-5.4 / 设计 §3-5）。

指令序列（关节运动/直线/圆弧）经 lowering（自适应采样 + 连续解 IK + 链式法则
求导）翻译为 copp 求解所需的 PathDerivatives。与 `trajectory-planning/copp/`
（纯数值求解核心）、顶层 `robot/`（运动学动力学本体）相互独立，经
`trajectory-planning/planner/`（M2+）编排串联。

- types.py   CartesianSamples / CartesianPath / JointSpacePath / PathDerivatives
- errors.py  PathError 异常层次
- commands/  指令层（M2 已实现；段间 G0 衔接，G2 blending 属 M3）
- lowering/  降维层（M2 已实现）
- blending/  最优 Hermite blending（M3，未实现）
"""

from .types import CartesianSamples, CartesianPath, JointSpacePath, PathDerivatives
from .errors import (
    PathError, ZeroLengthCommandError, DegenerateArcError,
    UnreachablePoseError, IkJumpError, JunctionMismatchError,
)

__all__ = [
    "CartesianSamples", "CartesianPath", "JointSpacePath", "PathDerivatives",
    "PathError", "ZeroLengthCommandError", "DegenerateArcError",
    "UnreachablePoseError", "IkJumpError", "JunctionMismatchError",
]
