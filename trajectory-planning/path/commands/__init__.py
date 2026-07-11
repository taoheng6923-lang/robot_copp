"""指令层：三类运动指令 → Section → 逐段 PathDerivatives（framework §5.2 / 设计 §3，M2 已实现）。

- base:          Section 数据结构 + MotionCommand 协议 + SLERP 工具
- joint_move:    JointMoveCommand（线性关节几何，见其 docstring 与设计 §3.1 的差异说明）
- linear_move:   LinearMoveCommand（笛卡尔直线 + SLERP，纯姿态退化安全）
- circular_move: CircularMoveCommand（三点定圆 / 圆心+法向，退化显式报错）
- assemble:      build_sections（翻译）+ lower_sections（G0 衔接校验 + seed 链逐段降维）

M2 无 blending：段间角点以停顿衔接（每段独立 rest-to-rest 求解）；G2 过渡属 M3。
"""

from .base import Section, MotionCommand
from .joint_move import JointMoveCommand
from .linear_move import LinearMoveCommand
from .circular_move import CircularMoveCommand
from .assemble import build_sections, lower_sections

__all__ = [
    "Section", "MotionCommand",
    "JointMoveCommand", "LinearMoveCommand", "CircularMoveCommand",
    "build_sections", "lower_sections",
]
