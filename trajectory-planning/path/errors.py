"""路径构造层异常层次（framework §5.10 errors 的 path 子集）。"""

from __future__ import annotations


class PathError(Exception):
    """路径构造 / 降维层错误基类。"""


class ZeroLengthCommandError(PathError):
    """指令起止重合（位移与转角均为零），无法参数化。"""


class DegenerateArcError(PathError):
    """圆弧退化：三点共线 / 起终点重合 / 端点不在给定圆上或圆平面外。"""


class UnreachablePoseError(PathError):
    """IK 解算后 FK 回代与目标位姿不符（超出工作空间或全支路退化）。"""


class IkJumpError(PathError):
    """相邻站点关节解跳变超阈值（解分支切换 / 采样过疏 / 过奇异点）。"""


class JunctionMismatchError(PathError):
    """相邻指令段端点位姿 / 关节角不衔接（M2 无 blending，要求 G0 精确衔接）。"""
