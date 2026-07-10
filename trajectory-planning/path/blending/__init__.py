"""最优 Hermite blending 层（M3，尚未实现）。

相邻 PoseSegment 交接处在笛卡尔空间做 G2 五次 Hermite 过渡，装配出
C²/分段C³ 连续的 BlendedPath。见 docs/python_framework.md §5.3、
docs/robot_copp_design.md §4。
"""
