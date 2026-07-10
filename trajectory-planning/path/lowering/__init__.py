"""降维到关节空间层（M2，尚未实现）。

BlendedPath（笛卡尔）→ PathDerivatives（关节 q,q',q'',q'''，自适应采样 +
连续解 IK + Jacobian 链式法则 + 奇异处理）。见 docs/python_framework.md
§5.4、docs/robot_copp_design.md §5。
"""
