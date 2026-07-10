"""路径构造模块（指令 → 笛卡尔路径 → 关节路径，M2/M3，尚未实现）。

指令序列（关节运动/直线/圆弧）经 blending（G2 Hermite 过渡）、lowering
（自适应采样 + IK + 链式法则求导）翻译为 copp 求解所需的 PathDerivatives。
与 `trajectory-planning/copp/`（纯数值求解核心）、顶层 `robot/`（运动学
动力学本体）相互独立，经 `trajectory-planning/planner/planner.py` 编排串联。
见 docs/python_framework.md §5.2-5.4、docs/robot_copp_design.md §3-5。

- commands/  指令层：JointMove / LinearMove / CircularMove → PoseSegment
- blending/  最优 Hermite blending（笛卡尔空间 G2 过渡）
- lowering/  降维：BlendedPath → 关节 PathDerivatives
"""
