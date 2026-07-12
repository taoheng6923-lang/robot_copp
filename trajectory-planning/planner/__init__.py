"""调度/合成层（framework §5.7-5.10 / 设计 §8-10）。

编排 `trajectory-planning/path/` → 顶层 `robot/` → `trajectory-planning/copp/`
全流程并产出可执行轨迹：

- planner.py  TrajectoryPlanner 门面（M2+ 已实现：指令 → 逐段降维 → 逐段
              rest-to-rest SPLP → 等时间栅格轨迹拼接 + R_v/D_v 校验）
- synth/      轨迹合成与验证（M2+ 已实现）
- hlaw/       长序列分层前瞻窗口调度（M5，未实现——目前每段整段离线求解）
"""

from .planner import TrajectoryPlanner, PlanOptions, PlanResult, SegmentPlan
from .synth import (
    TrajectoryResult, VerifyMetrics, synthesize, concatenate, verify_limits,
)

__all__ = [
    "TrajectoryPlanner", "PlanOptions", "PlanResult", "SegmentPlan",
    "TrajectoryResult", "VerifyMetrics", "synthesize", "concatenate", "verify_limits",
]
