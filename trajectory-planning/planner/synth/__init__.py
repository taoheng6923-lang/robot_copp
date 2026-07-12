"""轨迹合成与验证层（framework §5.8 / 设计 §9，M2+ 已实现）。

- resample: synthesize（单段：解析细剖面 → 等时间栅格 q,q̇,q̈,q⃛）+
            concatenate（多段拼接，段间 rest 停顿衔接）+ TrajectoryResult
- verify:   verify_limits（超限率 R_v / 超限时长比 D_v，论文 §6.1.2）+ VerifyMetrics
"""

from .resample import TrajectoryResult, synthesize, concatenate
from .verify import VerifyMetrics, verify_limits

__all__ = [
    "TrajectoryResult", "synthesize", "concatenate",
    "VerifyMetrics", "verify_limits",
]
