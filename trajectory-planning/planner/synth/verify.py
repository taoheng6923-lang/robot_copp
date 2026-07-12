"""约束满足度校验：超限率 R_v / 超限时长比 D_v（framework §5.8 / 设计 §9，论文 §6.1.2）。

对合成后的时间域轨迹逐样本检查六类约束利用率（|·|/上限），统计：
  - R_v：超限样本数 / 总样本数；
  - D_v：超限样本的时间权重和 / 总时长（时间权重取 np.gradient(t)，
    与 R_v 的差别在非等距样本——段 seam 与末点）。
论文目标 R_v, D_v < 0.1%（网格点严格满足 + 区间内 O(Δ²) 违约界）。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from copp import RobotLimits

from .resample import TrajectoryResult


@dataclass
class VerifyMetrics:
    """超限指标。max_util 按约束类别给最大利用率（1.0=贴边）。"""

    r_v: float
    d_v: float
    max_util: dict[str, float] = field(default_factory=dict)
    ok: bool = True                 # r_v 与 d_v 均低于 threshold
    threshold: float = 1e-3        # 论文 §6.1.2 的 0.1%

    def summary(self) -> str:
        util = ", ".join(f"{k}={v:.3f}" for k, v in self.max_util.items())
        return (f"R_v={self.r_v:.2e}, D_v={self.d_v:.2e} "
                f"({'OK' if self.ok else 'VIOLATED'}); max_util: {util}")


def verify_limits(
    result: TrajectoryResult,
    limits: RobotLimits,
    util_tol: float = 1e-6,
    threshold: float = 1e-3,
) -> VerifyMetrics:
    """逐样本约束校验。utilisation > 1+util_tol 记为超限样本。

    检查项：轴向速度/加速度/jerk（逐关节各自上限）+ TCP 位置速度模/姿态
    角速度模（若 limits 配置了对应上界）。力矩不在此层（需动力学系数，
    属 Topp3Data/求解侧，M2+ 无真实动力学）。
    """
    n = result.n_axis
    vmax, amax, jmax = limits.axis_arrays(n)

    utils: dict[str, np.ndarray] = {
        "velocity": np.max(np.abs(result.qd) / vmax[:, None], axis=0),
        "acceleration": np.max(np.abs(result.qdd) / amax[:, None], axis=0),
        "jerk": np.max(np.abs(result.qddd) / jmax[:, None], axis=0),
    }
    if limits.v_tcp_max is not None:
        utils["tcp_velocity"] = result.v_tcp / float(limits.v_tcp_max)
    if limits.w_tcp_max is not None:
        utils["tcp_angular_velocity"] = result.w_tcp / float(limits.w_tcp_max)

    worst = np.max(np.vstack(list(utils.values())), axis=0)   # (T,) 逐样本最坏利用率
    violated = worst > 1.0 + util_tol

    weights = np.gradient(result.t)
    r_v = float(np.count_nonzero(violated)) / result.t.size
    d_v = float(np.sum(weights[violated])) / max(result.t_final, 1e-12)

    metrics = VerifyMetrics(
        r_v=r_v, d_v=d_v,
        max_util={k: float(np.max(v)) for k, v in utils.items()},
        threshold=threshold,
    )
    metrics.ok = r_v <= threshold and d_v <= threshold
    return metrics
