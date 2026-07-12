"""TrajectoryPlanner 门面（framework §5.10 / 设计 §10，M2+ 已实现）。

串起全流程用户入口：

    add_command(...)  累积指令（JointMove / LinearMove / CircularMove）
        │ plan(q_seed)
        ▼
    path.commands.build_sections   指令 → Section（纯几何）
    path.commands.lower_sections   G0 衔接校验 + seed 链逐段降维 → PathDerivatives
    copp.solve_splp                逐段 rest-to-rest 时间最优求解（TOTP-SPLP）
    synth.synthesize/concatenate   解析细剖面 → 等时间栅格轨迹，多段拼接
    synth.verify_limits            超限率 R_v / 超限时长比 D_v 校验

M2 语义：段间角点以停顿衔接（每段独立 rest-to-rest），故多段规划要求
limits 的边界为静止（a_bnd=b_bnd=(0,0)）；G2 blending（不停顿过渡）属 M3。
力矩约束需真实动力学系数（DynamicsModel/RNE，未实现），本门面暂不摄入。
HLAW 长序列分窗（M5）未实现——每段整段离线求解。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from robot import KinematicsModel
from copp import RobotLimits, solve_splp, SolveOptions
from copp.types import Profile

from path.types import PathDerivatives
from path.commands import MotionCommand, build_sections, lower_sections
from path.lowering import SampleOptions, IkOptions, DerivativeOptions

from .synth import (
    TrajectoryResult, VerifyMetrics, synthesize, concatenate, verify_limits,
)


@dataclass
class PlanOptions:
    """plan() 的全部可调参数（各子层选项透传）。"""

    dt: float = 1e-3
    solve: SolveOptions = field(default_factory=lambda: SolveOptions(n_iter=3))
    sample: SampleOptions | None = None
    ik: IkOptions | None = None
    deriv: DerivativeOptions | None = None
    verify: bool = True


@dataclass
class SegmentPlan:
    """单段规划产物（局部时间，t∈[0, result.t_final]）。"""

    pd: PathDerivatives
    profile: Profile
    t_final: float
    splp_t_final: list[float]          # SPLP 各迭代的终止时间（单调不增）
    result: TrajectoryResult


@dataclass
class PlanResult:
    """plan() 返回：逐段产物 + 拼接轨迹 + 校验指标。"""

    segments: list[SegmentPlan]
    trajectory: TrajectoryResult
    metrics: VerifyMetrics | None

    @property
    def t_final(self) -> float:
        return self.trajectory.t_final


class TrajectoryPlanner:
    """指令序列 → 时间最优关节轨迹的用户门面。

    用法：
        planner = (TrajectoryPlanner(kin, limits)
                   .add_command(JointMoveCommand(q0, q1))
                   .add_command(LinearMoveCommand(poseA, poseB)))
        res = planner.plan(q_seed=q0)
        res.trajectory.q / .qd / .qdd / .qddd, res.t_final, res.metrics
    """

    def __init__(self, kin: KinematicsModel, limits: RobotLimits):
        self._kin = kin
        self._limits = limits
        self._commands: list[MotionCommand] = []

    def add_command(self, cmd: MotionCommand) -> "TrajectoryPlanner":
        """累积一条指令（链式调用）。真正的求解在 plan() 一次性完成。"""
        self._commands.append(cmd)
        return self

    def clear(self) -> "TrajectoryPlanner":
        self._commands.clear()
        return self

    def plan(self, q_seed: np.ndarray, opts: PlanOptions | None = None) -> PlanResult:
        """执行全流程规划。q_seed：首段 IK 种子 / 首关节段起点校验值。"""
        if not self._commands:
            raise ValueError("指令队列为空：先 add_command 再 plan")
        opts = opts or PlanOptions()

        sections = build_sections(self._commands)
        if len(sections) > 1:
            a_bnd = np.asarray(self._limits.a_bnd, dtype=float)
            b_bnd = np.asarray(self._limits.b_bnd, dtype=float)
            if np.max(np.abs(a_bnd)) > 1e-12 or np.max(np.abs(b_bnd)) > 1e-12:
                raise ValueError(
                    "多段规划要求静止边界 a_bnd=b_bnd=(0,0)（M2 段间停顿语义）；"
                    f"当前 a_bnd={tuple(a_bnd)}, b_bnd={tuple(b_bnd)}"
                )

        pds = lower_sections(
            sections, self._kin, q_seed=q_seed,
            sample_opts=opts.sample, ik_opts=opts.ik, deriv_opts=opts.deriv,
        )

        segments: list[SegmentPlan] = []
        for pd in pds:
            data = self._limits.to_topp3_data(
                pd.s_grid, pd.dq, pd.ddq, pd.dddq, tcp_geom=pd.tcp_geom(),
            )
            profile, hist = solve_splp(data, opts.solve)
            result = synthesize(pd, profile, opts.dt)
            segments.append(SegmentPlan(
                pd=pd, profile=profile,
                t_final=result.t_final,
                splp_t_final=list(hist.t_final),
                result=result,
            ))

        trajectory = concatenate([sp.result for sp in segments])
        metrics = verify_limits(trajectory, self._limits) if opts.verify else None
        return PlanResult(segments=segments, trajectory=trajectory, metrics=metrics)
