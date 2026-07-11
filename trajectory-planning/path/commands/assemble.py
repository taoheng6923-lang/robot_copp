"""指令序列装配与降维驱动（framework §5.2 / 设计 §3.5，M2 无 blending 版）。

M2 语义：相邻指令段只做 **G0 精确衔接校验**（位置/姿态/关节角必须对上，
不做几何过渡），每段独立降维、独立 rest-to-rest 求解——切向不连续的角点
以停顿衔接，物理上严格可行。G2 blending（不停顿的平滑过渡）属 M3。

IK seed 链：首个笛卡尔段用调用方 q_seed；其后各段用上一段末端关节角，
保证跨段解分支连续（设计 §12.2 的分支一致性）。
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from robot import KinematicsModel

from ..types import PathDerivatives
from ..errors import JunctionMismatchError
from ..lowering import (
    SampleOptions, IkOptions, DerivativeOptions, lower_cartesian, lower_joint,
)
from .base import MotionCommand, Section


def build_sections(commands: Sequence[MotionCommand]) -> list[Section]:
    """逐条翻译指令 → Section 列表（纯几何，无 IK）。"""
    return [cmd.to_section() for cmd in commands]


def lower_sections(
    sections: Sequence[Section],
    kin: KinematicsModel,
    q_seed: np.ndarray,
    sample_opts: SampleOptions | None = None,
    ik_opts: IkOptions | None = None,
    deriv_opts: DerivativeOptions | None = None,
    junction_tol_pos: float = 1e-6,
    junction_tol_rot: float = 1e-6,
    junction_tol_q: float = 1e-6,
) -> list[PathDerivatives]:
    """逐段降维，返回每段的 PathDerivatives（段间 G0 衔接校验 + seed 链）。

    q_seed：首段的 IK 种子（首段为关节段时校验其 q_start 与 q_seed 一致）。
    """
    results: list[PathDerivatives] = []
    cur_q = np.asarray(q_seed, dtype=float)

    for i, sec in enumerate(sections):
        if sec.native_space == "joint":
            gap = float(np.max(np.abs(cur_q - sec.q_start)))
            if gap > junction_tol_q:
                raise JunctionMismatchError(
                    f"段 {i}（关节段）起点关节角与上游不衔接：‖Δq‖∞={gap:.3e} rad"
                )
            pd = lower_joint(
                sec.path, kin,
                ds_max=sample_opts.ds_max if sample_opts is not None else None,
                deriv_opts=deriv_opts,
            )
            cur_q = np.asarray(sec.q_end, dtype=float).copy()
        else:
            if i > 0:
                fk = kin.fk(cur_q)
                pos_gap = float(np.linalg.norm(fk.position - sec.pose_start.position))
                rot_gap = float(np.linalg.norm(fk.rotation - sec.pose_start.rotation))
                if pos_gap > junction_tol_pos or rot_gap > junction_tol_rot:
                    raise JunctionMismatchError(
                        f"段 {i}（笛卡尔段）起点位姿与上游 FK 不衔接："
                        f"pos={pos_gap:.3e} m, rot={rot_gap:.3e}"
                    )
            pd = lower_cartesian(sec.path, kin, q_seed=cur_q,
                                 sample_opts=sample_opts, ik_opts=ik_opts,
                                 deriv_opts=deriv_opts)
            cur_q = pd.q[:, -1].copy()
        results.append(pd)
    return results
