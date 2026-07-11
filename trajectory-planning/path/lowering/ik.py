"""连续解 IK（framework §5.4 ik / 设计 §5.2）。

逐站点 `kin.ik(pose, seed=上一站点解)`，保证解分支连续（最近解选择由
KinematicsModel 协议约定）。两道防线：
  - FK 回代校验：IK 解的正运动学必须复现目标位姿（捕获不可达 / 退化解）；
  - 跳变校验：相邻站点 ‖Δq‖∞ 超阈值即报错（解分支切换 / 采样过疏 / 穿奇异点）。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from robot import KinematicsModel, Pose

from ..types import CartesianSamples
from ..errors import UnreachablePoseError, IkJumpError


@dataclass
class IkOptions:
    """连续解 IK 参数。

    jump_threshold : 相邻站点关节向量跳变上限 ‖Δq‖∞ [rad]。正常连续解在细网格
                     下步进 ≪ 此值；超限几乎必是解分支切换或穿越奇异点。
    fk_tol_pos     : FK 回代位置误差上限 [m]。
    fk_tol_rot     : FK 回代旋转矩阵误差上限（Frobenius 范数）。
    """

    jump_threshold: float = 0.5
    fk_tol_pos: float = 1e-6
    fk_tol_rot: float = 1e-6


def solve_ik_sequence(
    samples: CartesianSamples,
    s_grid: np.ndarray,
    kin: KinematicsModel,
    q_seed: np.ndarray,
    opts: IkOptions | None = None,
) -> np.ndarray:
    """逐站点连续解 IK，返回 q (n, N)。失败抛 UnreachablePoseError / IkJumpError。"""
    opts = opts or IkOptions()
    N = samples.n_grid
    q_prev = np.asarray(q_seed, dtype=float)
    q = np.zeros((q_prev.size, N))

    for k in range(N):
        pose = Pose(position=samples.p[:, k].copy(), rotation=samples.R[k].copy())
        qk = kin.ik(pose, seed=q_prev)

        chk = kin.fk(qk)
        pos_err = float(np.linalg.norm(chk.position - pose.position))
        rot_err = float(np.linalg.norm(chk.rotation - pose.rotation))
        if pos_err > opts.fk_tol_pos or rot_err > opts.fk_tol_rot:
            raise UnreachablePoseError(
                f"s={s_grid[k]:.6g} 处 IK 回代不符：pos_err={pos_err:.3e} m, "
                f"rot_err={rot_err:.3e}（目标位姿可能不可达）"
            )
        if k > 0:
            jump = float(np.max(np.abs(qk - q_prev)))
            if jump > opts.jump_threshold:
                raise IkJumpError(
                    f"s={s_grid[k]:.6g} 处关节解跳变 ‖Δq‖∞={jump:.3f} rad "
                    f"> {opts.jump_threshold}（解分支切换或采样过疏）"
                )
        q[:, k] = qk
        q_prev = qk
    return q
