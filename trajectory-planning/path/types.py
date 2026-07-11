"""路径构造层共享数据类型（framework §4 的 path 子集）。

约定：
- 全局路径参数 s ∈ [0, s_total]，近似弧长（每条指令段各自归一）。
- 姿态导数用**世界系角速度密度** ω̂(s) = ω/ṡ 表达（Ṙ Rᵀ = [ω]×，ω̂ = vee(R'Rᵀ)），
  与几何 Jacobian 的角速度行同一坐标系，链式法则 J(q)·q' = [p'; ω̂] 才成立。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np


@dataclass
class CartesianSamples:
    """位姿路径在一批站点 s 上的值与 1~3 阶导（全部对全局参数 s）。

    p / dp / ddp / dddp : (3, N)   TCP 位置及导数
    R                   : (N,3,3)  TCP 姿态旋转矩阵
    w / dw / ddw        : (3, N)   世界系 ω̂ = vee(R'Rᵀ) 及其 1/2 阶 s 导
    """

    p: np.ndarray
    dp: np.ndarray
    ddp: np.ndarray
    dddp: np.ndarray
    R: np.ndarray
    w: np.ndarray
    dw: np.ndarray
    ddw: np.ndarray

    @property
    def n_grid(self) -> int:
        return self.p.shape[1]


@runtime_checkable
class CartesianPath(Protocol):
    """笛卡尔位姿路径：可在任意 s 批量求值（含 1~3 阶导）。

    s_breaks 为段边界（含 0 与 s_total）；自适应采样网格必须包含它们
    （论文 Theorem 1 的 O(Δ²) 误差界前提，设计 §5.1）。
    """

    s_total: float
    s_breaks: np.ndarray

    def eval(self, s: np.ndarray) -> CartesianSamples: ...


@runtime_checkable
class JointSpacePath(Protocol):
    """关节空间原生路径（JointMove 快路径：无需 IK，解析给出 q 及导数）。"""

    s_total: float

    def eval_joint(
        self, s: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """返回 (q, dq, ddq, dddq)，各 (n, N)。"""
        ...


@dataclass
class PathDerivatives:
    """关节路径导数 + TCP 速度模系数——喂给 copp 求解的最终几何量（设计 §5.5）。

    s_grid          : (N,)   路径参数网格（严格递增）
    q/dq/ddq/dddq   : (n,N)  关节角及 1~3 阶 s 导
    singular        : (N,)   bool，该站点 Jacobian 奇异（已用 DLS 求解）标记
    cv              : (N,)   ‖p'(s)‖   —— TCP 位置速度模系数（tcp_geom 的 cv）
    cw              : (N,)   ‖ω̂(s)‖  —— TCP 姿态角速度模系数（tcp_geom 的 cw）
    """

    s_grid: np.ndarray
    q: np.ndarray
    dq: np.ndarray
    ddq: np.ndarray
    dddq: np.ndarray
    singular: np.ndarray
    cv: np.ndarray
    cw: np.ndarray

    @property
    def n_axis(self) -> int:
        return self.q.shape[0]

    @property
    def n_grid(self) -> int:
        return self.s_grid.size

    def tcp_geom(self) -> dict:
        """直接可传 `RobotLimits.to_topp3_data(tcp_geom=...)` 的 {cv, cw}。"""
        return {"cv": self.cv, "cw": self.cw}
