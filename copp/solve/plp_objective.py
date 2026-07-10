"""PLP 分段线性目标（framework §5.6 plp_objective.py / 论文 eq.27、29d、Prop.3）。

真实时间被积函数 1/√a 用割线上包络逼近。每个 a_k 取采样点
δ_{k,0}<δ_{k,1}<…<δ_{k,P}，引入辅助变量 J_k 与约束（eq.29d）：

    J_k ≥ [δ_{i-1} + √(δ_{i-1}δ_i) + δ_i - a_k] / [(√δ_{i-1}+√δ_i)·√(δ_{i-1}δ_i)]

对每条割线 i=1..P。目标 min Σ w_k J_k；再加下界 a_k ≥ δ_{k,0}（Prop.3，
δ_{k,0} 足够小 → 不影响最优、却在 a→0⁺ 施加无穷惩罚，根除零进给奇异）。

默认采样点 δ_{k,l}=10^{l-4}·a_seed[k]（对齐论文 §6.1 实验设置）。
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import cvxpy as cp


@dataclass
class PlpObjective:
    deltas: np.ndarray   # (N, P+1) 每点采样点，列递增
    weights: np.ndarray  # (N,) 时间目标梯形权重

    @property
    def a_floor(self) -> np.ndarray:
        """下界 δ_{k,0}（Prop.3）。"""
        return self.deltas[:, 0]


def default_delta_samples(
    a_seed: np.ndarray, levels=(1e-4, 1e-3, 1e-2, 1e-1)
) -> np.ndarray:
    """δ_{k,l}=level·a_seed[k]，返回 (N, P+1)。"""
    a = np.maximum(a_seed, 1e-9)
    return np.column_stack([lvl * a for lvl in levels])


def build_plp(
    plp: PlpObjective, a: cp.Variable, jvar: cp.Variable
) -> tuple[list[cp.Constraint], cp.Expression]:
    """返回 (割线约束列表 + a≥δ0 下界, 目标表达式 Σ w_k J_k)。"""
    deltas = plp.deltas
    cons: list[cp.Constraint] = []

    # 割线约束 eq.29d：对每条割线 i，J ≥ (仿射于 a) 逐点成立
    for i in range(1, deltas.shape[1]):
        d0, d1 = deltas[:, i - 1], deltas[:, i]        # (N,)
        den = (np.sqrt(d0) + np.sqrt(d1)) * np.sqrt(d0 * d1)  # (N,) 常数
        num_const = d0 + np.sqrt(d0 * d1) + d1          # (N,)
        secant = cp.multiply(num_const - a, 1.0 / den)  # (num - a)/den
        cons.append(jvar >= secant)

    # 下界 δ_{k,0}（仅内部点；端点 a 由边界等式固定）
    cons.append(a[1:-1] >= plp.a_floor[1:-1])

    objective = plp.weights @ jvar
    return cons, objective
