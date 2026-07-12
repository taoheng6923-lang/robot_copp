"""TOTP-SPLP 迭代（framework §5.6 splp.py / 论文 Algorithm 2）。

从种子 a⁽⁰⁾ 出发，第 p 次在 a^{(p-1)} 处线性化 jerk（eq.32）+ 解一次 PLP-LP，
得 a^{(p)}；直到停止准则（eq.30）：
    |t_f^{(p)} - t_f^{(p-1)}| < eps_t  或  ‖a^{(p)}-a^{(p-1)}‖ < eps_a  或  p ≥ n_iter。

PLP 采样点 δ 固定自种子（论文做法）；仅 jerk 线性化点随迭代更新（序列线性化）。
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

from ..types import Topp3Data, Profile
from ..options import ConstraintFlags
from .seed import compute_seed
from .plp_objective import PlpObjective, default_delta_samples
from .state import trapz_weights, resolve_num_stationary
from .lp_problem import build_and_solve
from .interp import s_to_t


@dataclass
class SolveOptions:
    n_iter: int = 4          # SPLP 最大迭代数（论文实践 2~4 足够）
    eps_t: float = 1e-4      # 终止时间收敛容差
    eps_a: float = 1e-6      # a 剖面收敛容差
    mode: str = "plp_lp"     # 默认 PLP+LP（算力最省）；"socp" 为备选（M1 未实现）
    n_stationary: int = 1    # 零进给端的静止段宽度 N_s（Box I / Prop.2）；仅对 a_bnd≈0 端生效
    flags: ConstraintFlags = field(default_factory=ConstraintFlags)  # 各约束启用开关（默认全开）
    smooth_c_weight: float = 0.0  # 非静止段 c 平滑惩罚权重 λ：目标 += λ·Σ|c_i−c_{i+1}|；0 关闭
    solver: str | None = None
    verbose: bool = False


@dataclass
class SplpHistory:
    t_final: list[float] = field(default_factory=list)
    converged: bool = False
    n_iter_run: int = 0


def solve_splp(
    data: Topp3Data,
    opts: SolveOptions | None = None,
    seed: np.ndarray | None = None,
) -> tuple[Profile, SplpHistory]:
    """论文 TOTP-SPLP。返回 (最终 Profile, 迭代历史)。"""
    data.validate()
    opts = opts or SolveOptions()
    if opts.mode != "plp_lp":
        raise NotImplementedError("M1 仅实现默认 PLP+LP 路线（mode='plp_lp'）")

    num_stat = resolve_num_stationary(data, opts.n_stationary)
    a_lin = (
        compute_seed(data, num_stat, opts.flags)
        if seed is None else np.asarray(seed).ravel()
    )
    plp = PlpObjective(
        deltas=default_delta_samples(a_lin),
        weights=trapz_weights(data.s_grid),
    )

    hist = SplpHistory()
    a_prev: np.ndarray | None = None
    tf_prev = np.inf
    profile: Profile | None = None

    for p in range(1, opts.n_iter + 1):
        profile = build_and_solve(
            data, a_lin, plp, solver=opts.solver, num_stat=num_stat, flags=opts.flags,
            smooth_c_weight=opts.smooth_c_weight,
        )
        tf, _ = s_to_t(data.s_grid, profile)
        hist.t_final.append(tf)
        hist.n_iter_run = p
        if opts.verbose:
            print(f"  SPLP iter {p}: t_final = {tf:.6f}")

        if a_prev is not None and (
            abs(tf - tf_prev) < opts.eps_t
            or np.linalg.norm(profile.a - a_prev) < opts.eps_a
        ):
            hist.converged = True
            break

        a_lin = np.maximum(profile.a, 1e-9)  # 重线性化（序列线性化）
        a_prev, tf_prev = profile.a, tf

    assert profile is not None
    return profile, hist
