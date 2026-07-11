"""关节导数链式法则 + 降维驱动（framework §5.4 derivatives / 设计 §5.3）。

几何 Jacobian 关系 J(q(s))·q'(s) = r₁(s) ≡ [p'(s); ω̂(s)] 逐阶对 s 求导：

    J q'   = r₁
    J q''  = r₂ − J' q'                    （r₂ = [p''; ω̂']）
    J q''' = r₃ − 2·J' q'' − J'' q'        （r₃ = [p'''; ω̂'']）

注意三阶式中 J'q'' 的系数是 **2**（乘积求导两次产生两项）——设计文档 §5.3
原式漏了该系数，本实现以逐阶求导推导为准，并在自测中用 q 序列的有限差分
交叉验证。

J'(s) = DJ[q']（J 对 q 的方向导数，方向 q'），J''(s) = D²J[q',q'] + DJ[q'']；
两者按方向有限差分计算（J 解析、光滑，步长见 DerivativeOptions）。

奇异站点（σ_min/σ_max < 阈值）改用阻尼最小二乘解并打 singular 标记（设计 §5.4）。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from robot import KinematicsModel

from ..types import CartesianSamples, CartesianPath, JointSpacePath, PathDerivatives
from .sampling import SampleOptions, adaptive_sample, uniform_sample
from .ik import IkOptions, solve_ik_sequence
from .singularity import min_singular_ratio, damped_inverse_solve


@dataclass
class DerivativeOptions:
    """链式法则求导参数。

    h1 : J 一阶方向差分的绝对步长（‖δq‖，rad）。
    h2 : J 二阶方向差分的绝对步长（二阶差分对舍入误差敏感，取更大）。
    sigma_ratio_singular : σ_min/σ_max 低于此值视为奇异 → DLS + 标记。
    dls_lambda           : DLS 阻尼系数 λ（设计 §5.4 建议 0.05）。
    """

    h1: float = 1e-6
    h2: float = 1e-4
    sigma_ratio_singular: float = 1e-6
    dls_lambda: float = 0.05


def _dir_diff(kin: KinematicsModel, q: np.ndarray, v: np.ndarray, h_abs: float) -> np.ndarray:
    """DJ[v]：J 沿方向 v 的一阶方向导数（中心差分，步长按 ‖v‖ 归一）。"""
    nv = float(np.linalg.norm(v))
    if nv < 1e-14:
        return np.zeros_like(kin.jacobian(q))
    h = h_abs / nv
    return (kin.jacobian(q + h * v) - kin.jacobian(q - h * v)) / (2.0 * h)


def _dir_diff2(
    kin: KinematicsModel, q: np.ndarray, v: np.ndarray, h_abs: float, J0: np.ndarray
) -> np.ndarray:
    """D²J[v,v]：J 沿方向 v 的二阶方向导数（中心二阶差分）。"""
    nv = float(np.linalg.norm(v))
    if nv < 1e-14:
        return np.zeros_like(J0)
    h = h_abs / nv
    return (kin.jacobian(q + h * v) - 2.0 * J0 + kin.jacobian(q - h * v)) / (h * h)


def joint_derivatives(
    samples: CartesianSamples,
    s_grid: np.ndarray,
    q: np.ndarray,
    kin: KinematicsModel,
    opts: DerivativeOptions | None = None,
) -> PathDerivatives:
    """由笛卡尔导数 + 已解 IK 的 q 序列，链式法则求 q', q'', q'''。"""
    opts = opts or DerivativeOptions()
    n, N = q.shape
    dq = np.zeros((n, N))
    ddq = np.zeros((n, N))
    dddq = np.zeros((n, N))
    singular = np.zeros(N, dtype=bool)

    r1 = np.vstack([samples.dp, samples.w])       # (6, N)
    r2 = np.vstack([samples.ddp, samples.dw])
    r3 = np.vstack([samples.dddp, samples.ddw])

    for k in range(N):
        qk = q[:, k]
        J = kin.jacobian(qk)
        sing = min_singular_ratio(J) < opts.sigma_ratio_singular
        singular[k] = sing

        if sing:
            def solve(b, _J=J):
                return damped_inverse_solve(_J, b, opts.dls_lambda)
        else:
            def solve(b, _J=J):
                return np.linalg.solve(_J, b)

        qp = solve(r1[:, k])
        dJ = _dir_diff(kin, qk, qp, opts.h1)                       # J'(s) = DJ[q']
        qpp = solve(r2[:, k] - dJ @ qp)
        ddJ = _dir_diff2(kin, qk, qp, opts.h2, J) + _dir_diff(kin, qk, qpp, opts.h1)
        qppp = solve(r3[:, k] - 2.0 * dJ @ qpp - ddJ @ qp)

        dq[:, k], ddq[:, k], dddq[:, k] = qp, qpp, qppp

    return PathDerivatives(
        s_grid=np.asarray(s_grid, dtype=float),
        q=q, dq=dq, ddq=ddq, dddq=dddq,
        singular=singular,
        cv=np.linalg.norm(samples.dp, axis=0),
        cw=np.linalg.norm(samples.w, axis=0),
    )


def lower_cartesian(
    path: CartesianPath,
    kin: KinematicsModel,
    q_seed: np.ndarray,
    sample_opts: SampleOptions | None = None,
    ik_opts: IkOptions | None = None,
    deriv_opts: DerivativeOptions | None = None,
    s_grid: np.ndarray | None = None,
) -> PathDerivatives:
    """笛卡尔路径 → 关节 PathDerivatives 全流程（采样 → 连续 IK → 链式求导）。

    s_grid 显式给定时跳过自适应采样（测试 / 网格复用场景）。
    """
    if s_grid is None:
        s_grid = adaptive_sample(path, sample_opts)
    samples = path.eval(np.asarray(s_grid, dtype=float))
    q = solve_ik_sequence(samples, s_grid, kin, q_seed, ik_opts)
    return joint_derivatives(samples, s_grid, q, kin, deriv_opts)


def lower_joint(
    jpath: JointSpacePath,
    kin: KinematicsModel,
    ds_max: float | None = None,
    deriv_opts: DerivativeOptions | None = None,
    s_grid: np.ndarray | None = None,
) -> PathDerivatives:
    """关节原生路径 → PathDerivatives 快路径（无 IK；设计 §4.5 快路径）。

    q 及导数解析已知；TCP 速度模系数经真实 Jacobian 求 [p'; ω̂] = J q'。
    """
    opts = deriv_opts or DerivativeOptions()
    L = float(jpath.s_total)
    if s_grid is None:
        s_grid = uniform_sample(L, ds_max if ds_max is not None else L / 60.0)
    q, dq, ddq, dddq = jpath.eval_joint(np.asarray(s_grid, dtype=float))

    N = s_grid.size
    cv = np.zeros(N)
    cw = np.zeros(N)
    singular = np.zeros(N, dtype=bool)
    for k in range(N):
        J = kin.jacobian(q[:, k])
        singular[k] = min_singular_ratio(J) < opts.sigma_ratio_singular
        twist = J @ dq[:, k]
        cv[k] = float(np.linalg.norm(twist[:3]))
        cw[k] = float(np.linalg.norm(twist[3:]))

    return PathDerivatives(
        s_grid=np.asarray(s_grid, dtype=float),
        q=q, dq=dq, ddq=ddq, dddq=dddq,
        singular=singular, cv=cv, cw=cw,
    )
