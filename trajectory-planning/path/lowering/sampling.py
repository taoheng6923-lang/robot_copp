"""曲率驱动自适应离散化（framework §5.4 sampling / 设计 §5.1）。

弦高误差模型：区间 [s, s+Δs] 内曲线偏离弦线 e ≈ ‖r''(s)‖·Δs²/8，令 e ≤ ε 得
Δs ≤ √(8ε/‖r''‖)。位置项用 ‖p''(s)‖（含切向分量，比只取法向分量保守），姿态项
用 ‖ω̂'(s)‖（姿态角对弦线的偏差，rad）。二者与 ds_max 取 min 得逐点允许步长，
再从左到右贪心铺网格。

网格**必须包含 path.s_breaks**（约束域 / 指令段边界）——这是论文 Theorem 1
的 O(Δ²) 离散化误差界成立的前提（设计 §5.1）。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..types import CartesianPath
from ..errors import PathError


@dataclass
class SampleOptions:
    """自适应采样参数。

    ds_max  : 最大步长；None 取 s_total/60。
    eps_pos : 位置弦高容差 [m]。
    eps_ori : 姿态弦差容差 [rad]。
    n_scan  : 曲率预扫描密度（均匀站点数）。
    n_min   : 全程最少站点数（不足则均匀加密补齐）。
    n_max   : 站点数上限（防曲率异常导致爆炸）。
    """

    ds_max: float | None = None
    eps_pos: float = 1e-4
    eps_ori: float = 1e-3
    n_scan: int = 2001
    n_min: int = 24
    n_max: int = 20000


def uniform_sample(s_total: float, ds_max: float) -> np.ndarray:
    """均匀网格（关节原生段用；站点数 = ceil(s_total/ds_max)+1，至少 2）。"""
    n = max(int(np.ceil(s_total / ds_max)) + 1, 2)
    return np.linspace(0.0, s_total, n)


def adaptive_sample(path: CartesianPath, opts: SampleOptions | None = None) -> np.ndarray:
    """曲率驱动采样，返回严格递增网格 (N,)，含 0、s_total 与全部 s_breaks。"""
    opts = opts or SampleOptions()
    L = float(path.s_total)
    if L <= 0.0:
        raise PathError(f"path.s_total={L} 必须为正")
    ds_max = opts.ds_max if opts.ds_max is not None else L / 60.0
    ds_floor = L / opts.n_max          # 步长硬下限（同时限制总站点数）

    # ── 预扫描逐点允许步长 ──────────────────────────────────────────────
    ss = np.linspace(0.0, L, opts.n_scan)
    smp = path.eval(ss)
    kappa_pos = np.linalg.norm(smp.ddp, axis=0)
    kappa_ori = np.linalg.norm(smp.dw, axis=0)
    with np.errstate(divide="ignore"):
        ds_p = np.where(kappa_pos > 1e-12, np.sqrt(8.0 * opts.eps_pos / kappa_pos), np.inf)
        ds_o = np.where(kappa_ori > 1e-12, np.sqrt(8.0 * opts.eps_ori / kappa_ori), np.inf)
    ds_allow = np.minimum(ds_max, np.minimum(ds_p, ds_o))

    # ── 分段边界（必含 0 与 L）────────────────────────────────────────────
    breaks = np.unique(np.clip(np.asarray(path.s_breaks, dtype=float), 0.0, L))
    if breaks.size == 0 or breaks[0] > 0.0:
        breaks = np.concatenate([[0.0], breaks])
    if breaks[-1] < L:
        breaks = np.concatenate([breaks, [L]])

    # ── 逐区间贪心铺点（步长取 [s, s+Δs] 两端允许值的 min，保守）─────────
    grid: list[float] = [0.0]
    for b0, b1 in zip(breaks[:-1], breaks[1:]):
        s = b0
        if s > grid[-1]:
            grid.append(s)
        while True:
            ds = float(np.interp(s, ss, ds_allow))
            ds = min(ds, float(np.interp(min(s + ds, b1), ss, ds_allow)))
            ds = max(ds, ds_floor)
            if s + ds >= b1 - 0.25 * ds:      # 尾步并入边界，避免残留极小步
                grid.append(b1)
                break
            s += ds
            grid.append(s)
            if len(grid) > opts.n_max:
                raise PathError(f"自适应采样超过 n_max={opts.n_max}（曲率异常？）")

    g = np.unique(np.asarray(grid))

    # ── 站点不足则均匀加密补齐（仍保留 breaks）──────────────────────────
    if g.size < opts.n_min:
        g = np.unique(np.concatenate([np.linspace(0.0, L, opts.n_min), breaks]))

    # 容差去重（防浮点近重合站点导致 Δs≈0）；若端点 L 被并掉则用 L 顶替末点
    keep = np.concatenate([[True], np.diff(g) > 1e-9 * max(L, 1.0)])
    g = g[keep]
    if g[-1] < L:
        g[-1] = L
    if not np.all(np.diff(g) > 0.0):
        raise PathError("采样网格非严格递增（去重后仍有重合站点）")
    return g
