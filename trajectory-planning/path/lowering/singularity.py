"""奇异检测与阻尼最小二乘逆（framework §5.4 singularity / 设计 §5.4）。"""

from __future__ import annotations

import numpy as np


def min_singular_ratio(J: np.ndarray) -> float:
    """σ_min/σ_max（Jacobian 条件数的倒数）。0 = 完全奇异，1 = 完美条件。"""
    s = np.linalg.svd(J, compute_uv=False)
    if s[0] <= 0.0:
        return 0.0
    return float(s[-1] / s[0])


def damped_inverse_solve(J: np.ndarray, b: np.ndarray, lam: float = 0.05) -> np.ndarray:
    """阻尼最小二乘解 x = Jᵀ(JJᵀ + λ²I)⁻¹ b（奇异邻域的良态替代，设计 §5.4）。"""
    m = J.shape[0]
    return J.T @ np.linalg.solve(J @ J.T + (lam ** 2) * np.eye(m), b)
