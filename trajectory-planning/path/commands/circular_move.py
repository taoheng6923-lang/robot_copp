"""CircularMove 圆弧指令（framework §5.2 / 设计 §3.3）。

几何两种给法：
  - 三点定圆（via）：起点/终点/途经点确定唯一圆，扫角方向取"经过 via"的一侧；
  - 显式 (center, normal)：圆心 + 平面法向，direction 选 ccw/cw/shortest。

弧按弧长参数化（s ∈ [0, L]，L = r·|Φ|，φ(s) = Φ·s/L，σ = sign(Φ)）：
    p    = c + r(cosφ·e1 + sinφ·e2)
    p'   = σ(−sinφ·e1 + cosφ·e2)          （单位速度）
    p''  = −(cosφ·e1 + sinφ·e2)/r          （‖p''‖ = 1/r = 曲率）
    p''' = σ(sinφ·e1 − cosφ·e2)/r²
姿态 SLERP 同 LinearMove（ω̂ 恒定）。

退化情形显式报错（设计 §3.3）：三点共线 / 起终点重合 / 端点不在圆上或平面外。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from robot import Pose

from ..types import CartesianSamples
from ..errors import DegenerateArcError
from .base import Section, rotvec_between, slerp_frames

_TWO_PI = 2.0 * np.pi


@dataclass
class _ArcPath:
    """圆弧 + SLERP 的 CartesianPath 实现。"""

    center: np.ndarray
    e1: np.ndarray               # (p0−c)/r
    e2: np.ndarray               # n̂×e1
    radius: float
    sweep: float                 # 带符号扫角 Φ（>0 沿 e1→e2 方向）
    R0: np.ndarray
    rotvec: np.ndarray
    s_total: float

    def __post_init__(self):
        self.s_breaks = np.array([0.0, self.s_total])
        self._w = (self.R0 @ self.rotvec) / self.s_total

    def eval(self, s: np.ndarray) -> CartesianSamples:
        s = np.atleast_1d(np.asarray(s, dtype=float))
        N = s.size
        r, sg = self.radius, float(np.sign(self.sweep))
        phi = self.sweep * s / self.s_total
        c1, s1 = np.cos(phi), np.sin(phi)
        E1, E2 = self.e1[:, None], self.e2[:, None]
        return CartesianSamples(
            p=self.center[:, None] + r * (c1 * E1 + s1 * E2),
            dp=sg * (-s1 * E1 + c1 * E2),
            ddp=-(c1 * E1 + s1 * E2) / r,
            dddp=sg * (s1 * E1 - c1 * E2) / r**2,
            R=slerp_frames(self.R0, self.rotvec, s / self.s_total),
            w=np.repeat(self._w[:, None], N, axis=1),
            dw=np.zeros((3, N)), ddw=np.zeros((3, N)),
        )


def _circumcenter(p0: np.ndarray, via: np.ndarray, p1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """三点定圆：返回 (圆心 c, 单位法向 n̂)。共线/重合 → DegenerateArcError。

    圆心 c = p0 + αa + βb（a=via−p0, b=p1−p0），由 (c−p0)·a=‖a‖²/2、
    (c−p0)·b=‖b‖²/2 的 2×2 方程解出（垂直平分条件）。
    """
    a, b = via - p0, p1 - p0
    n = np.cross(a, b)
    scale = max(np.linalg.norm(a), np.linalg.norm(b))
    if scale < 1e-12 or np.linalg.norm(n) < 1e-9 * scale**2:
        raise DegenerateArcError("三点共线或重合，无法定圆")
    A = np.array([[a @ a, a @ b], [a @ b, b @ b]])
    rhs = 0.5 * np.array([a @ a, b @ b])
    alpha, beta = np.linalg.solve(A, rhs)
    return p0 + alpha * a + beta * b, n / np.linalg.norm(n)


def _angle_of(v: np.ndarray, e1: np.ndarray, e2: np.ndarray) -> float:
    """v 在 (e1,e2) 平面内的角，折到 [0, 2π)。"""
    return float(np.arctan2(v @ e2, v @ e1)) % _TWO_PI


@dataclass
class CircularMoveCommand:
    """圆弧指令：pose_start → pose_end，经 via 三点定圆，或显式 (center, normal)。

    via       : (3,) 途经点（与 center/normal 二选一）
    center    : (3,) 圆心（配合 normal）
    normal    : (3,) 圆弧平面法向（配合 center；扫角正向 = 绕 normal 右手 CCW）
    direction : 仅 (center,normal) 模式有效；via 模式方向由"经过 via"唯一确定
    rot_scale : 预留（弧长 L 由几何决定，纯姿态不可能是圆弧）
    """

    pose_start: Pose
    pose_end: Pose
    via: np.ndarray | None = None
    center: np.ndarray | None = None
    normal: np.ndarray | None = None
    direction: Literal["shortest", "ccw", "cw"] = "shortest"

    def to_section(self) -> Section:
        p0 = np.asarray(self.pose_start.position, dtype=float)
        p1 = np.asarray(self.pose_end.position, dtype=float)
        R0 = np.asarray(self.pose_start.rotation, dtype=float)
        R1 = np.asarray(self.pose_end.rotation, dtype=float)

        if self.via is not None:
            c, n_hat = _circumcenter(p0, np.asarray(self.via, dtype=float), p1)
            r = float(np.linalg.norm(p0 - c))
            e1 = (p0 - c) / r
            e2 = np.cross(n_hat, e1)
            phi_end = _angle_of(p1 - c, e1, e2)
            phi_via = _angle_of(np.asarray(self.via, dtype=float) - c, e1, e2)
            if phi_end < 1e-9:
                raise DegenerateArcError("圆弧起终点重合（不支持整圆）")
            sweep = phi_end if phi_via < phi_end else phi_end - _TWO_PI
        else:
            if self.center is None or self.normal is None:
                raise DegenerateArcError("需给 via 或 (center, normal) 之一")
            c = np.asarray(self.center, dtype=float)
            n_hat = np.asarray(self.normal, dtype=float)
            n_hat = n_hat / np.linalg.norm(n_hat)
            r = float(np.linalg.norm(p0 - c))
            if r < 1e-12:
                raise DegenerateArcError("起点与圆心重合")
            if abs((p0 - c) @ n_hat) > 1e-6 * r or abs((p1 - c) @ n_hat) > 1e-6 * r:
                raise DegenerateArcError("端点不在圆弧平面内")
            if abs(np.linalg.norm(p1 - c) - r) > 1e-6 * r:
                raise DegenerateArcError("终点不在圆上（|p1−c| ≠ r）")
            v0 = p0 - c
            e1 = v0 / r
            e2 = np.cross(n_hat, e1)
            phi_end = _angle_of(p1 - c, e1, e2)
            if phi_end < 1e-9:
                raise DegenerateArcError("圆弧起终点重合（不支持整圆）")
            if self.direction == "ccw":
                sweep = phi_end
            elif self.direction == "cw":
                sweep = phi_end - _TWO_PI
            else:
                sweep = phi_end if phi_end <= np.pi else phi_end - _TWO_PI

        L = r * abs(sweep)
        return Section(
            path=_ArcPath(center=c, e1=e1, e2=e2, radius=r, sweep=float(sweep),
                          R0=R0, rotvec=rotvec_between(R0, R1), s_total=L),
            native_space="cartesian",
            pose_start=self.pose_start, pose_end=self.pose_end,
        )
