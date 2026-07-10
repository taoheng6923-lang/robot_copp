"""统一的机器人约束配置（framework §5.5 的轻量前身）。

把散落在 demo/测试里的 vmax/amax/jmax/边界/TCP 限值收拢到一个 dataclass，
改约束只改一处；再由 `to_topp3_data` 结合路径几何（q',q'',q'''）产出求解输入。

- 标量或逐轴数组皆可（标量自动广播到 n 轴）。
- `from_ratios` 支持"加速度=速度×k_a、jerk=加速度×k_j"这类比例设定。
- TCP 限值（位置速度模 / 姿态角速度模）目前仅供可视化；M4 接入求解器约束后复用同一入口。
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from ..types import Topp3Data, TcpConstraint, TorqueConstraint

Scalar = float
Limit = "float | np.ndarray"


def _as_axis(val, n: int) -> np.ndarray:
    """标量 → 广播为 (n,)；数组 → 校验长度后原样返回（float）。"""
    arr = np.atleast_1d(np.asarray(val, dtype=float))
    if arr.size == 1:
        return np.full(n, float(arr[0]))
    if arr.size != n:
        raise ValueError(f"约束长度 {arr.size} 与轴数 {n} 不符")
    return arr


@dataclass
class RobotLimits:
    """机器人运动学约束配置。

    vmax/amax/jmax : 轴向速度/加速度/jerk 上界（对称）。标量或 (n,) 数组。
    a_bnd/b_bnd    : 两端 a=ṡ²、b=s̈ 边界。
    v_tcp_max      : TCP 位置速度模上界（可选，可视化/未来求解器约束）。
    w_tcp_max      : TCP 姿态角速度模上界（可选）。
    """

    vmax: Limit
    amax: Limit
    jmax: Limit
    a_bnd: tuple[float, float] = (0.0, 0.0)
    b_bnd: tuple[float, float] = (0.0, 0.0)
    v_tcp_max: float | None = None
    w_tcp_max: float | None = None
    tau_max: Limit | None = None          # M4：关节力矩上界
    tau_min: Limit | None = None          # M4：关节力矩下界（缺省取 -tau_max）

    @classmethod
    def from_ratios(
        cls,
        vmax: Limit,
        acc_ratio: float,
        jerk_ratio: float,
        a_bnd: tuple[float, float] = (0.0, 0.0),
        b_bnd: tuple[float, float] = (0.0, 0.0),
        v_tcp_max: float | None = None,
        w_tcp_max: float | None = None,
        tau_max: Limit | None = None,
        tau_min: Limit | None = None,
    ) -> "RobotLimits":
        """按比例设定：amax = acc_ratio·vmax，jmax = jerk_ratio·amax。"""
        v = np.asarray(vmax, dtype=float)
        a = acc_ratio * v
        j = jerk_ratio * a
        return cls(vmax=v, amax=a, jmax=j, a_bnd=a_bnd, b_bnd=b_bnd,
                   v_tcp_max=v_tcp_max, w_tcp_max=w_tcp_max,
                   tau_max=tau_max, tau_min=tau_min)

    def axis_arrays(self, n_axis: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """把 (vmax, amax, jmax) 广播/校验为 (n,) 数组三元组。"""
        return (
            _as_axis(self.vmax, n_axis),
            _as_axis(self.amax, n_axis),
            _as_axis(self.jmax, n_axis),
        )

    def to_topp3_data(
        self,
        s_grid: np.ndarray,
        dq: np.ndarray,
        ddq: np.ndarray,
        dddq: np.ndarray,
        tcp_geom: dict | None = None,
        torque_coeffs: dict | None = None,
        validate: bool = True,
    ) -> Topp3Data:
        """结合路径几何（q', q'', q'''）构造 Topp3Data。

        tcp_geom     : {"cv": (N,), "cw": (N,)}，配合 v_tcp_max/w_tcp_max 生成 TCP 约束。
        torque_coeffs: {"n_tor":(n,N), "m_tor":(n,N), "g_tor":(n,N)}，配合 tau_max/min 生成力矩约束。
        """
        n = dq.shape[0]
        vmax, amax, jmax = self.axis_arrays(n)

        tcp = None
        if tcp_geom is not None and self.v_tcp_max is not None and self.w_tcp_max is not None:
            tcp = TcpConstraint(
                cv=np.asarray(tcp_geom["cv"], float),
                cw=np.asarray(tcp_geom["cw"], float),
                v_max=float(self.v_tcp_max), w_max=float(self.w_tcp_max),
            )

        torque = None
        if torque_coeffs is not None and self.tau_max is not None:
            tau_max = _as_axis(self.tau_max, n)
            tau_min = _as_axis(self.tau_min, n) if self.tau_min is not None else -tau_max
            torque = TorqueConstraint(
                n_tor=np.asarray(torque_coeffs["n_tor"], float),
                m_tor=np.asarray(torque_coeffs["m_tor"], float),
                g_tor=np.asarray(torque_coeffs["g_tor"], float),
                tau_max=tau_max, tau_min=tau_min,
            )

        data = Topp3Data(
            s_grid=s_grid, dq=dq, ddq=ddq, dddq=dddq,
            vmax=vmax, amax=amax, jmax=jmax,
            a_bnd=self.a_bnd, b_bnd=self.b_bnd,
            tcp=tcp, torque=torque,
        )
        if validate:
            data.validate()
        return data
