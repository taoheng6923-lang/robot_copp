"""横切配置项（framework §5.10）。

当前只有 ConstraintFlags；BlendOptions / SampleOptions / WindowSpec / PlannerOptions
等留待对应里程碑（M3/M4/M5）落地时加入。

ConstraintFlags：每类约束都可单独启用/关闭；开关设置放在全局参数文件
configs/comm_paras.yaml 的 `constraints` 节，经 config.load_constraint_flags() 读取，
传入 SolveOptions.flags。关闭某约束仅表示求解时不施加它（其上/下界数据仍可存在，
只是不生效）。
"""

from __future__ import annotations

from dataclasses import dataclass, fields


@dataclass
class ConstraintFlags:
    """六类约束的启用开关（默认全开）。

    velocity             轴向速度 |q̇|≤vmax
    acceleration         轴向加速度 |q̈|≤amax
    jerk                 轴向 jerk |q⃛|≤jmax（3 阶）
    torque               关节力矩 τ_min≤τ≤τ_max
    tcp_velocity         TCP 位置速度模 ‖ṗ‖≤v_tcp_max
    tcp_angular_velocity TCP 姿态角速度模 ‖ω‖≤w_tcp_max
    """

    velocity: bool = True
    acceleration: bool = True
    jerk: bool = True
    torque: bool = True
    tcp_velocity: bool = True
    tcp_angular_velocity: bool = True

    @classmethod
    def from_dict(cls, d: dict | None) -> "ConstraintFlags":
        """由 dict 构造；未列出的项保持默认 True，未知键报错。"""
        d = d or {}
        known = {f.name for f in fields(cls)}
        unknown = set(d) - known
        if unknown:
            raise ValueError(f"未知约束开关：{sorted(unknown)}；可用：{sorted(known)}")
        return cls(**{k: bool(v) for k, v in d.items()})
