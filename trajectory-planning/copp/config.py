"""从 YAML 配置文件加载机器人约束（→ RobotLimits）。

配置文件描述机器人本体的逐关节轴向约束（速度/加速度/jerk/力矩）与路径两端边界；
TCP 速度模上界（v_tcp_max/w_tcp_max）不在文件中，作为“给定”参数传入本函数
（属任务/工艺侧设定）。见 configs/robot_ur5.yaml 的字段说明。
"""

from __future__ import annotations

import os

import numpy as np

from .constraints import RobotLimits

_CONFIG_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "configs"
)
# 机器人本体约束（逐关节，默认 UR5；robot_3axis.yaml 保留作通用/教学示例）
DEFAULT_CONFIG = os.path.join(_CONFIG_DIR, "robot_ur5.yaml")
# 全局通用参数（跨模块 / 演示 / 示意等共享设置）
DEFAULT_COMM_CONFIG = os.path.join(_CONFIG_DIR, "comm_paras.yaml")


def _joint_field(joints: list[dict], key: str, required: bool = True):
    """抽取每关节的某字段为 (n,) 数组；required=False 且全缺时返回 None。"""
    present = [key in j for j in joints]
    if not any(present):
        if required:
            raise ValueError(f"配置每个关节都需字段 '{key}'")
        return None
    if not all(present):
        raise ValueError(f"字段 '{key}' 必须所有关节都给或都不给")
    return np.asarray([float(j[key]) for j in joints], dtype=float)


def _scale(cfg: dict, key: str) -> float:
    """读顶层倍率字段（实验旋钮，逐类独立；缺省 1.0，须为正）。"""
    val = float(cfg.get(key, 1.0))
    if val <= 0:
        raise ValueError(f"{key} 须为正，得到 {val}")
    return val


def load_robot_limits(
    path: str | os.PathLike | None = None,
    v_tcp_max: float | None = 0.6,
    w_tcp_max: float | None = 0.9,
) -> RobotLimits:
    """读取 YAML 机器人配置，返回 RobotLimits。

    path      : 配置文件路径；None 用 DEFAULT_CONFIG。
    v_tcp_max : TCP 位置速度模上界（给定，不在配置文件里）。None 表示不设。
    w_tcp_max : TCP 姿态角速度模上界（给定）。None 表示不设。

    读取字段：顶层 {n_axis?,vmax_scale?,amax_scale?,jmax_scale?,tau_scale?}、
    joints[].{vmax,amax,jmax,tau_max?,tau_min?,noload_speed?,viscous?,coulomb?}、
    boundary.{a_bnd,b_bnd}。tau_min 缺省取 -tau_max；力矩字段整体可缺（则不启用力矩上下界）。
    noload_speed（空载转速 ω0）给全才启用速度相关力矩（t–n）：τ0 复用 tau_max、拐点 ω_c 取 vmax、
    viscous/coulomb 缺省 0。

    四个 *_scale 倍率（缺省 1.0，须为正）是**各自独立**的实验旋钮，只缩放本类限值：
      vmax_scale → vmax（同时即 t–n 拐点 ω_c；**不**缩放 ω0/摩擦，故须 vmax_scale < min(ω0/vmax)）
      amax_scale → amax；jmax_scale → jmax（无耦合）
      tau_scale  → tau_max/tau_min（及复用 tau_max 的 t–n 平台 τ0）
    注意：t–n 启用时轴速上界取 ω0 而非 vmax（见 solve/state.velocity_upper_bound），
    此时 vmax_scale 只移动梯形拐点、并不放开速度上限——要放开速度须改 noload_speed。
    """
    import yaml

    path = os.fspath(path) if path is not None else DEFAULT_CONFIG
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    joints = cfg.get("joints")
    if not joints:
        raise ValueError(f"配置缺少非空 'joints' 列表：{path}")
    if "n_axis" in cfg and int(cfg["n_axis"]) != len(joints):
        raise ValueError(f"n_axis={cfg['n_axis']} 与 joints 数 {len(joints)} 不符")

    # 顶层倍率（各类独立的实验旋钮，缺省 1.0）：只缩放本类限值，互不牵连。
    vmax_raw = _joint_field(joints, "vmax")
    vmax_scale = _scale(cfg, "vmax_scale")
    vmax = vmax_raw * vmax_scale
    amax = _joint_field(joints, "amax") * _scale(cfg, "amax_scale")
    jmax = _joint_field(joints, "jmax") * _scale(cfg, "jmax_scale")
    tau_max = _joint_field(joints, "tau_max", required=False)
    tau_min = _joint_field(joints, "tau_min", required=False)
    tau_scale = _scale(cfg, "tau_scale")  # 统一乘 tau_max（及复用它的 t–n 平台 τ0）与 tau_min
    if tau_max is not None:
        tau_max = tau_max * tau_scale
        tau_min = -tau_max if tau_min is None else tau_min * tau_scale  # 缺省对称，否则同步缩放
    # 速度相关力矩（t–n）：空载转速 ω0 / 粘滞 Fv / 库仑 Fc（整体可缺；τ0=tau_max、ω_c=vmax）
    noload_speed = _joint_field(joints, "noload_speed", required=False)
    viscous = _joint_field(joints, "viscous", required=False)
    coulomb = _joint_field(joints, "coulomb", required=False)
    # vmax_scale 不缩放 ω0（用户约定），故放大 vmax 会把 t–n 拐点 ω_c(=vmax) 推向/推过 ω0。
    # 梯形要求 ω_c < ω0（types.Topp3Data.validate），此处提前拦下并给出可用上限。
    if noload_speed is not None and np.any(vmax >= noload_speed):
        i = int(np.argmax(vmax / noload_speed))
        cap = float(np.min(noload_speed / vmax_raw))
        raise ValueError(
            f"vmax_scale={vmax_scale:g} 过大：关节 '{joints[i].get('name', i)}' 缩放后 "
            f"vmax={vmax[i]:.4g} ≥ 空载转速 ω0={noload_speed[i]:.4g}，而 t–n 梯形要求拐点 "
            f"ω_c(=vmax) < ω0。请把 vmax_scale 降到 < {cap:.4g}，或调大各关节 noload_speed。"
            f"（注意 vmax_scale 只缩放 vmax/ω_c；t–n 启用时轴速上界取 ω0，不受 vmax_scale 影响）"
        )

    bnd = cfg.get("boundary", {}) or {}
    a_bnd = tuple(float(x) for x in bnd.get("a_bnd", (0.0, 0.0)))
    b_bnd = tuple(float(x) for x in bnd.get("b_bnd", (0.0, 0.0)))

    return RobotLimits(
        vmax=vmax, amax=amax, jmax=jmax, a_bnd=a_bnd, b_bnd=b_bnd,
        v_tcp_max=v_tcp_max, w_tcp_max=w_tcp_max,
        tau_max=tau_max, tau_min=tau_min,
        st_noload_speed=noload_speed, st_viscous=viscous, st_coulomb=coulomb,
    )


def load_comm_paras(path: str | os.PathLike | None = None) -> dict:
    """读取全局通用参数文件（configs/comm_paras.yaml），返回整个 dict（按节组织）。"""
    import yaml

    path = os.fspath(path) if path is not None else DEFAULT_COMM_CONFIG
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_constraint_flags(path: str | os.PathLike | None = None):
    """从全局参数取 `constraints` 节，返回 ConstraintFlags（各约束启用开关）。"""
    from .options import ConstraintFlags

    return ConstraintFlags.from_dict(load_comm_paras(path).get("constraints", {}))


def load_smooth_c_weight(path: str | os.PathLike | None = None) -> float:
    """从全局参数取 `objective.smooth_c_weight`（非静止段 c 平滑惩罚权重 λ）；缺省 0.0。"""
    return float(load_comm_paras(path).get("objective", {}).get("smooth_c_weight", 0.0))
