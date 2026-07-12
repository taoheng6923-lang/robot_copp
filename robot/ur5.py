"""UR5 六轴机器人模型：真实 DH 运动学 + 力矩/速度限值（M1→真实本体的第一步）。

运动学（`UR5Kinematics`）严格按 Universal Robots 官方标准 DH 参数表实现
（https://www.universal-robots.com/articles/ur/application-installation/
dh-parameters-for-calculations-of-kinematics-and-dynamics/，标准 DH 约定）：
`fk`/`jacobian` 为解析解，`ik` 用闭式解析逆解枚举 UR5 全部 8 支路并取离 seed
最近者（O(1)，退化位姿回退 DLS 兜底），`jacobian_derivative` 按
`robot.base.KinematicsModel` 协议用有限差分实现。实现 `robot.base` 的
`KinematicsModel` 协议。

动力学（`UR5RobotModel.torque_coeffs`）仍是**对角/无耦合近似**（M2+ 的
`DynamicsModel`/真实 RNE 落地前的 stand-in）：用 ROS-Industrial
`ur5.urdf.xacro` 公开的连杆质量/长度做"下游集中质量单摆臂"式估算，量级贴近
真实 UR5，但不是精确 RNE（无科氏力/惯量耦合项）。`joint_path(s)` 仍是合成
随机轨迹（无路径生成层前的占位），只是关节数固定为 6、幅值裁到 UR5 关节
限位附近且绕一个非奇异"家位姿"摆动；`tcp_geometry(s)` 改为**真实**调用
`UR5Kinematics.jacobian(q)` 沿该合成路径求 TCP 线/角速度系数，不再是与 q
无关的虚构公式（这是相对旧 `SyntheticRobotModel` 的关键差异）。

参考数据来源（准确度分三档，见各常量注释）：
- DH 参数：Universal Robots 官网 DH 参数页（标准 DH 约定）—— **官方权威**。
- 关节最大力矩 TAU_MAX：Universal Robots 官网 "Max. joint torques CB3 and
  e-Series"（[54, 150, 150, 28, 28, 9] Nm）—— **官方权威**。
- 关节最大速度 V_MAX：ROS-Industrial `universal_robot` 仓库 `ur5.urdf.xacro`
  的 `<limit velocity=.../>`（[3.15,3.15,3.15,3.2,3.2,3.2] rad/s）—— **社区
  维护、驱动真实 UR5 硬件的描述文件，可信但非 UR 官方数据表原文**。
- 连杆质量/长度：同一 `ur5.urdf.xacro` 的 mass/length 属性 —— 同上可信度，
  仅用于 torque_coeffs 的近似估算（非精确 RNE 输入）。
- 关节最大加速度/jerk AMAX_RATIO/JMAX_RATIO：**无任何官方或社区公开数据**
  （UR 数据表、`ur5.urdf.xacro` 均只给速度/力矩限，不含加速度/jerk），此处
  按 vmax 的 10×/200× 经验比例**臆造**，仅供数值示例，不代表真实机器人规格。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .base import Pose

# ── UR5 标准 DH 参数（Universal Robots 官网，标准 DH 约定） ─────────────
# 关节顺序：base(1) shoulder(2) elbow(3) wrist1(4) wrist2(5) wrist3(6)
DH_A = np.array([0.0, -0.425, -0.39225, 0.0, 0.0, 0.0])
DH_D = np.array([0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823])
DH_ALPHA = np.array([np.pi / 2, 0.0, 0.0, np.pi / 2, -np.pi / 2, 0.0])
N_AXIS = 6

# ── 关节限值（TAU_MAX 官方权威；V_MAX 来自 ur5.urdf.xacro，逐关节不同） ──
TAU_MAX = np.array([54.0, 150.0, 150.0, 28.0, 28.0, 9.0])                    # Nm
V_MAX = np.array([3.15, 3.15, 3.15, 3.2, 3.2, 3.2])                          # rad/s

# ── 连杆质量/长度（ROS-Industrial ur5.urdf.xacro 公开参数，供 torque_coeffs 近似估算）
LINK_MASS = np.array([3.7, 8.393, 2.275, 1.219, 1.219, 0.1879])       # kg，joint i 之后的连杆
LINK_LENGTH = np.array([0.089159, 0.425, 0.39225, 0.09540, 0.09465, 0.0823])  # m，近似臂长
_GRAVITY = 9.81


def _dh_transform(theta: float, d: float, a: float, alpha: float) -> np.ndarray:
    """标准 DH 单步齐次变换 A = Rot_z(theta)·Trans_z(d)·Trans_x(a)·Rot_x(alpha)。"""
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [0.0,       sa,       ca,      d],
        [0.0,      0.0,      0.0,    1.0],
    ])


@dataclass
class UR5Kinematics:
    """UR5 正/逆运动学 + 几何 Jacobian（实现 `robot.base.KinematicsModel` 协议）。

    `fk`/`jacobian` 为标准 DH 解析解；`ik` 为**闭式解析逆解**：枚举 UR5 全部
    ≤8 支路（Andersen 标准 DH 推导），按逐关节 2π 折叠后取离 seed 最近者，
    天然满足协议要求的"连续解选择"，O(1) 求解、比旧版 DLS 数值迭代快 1~2 个
    数量级；极少数退化位姿回退 `_ik_dls` 兜底。
    """

    dof: int = N_AXIS

    def _frames(self, q: np.ndarray) -> list[np.ndarray]:
        """T[0..6]：T[0]=基座（单位阵），T[k] 为到第 k 个 DH 坐标系的齐次变换。"""
        T = [np.eye(4)]
        for i in range(self.dof):
            A = _dh_transform(q[i], DH_D[i], DH_A[i], DH_ALPHA[i])
            T.append(T[-1] @ A)
        return T

    def fk(self, q: np.ndarray) -> Pose:
        T = self._frames(q)[-1]
        return Pose(position=T[:3, 3].copy(), rotation=T[:3, :3].copy())

    def jacobian(self, q: np.ndarray) -> np.ndarray:
        """几何 Jacobian (6, 6)：行 0-2 线速度、行 3-5 角速度。"""
        T = self._frames(q)
        o_n = T[-1][:3, 3]
        J = np.zeros((6, self.dof))
        for i in range(self.dof):
            z = T[i][:3, 2]
            o = T[i][:3, 3]
            J[:3, i] = np.cross(z, o_n - o)
            J[3:, i] = z
        return J

    def jacobian_derivative(self, q: np.ndarray, dq: np.ndarray, h: float = 1e-6) -> np.ndarray:
        """dJ/ds 沿方向 dq 的中心差分，(6, n)（见协议文档）。"""
        norm = np.linalg.norm(dq)
        if norm < 1e-12:
            return np.zeros((6, self.dof))
        step = dq / norm * h
        return (self.jacobian(q + step) - self.jacobian(q - step)) / (2.0 * h)

    def ik(self, pose: Pose, seed: np.ndarray) -> np.ndarray:
        """UR5 解析 IK：闭式枚举全部 ≤8 支路，返回离 seed 最近的关节解。

        相比旧版阻尼最小二乘（DLS）逐点牛顿迭代（≤200 步），解析法为
        O(1) 闭式求解，快 1~2 个数量级且不受初值/奇异收敛性影响。解按
        "逐关节 2π 折叠到最接近 seed" 后取整体 ‖Δq‖ 最小者，天然满足协议
        要求的连续解选择。极少数无解析解的退化位姿（不可达/腕奇异导致
        全支路 NaN）回退到 DLS 兜底。
        """
        seed = np.asarray(seed, dtype=float)
        sols = self._ik_analytic(pose)
        if sols:
            wrapped = [self._wrap_to_seed(q, seed) for q in sols]
            return min(wrapped, key=lambda q: np.linalg.norm(q - seed))
        return self._ik_dls(pose, seed)

    @staticmethod
    def _wrap_to_seed(q: np.ndarray, seed: np.ndarray) -> np.ndarray:
        """逐关节加减 2π 的整数倍，折叠到离 seed 最近（连续解、避免整圈跳变）。"""
        return q + 2.0 * np.pi * np.round((seed - q) / (2.0 * np.pi))

    def _ik_analytic(self, pose: Pose) -> list[np.ndarray]:
        """闭式枚举 UR5 全部逆解（Andersen 标准 DH 推导），返回有效解列表（≤8）。

        约定与本类 `fk`/DH 参数完全一致：θ1(2 支)×θ5(2 支)×θ3(2 支)=8 支，
        θ6/θ2/θ4 由对应支路唯一确定。腕奇异（sinθ5≈0）时 θ6 不定，置 0
        （随后由 `_wrap_to_seed`/最近解选择消化）。域外/退化支路（acos 越界
        产生 NaN）被过滤。
        """
        d1, d4, d5, d6 = DH_D[0], DH_D[3], DH_D[4], DH_D[5]
        a2, a3 = DH_A[1], DH_A[2]

        T06 = np.eye(4)
        T06[:3, :3] = pose.rotation
        T06[:3, 3] = pose.position
        p06x, p06y = T06[0, 3], T06[1, 3]

        sols: list[np.ndarray] = []

        # ── θ1：肩关节两支 ─────────────────────────────────────────────
        p05 = T06 @ np.array([0.0, 0.0, -d6, 1.0])
        r1 = np.hypot(p05[0], p05[1])
        if r1 < abs(d4):                        # d4 圆柱外：θ1 无实解
            return sols
        psi, phi = np.arctan2(p05[1], p05[0]), np.arccos(np.clip(d4 / r1, -1.0, 1.0))
        for t1 in (psi + phi + np.pi / 2.0, psi - phi + np.pi / 2.0):
            s1, c1 = np.sin(t1), np.cos(t1)

            # ── θ5：腕两支 ─────────────────────────────────────────────
            c5 = (p06x * s1 - p06y * c1 - d4) / d6
            if abs(c5) > 1.0:                   # 该 θ1 下 θ5 无实解
                continue
            for t5 in (np.arccos(c5), -np.arccos(c5)):
                s5 = np.sin(t5)

                # ── θ6：腕奇异时不定，置 0 ───────────────────────────
                if abs(s5) < 1e-8:
                    t6 = 0.0
                else:
                    t6 = np.arctan2((-T06[0, 1] * s1 + T06[1, 1] * c1) / s5,
                                    (T06[0, 0] * s1 - T06[1, 0] * c1) / s5)

                # ── 由 θ1,θ5,θ6 反解出平面两连杆(a2,a3)问题 ─────────────
                T01 = _dh_transform(t1, d1, DH_A[0], DH_ALPHA[0])
                T45 = _dh_transform(t5, d5, DH_A[4], DH_ALPHA[4])
                T56 = _dh_transform(t6, d6, DH_A[5], DH_ALPHA[5])
                # T14 = A2·A3·A4（α2=α3=0）：原点退化为 x-y 平面两连杆(a2,a3)，
                # z≡d4 常量；旋转部 R14 = Rz(θ2+θ3+θ4)·Rx(π/2) → θ234 可直接取出
                T14 = np.linalg.inv(T01) @ T06 @ np.linalg.inv(T56) @ np.linalg.inv(T45)
                p14x, p14y = T14[0, 3], T14[1, 3]
                r2 = p14x ** 2 + p14y ** 2

                # ── θ3：肘两支（elbow up/down），平面两连杆余弦定理 ─────
                c3 = (r2 - a2 ** 2 - a3 ** 2) / (2.0 * a2 * a3)
                if abs(c3) > 1.0:               # 不可达
                    continue
                t234 = np.arctan2(T14[1, 0], T14[0, 0])   # θ2+θ3+θ4
                for t3 in (np.arccos(c3), -np.arccos(c3)):
                    # ── θ2（标准 2R 逆解）、θ4（由 θ234 反推） ──────────
                    t2 = np.arctan2(p14y, p14x) - np.arctan2(a3 * np.sin(t3), a2 + a3 * np.cos(t3))
                    t4 = t234 - t2 - t3

                    q = np.array([t1, t2, t3, t4, t5, t6])
                    if np.all(np.isfinite(q)):
                        sols.append(q)
        return sols

    def _ik_dls(
        self, pose: Pose, seed: np.ndarray,
        max_iter: int = 200, tol: float = 1e-8, damping: float = 1e-2,
    ) -> np.ndarray:
        """阻尼最小二乘数值 IK（解析法的退化位姿兜底）：从 seed 牛顿迭代。"""
        q = np.asarray(seed, dtype=float).copy()
        for _ in range(max_iter):
            cur = self.fk(q)
            pos_err = pose.position - cur.position
            R_err = pose.rotation @ cur.rotation.T
            rot_err = 0.5 * np.array([
                R_err[2, 1] - R_err[1, 2],
                R_err[0, 2] - R_err[2, 0],
                R_err[1, 0] - R_err[0, 1],
            ])
            err = np.concatenate([pos_err, rot_err])
            if np.linalg.norm(err) < tol:
                break
            J = self.jacobian(q)
            JJt = J @ J.T + (damping ** 2) * np.eye(6)
            q = q + J.T @ np.linalg.solve(JJt, err)
        return q


@dataclass
class UR5RobotModel:
    """UR5 本体：真实 DH 运动学（TCP 几何）+ 近似对角动力学（力矩系数）。

    关节数固定为 6（UR5）。`joint_path` 仍是合成随机轨迹（供 SPLP 求解/
    可视化冒烟测试，幅值绕一个非奇异"家位姿"小幅摆动）；`tcp_geometry`/
    `tcp_coeffs` 真实调用 `UR5Kinematics.jacobian` 沿该路径求值；
    `torque_coeffs` 是对角近似（stand-in，见模块 docstring）。

    seed  合成关节路径的随机种子
    kin   UR5Kinematics 实例（默认新建）
    """

    seed: int = 3
    kin: UR5Kinematics = field(default_factory=UR5Kinematics)
    n_axis: int = field(default=N_AXIS, init=False)

    # UR5 常见非奇异"家位姿"（肘部弯曲，避免落在肩/肘完全伸直的奇异位形）
    _HOME = np.array([0.0, -np.pi / 2, np.pi / 2, -np.pi / 2, -np.pi / 2, 0.0])

    # 合成路径速度尺度（不代表真实 UR5；越大关节转得越快，用于让 t–n 的高速 rolloff 段被激活）
    path_amp_scale: float = 1.0
    path_freq: tuple = (1.0, 2.0)

    # ── 运动学：合成关节路径（lowering/IK 落地前的占位） ─────────────────
    def joint_path(self, s: np.ndarray):
        """返回 (q0, q1, q2, q3)，各 (6, N)：q(s) 及其 1/2/3 阶导。"""
        rng = np.random.default_rng(self.seed)
        amp = 0.15 * np.pi * self.path_amp_scale  # 幅值（默认取关节全量程一小部分，避奇异）
        A = rng.uniform(0.5, 1.0, self.n_axis) * amp
        w = rng.uniform(self.path_freq[0], self.path_freq[1], self.n_axis)
        phi = rng.uniform(0.0, 2.0 * np.pi, self.n_axis)
        th = w[:, None] * s[None, :] + phi[:, None]
        q0 = self._HOME[:, None] + A[:, None] * np.sin(th)
        q1 = A[:, None] * w[:, None] * np.cos(th)
        q2 = -A[:, None] * w[:, None] ** 2 * np.sin(th)
        q3 = -A[:, None] * w[:, None] ** 3 * np.cos(th)
        return q0, q1, q2, q3

    # ── 运动学：TCP 几何（真实 UR5 Jacobian，替代旧版与 q 无关的虚构公式） ──
    def tcp_geometry(self, s: np.ndarray) -> dict:
        """返回 {dp, wdir}，各 (3, N)：逐点用真实 UR5 Jacobian 求
        dp/ds = J_v(q)·dq/ds、dw/ds = J_ω(q)·dq/ds（沿合成关节路径）。
        """
        q0, q1, _, _ = self.joint_path(s)
        N = s.size
        dp = np.zeros((3, N))
        wdir = np.zeros((3, N))
        for k in range(N):
            J = self.kin.jacobian(q0[:, k])
            dp[:, k] = J[:3, :] @ q1[:, k]
            wdir[:, k] = J[3:, :] @ q1[:, k]
        return {"dp": dp, "wdir": wdir}

    def tcp_coeffs(self, s: np.ndarray) -> dict:
        """返回 {cv, cw}，各 (N,)：TCP 速度模系数（供 to_topp3_data(tcp_geom=...)）。"""
        g = self.tcp_geometry(s)
        return {
            "cv": np.linalg.norm(g["dp"], axis=0),
            "cw": np.linalg.norm(g["wdir"], axis=0),
        }

    # ── 动力学：力矩系数（对角近似，见模块 docstring） ──────────────────
    def torque_coeffs(self, q0: np.ndarray, q1: np.ndarray, q2: np.ndarray) -> dict:
        """τ ≈ n_tor·a + m_tor·b + g_tor 的对角近似系数（非精确 RNE）。

        用"下游集中质量单摆臂"估算各关节等效惯量/重力力矩量级：下游（含自身）
        连杆质量之和 × 下游连杆长度之和作为力臂，量级贴近真实 UR5 关节力矩
        （如此估算下肩关节 ≈140 Nm，对照官方 150 Nm 上限）。关节 1（base）
        绕竖直轴（yaw）转动，典型位形下重力不产生绕该轴的力矩，故其重力项
        置零；惯量项仍保留（yaw 运动仍需克服下游连杆的转动惯量）。
        """
        down_mass = np.cumsum(LINK_MASS[::-1])[::-1]      # (6,) 下游（含自身）连杆质量和
        down_reach = np.cumsum(LINK_LENGTH[::-1])[::-1]   # (6,) 下游连杆长度和（近似力臂）
        inertia = down_mass * down_reach ** 2
        gravity_scale = down_mass * _GRAVITY * down_reach
        gravity_scale[0] = 0.0  # 关节1为竖直 yaw 轴，重力不产生绕该轴的力矩
        return {
            "n_tor": inertia[:, None] * q2,
            "m_tor": inertia[:, None] * q1,
            "g_tor": gravity_scale[:, None] * np.sin(q0),
        }
