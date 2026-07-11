# 转矩–转速（t–n）约束的凸化处理 — Ardeshiri 论文精读

> **论文**：Tohid Ardeshiri, Mikael Norrlöf, Johan Löfberg, Anders Hansson.
> *Convex Optimization approach for Time-Optimal Path Tracking of Robots with Speed Dependent Constraints*.
> **报告号**：LiTH-ISY-R-2970，Linköping University，Division of Automatic Control，2010-10-08（投稿 IFAC World Congress 2011）。
>
> 本文档聚焦论文的**核心贡献**：如何把电机的**转矩–转速（torque–speed，简称 t–n）约束**——一种随关节速度变化的力矩上限——**加入时间最优路径跟踪的凸优化里而不破坏凸性**。按用户要求，只梳理论文思路与数学，不涉及代码与软件架构。
>
> ⚠️ **符号提醒**：本论文用 $b=\dot s^2$、$a=\ddot s$（速度平方是 $b$，加速度是 $a$），与本仓库 [`docs/paper_notes.md`](./paper_notes.md)（Wang/Hu 论文，$a=\dot u^2$、$b=\ddot u$）**恰好相反**。本文档一律采用**本论文自身**的记号。

---

## 目录

1. [一句话总览](#1-一句话总览)
2. [背景：Verscheure 的凸形式（基座）](#2-背景verscheure-的凸形式基座)
3. [为什么要 t–n 约束：电机物理动机](#3-为什么要-tn-约束电机物理动机)
4. [核心难点：仿射-于-速度的约束在 $(a,b)$ 参数化下非凸](#4-核心难点仿射-于-速度的约束在-ab-参数化下非凸)
5. [凸化处理：改写到 $(\tau,\dot q^2)$ 平面](#5-凸化处理改写到-tau-dot-q2-平面)
6. [离散化与并入优化问题](#6-离散化与并入优化问题)
7. [扩展后的完整问题（连续 + SOCP）](#7-扩展后的完整问题连续--socp)
8. [扩展：粘滞与库仑摩擦并入 t–n 约束](#8-扩展粘滞与库仑摩擦并入-tn-约束)
9. [数值实验结论](#9-数值实验结论)
10. [方法要点总结](#10-方法要点总结)

---

## 1. 一句话总览

给定关节空间的一条几何路径 $\mathbf q=\mathbf q(s),\ s\in[0,1]$，在满足机器人动力学与执行器约束的前提下求时间参数化 $s(t)$ 使运动时间最短。

- **Verscheure et al. (2009)** 已把该问题在 $b=\dot s^2$、$a=\ddot s$ 的参数化下写成**凸问题（SOCP）**——前提是所有约束在优化变量 $(a,b,\tau)$ 中**仿射**。
- 但真实电机的可用力矩**随转速下降**（反电动势 + 直流母线电压限制），再叠加**粘滞摩擦**，得到一条"转矩–转速"特性曲线（t–n 曲线）。这种**速度相关力矩约束**若按"仿射于关节速度 $\dot q$"直接写入，会因为 $\dot q=\mathbf q'(s)\sqrt{b}$ 中的 $\sqrt b$ 而**破坏凸性**。
- **本文贡献**：证明只要把该约束改写成**仿射于 $(\tau_i,\ \dot q_i^{\,2})$**（转矩与关节速度**平方**）的形式，就能在 $b=\dot s^2$ 的参数化里保持仿射、从而保凸；对于真实（仿射于 $\dot q$ 的）可行域，用一组**仿射不等式做保守内逼近**（切掉非凸的一角）。
- 数值例子：循环时间 9.83 s → 9.38 s（**−4.6%**），顶速提高 **20%**；在最大速度约束活跃的轨迹上，差距可接近 **20%**。

---

## 2. 背景：Verscheure 的凸形式（基座）

$n$ 自由度机械臂动力学（式 1）：

$$
\boldsymbol\tau=\mathbf M(\mathbf q)\ddot{\mathbf q}+\mathbf C(\mathbf q,\dot{\mathbf q})\dot{\mathbf q}+\mathbf F_s(\mathbf q,\dot{\mathbf q})\,\mathrm{sgn}(\dot{\mathbf q})+\mathbf G(\mathbf q).
$$

沿给定路径 $\mathbf q(s)$，令时间进程由标量 $s(t)$（约束 $\dot s\ge0$）承载。引入**两个核心状态量**：

$$
b(s)=\dot s^2\quad(\text{路径速度平方}),\qquad a(s)=\ddot s\quad(\text{路径加速度}).
$$

则关节速度/加速度对 $(a,b)$ **线性**（式 2）：

$$
\dot{\mathbf q}=\mathbf q'(s)\,\dot s=\mathbf q'(s)\sqrt{b(s)},\qquad
\ddot{\mathbf q}=\mathbf q'(s)\,a(s)+\mathbf q''(s)\,b(s).
$$

代回动力学，力矩也对 $(a,b)$ **仿射**：

$$
\boldsymbol\tau(s)=\mathbf m(s)\,a(s)+\mathbf c(s)\,b(s)+\mathbf g(s).\tag{4}
$$

**时间最优凸问题**（式 3–11）：

$$
\min_{a,\,b,\,\tau}\ \int_0^1\frac{1}{\sqrt{b(s)}}\,\mathrm ds
$$

$$
\begin{aligned}
&\boldsymbol\tau(s)=\mathbf m(s)a(s)+\mathbf c(s)b(s)+\mathbf g(s) &&(4)\\
&b(0)=\dot s_0^2,\quad b(1)=\dot s_T^2 &&(5,6)\\
&b'(s)=2\,a(s) &&(7)\quad[\,\tfrac{\mathrm d}{\mathrm ds}\dot s^2=2\ddot s\,]\\
&b(s)\ge0 &&(8)\\
&b(s)\le\bar b(s) &&(9)\\
&\underline{\boldsymbol\tau}(s)\le\boldsymbol\tau(s)\le\bar{\boldsymbol\tau}(s) &&(10)\\
&\underline{\mathbf f}(s)\le\mathbf f(s)a(s)+\mathbf h(s)b(s)\le\bar{\mathbf f}(s) &&(11)
\end{aligned}
$$

**为何是凸的**：目标 $1/\sqrt b$ 对 $b>0$ **凸**；其余约束在 $(a,b,\tau)$ 中全部**仿射**。离散化后可写成**二阶锥规划（SOCP）**，用 SeDuMi / SDPT3 高效求解。

> 关键点（贯穿全文）：**"仿射于优化变量" 是保凸的护身符**。式(9)/(10)/(11) 之所以能进来，正因为它们仿射于 $a,b,\tau$。笛卡尔速度/加速度约束、以及速度盒式约束 $\underline{\dot q}\le \mathbf q'(s)\dot s\le\overline{\dot q}$（式 13）都能等价整理成 $b$ 的上界(9)。t–n 约束的麻烦就在于它**做不到**这一点——见第 4 节。

---

## 3. 为什么要 t–n 约束：电机物理动机

驱动系统有**两条**物理限制（论文 Fig. 1、Fig. 2）：

| 限制 | 来源 | 对 t–n 曲线的影响 |
|------|------|-------------------|
| **热约束** | 电枢电流受发热限制（热耗散能力） | 低速段：**恒定最大转矩**平台 |
| **电压约束** | 直流母线电压有限；**反电动势**（counter-EMF）$\propto$ 角速度 | 高速段：可用转矩随速度**近似线性下降** |

于是执行器真实可用力矩是一条**梯形/折线**特性：低速恒转矩、高速线性收窄。此外还有**粘滞（动）摩擦**：

- 加速时，摩擦**吃掉**一部分可用力矩；
- 减速时，摩擦**帮忙**，可用减速力矩反而**增大**。

所以真实可行域在加速/减速方向**不对称**（减速性能可以更高）。论文例子取对称假设，但方法可直接推广到非对称。

**恒转矩盒式约束的代价**：式(10) 那种与速度无关的常数上下界（"box"），必须把盒子**塞进**真实梯形区域内部（Fig. 2 左）才安全，于是**白白浪费**了高速外扩的能力。改用**速度相关**约束（Fig. 2 右）就能吃满整块梯形——**这就是引入 t–n 约束的收益来源**。

---

## 4. 核心难点：仿射-于-速度的约束在 $(a,b)$ 参数化下非凸

把 t–n 特性写成"仿射于**转矩与关节速度**"的一条约束（式 14）：

$$
\tilde T\,\tau_i(s)+W_i\,\dot q_i(s)\le P_{i1}.
$$

用 $\dot q_i=\mathbf q_i'(s)\,\dot s=\mathbf q_i'(s)\sqrt{b(s)}$ 代入（式 15）：

$$
\tilde T\,\tau_i(s)+W_i\,\mathbf q_i'(s)\,\sqrt{b(s)}\le P_{i1}.\tag{15b}
$$

**问题出在 $\sqrt b$**。这条约束把优化变量 $\tau_i$ 与 $\sqrt{b}$ **耦合**在一起：

$$
g(b,\tau_i)=\tilde T\,\tau_i+W_i\mathbf q_i'\sqrt b .
$$

- $\sqrt b$ 是**凹**函数；当 $W_i\mathbf q_i'>0$ 时 $g$ 整体**凹**，而**凹函数的下水平集 $\{g\le P\}$ 一般非凸**。
- 是否非凸还取决于 $W_i\mathbf q_i'(s)$ 的**符号**（随路径点 $s$ 变化）。

论文原话：*"constraints that are affine in joint velocity are no longer convex in the parameterization where the joint speed square is used"*——**仿射于 $\dot q$ 的约束，一旦搬到用 $\dot s^2$（即 $b$）作变量的参数化里，就不再凸**。因此**不能**把物理 t–n 约束原样塞进 Verscheure 的凸问题。

---

## 5. 凸化处理：改写到 $(\tau,\dot q^2)$ 平面

**核心思想**：不要在 $(\tau_i,\dot q_i)$ 平面描述可行域，改到 $(\tau_i,\ \dot q_i^{\,2})$ 平面——把横轴换成关节速度的**平方**。

### 5.1 为什么换成 $\dot q^2$ 就保凸：线性映射保凸

关键在于 $\dot q_i^2$ 与优化变量 $b$ 是**线性**关系：$\dot q_i^2=\mathbf q_i'^2(s)\,b$，即 $b=\mathbf q_i'^{-2}(s)\,\dot q_i^2$。写成映射（式 16）：

$$
\begin{bmatrix}b(s)\\[2pt]\tau_i\end{bmatrix}
=F\!\left(\dot q_i^2(s),\ \tau_i\right)
=\begin{bmatrix}\mathbf q_i'^{-2}(s)&0\\[2pt]0&1\end{bmatrix}
\begin{bmatrix}\dot q_i^2(s)\\[2pt]\tau_i\end{bmatrix}.
$$

$F$ 是**仿射（这里是线性）映射**。凸分析基本事实（Boyd & Vandenberghe 2004）：**仿射映射把凸集映成凸集**——对任意凸集 $S\subset\mathbb R^2$，$F(S)$ 仍凸。

于是：只要可行域在 $(\dot q_i^2,\tau_i)$ 平面里是**凸集**，它经 $F$ 映到 $(b,\tau_i)$ 平面后**仍是凸集**，并且在真正的优化变量 $b$ 下保持凸——凸性在非线性变换 (12a) 下被**保住**了。

### 5.2 保守仿射内逼近（切掉非凸的一角）

真实可行域在 $(\dot q_i,\tau_i)$ 平面是那块**梯形**（第 3 节，边是仿射的直线）。但横轴从 $\dot q_i$ 换成 $\dot q_i^2$ 后是一次**非线性拉伸**：原来的**直边被弯成曲线**，于是真实区域在 $(\dot q_i^2,\tau_i)$ 平面里出现**非凸**的部分。

**做法**（论文 Fig. 3）：用一组**仿射不等式**在 $(\dot q_i^2,\tau_i)$ 平面里对真实可行域做**保守（内）逼近**——即取真实区域的一个**凸子集**，把弯曲造成的非凸尖角**切掉**（图中阴影 = 逼近后的凸可行集，虚线 = 真实执行器边界）。

- **代价**：被切掉的那一小块力矩/速度能力无法利用（保守）。
- **收益**：可行集**保证凸**，且仍显著大于恒转矩盒子。

第 $j$ 条仿射约束（对第 $i$ 个执行器，式 17，仿射于 $\tau_i$ 与 $\dot q_i^2$）：

$$
T_{ij}\,\tau_i(s)+\bar U_{ij}\,\dot q_i^{\,2}(s)\le P_{ij}.
$$

代入 $\dot q_i^{\,2}=\mathbf q_i'^2(s)\,b$，得到**仿射于优化变量 $(b,\tau)$** 的约束（式 18–19）：

$$
\boxed{\,T_{ij}\,\tau_i(s)+U_{ij}(s)\,b(s)\le P_{ij}\,},\qquad
U_{ij}(s)=\bar U_{ij}\,\mathbf q_i'^{\,2}(s).
$$

**速度相关性去哪了？** 全部吸收进**随路径变化的系数 $U_{ij}(s)=\bar U_{ij}\,\mathbf q_i'^2(s)$**——它只依赖已知的路径几何 $\mathbf q_i'(s)$，可**离线预计算**，对优化器而言只是一条普通的仿射不等式。这正是"速度相关约束"能进凸问题的关键：**把非凸性外包给了几何系数，优化变量里只留仿射**。

---

## 6. 离散化与并入优化问题

在网格 $\{s_k\}$ 上，于**区间中点** $s_{k+1/2}=(s_k+s_{k+1})/2$ 处施加约束，并取 $b_{k+1/2}=(b_k+b_{k+1})/2$（式 20–22）：

$$
T_{ij}\,\tau_i^{\,k}+U_{ij}(s_{k+1/2})\,b_{k+1/2}\le P_{ij},\qquad
U_{ij}(s_{k+1/2})=\bar U_{ij}\,\mathbf q_i'^{\,2}(s_{k+1/2}).
$$

把所有执行器 $i=1..n$ 及其各自的 $m_i$ 条仿射边 $j$ 堆叠，写成矩阵形式（式 23–24）：

$$
\mathbf T\,\boldsymbol\tau^{\,k}+\mathbf U^{\,k+1/2}\,b_{k+1/2}\le \mathbf P,
$$

其中 $\mathbf T$ 是按关节分块的（块对角）转矩系数矩阵，$\mathbf U^{k+1/2}$ 收集各边的 $U_{ij}(s_{k+1/2})$，$\mathbf P$ 收集各边右端 $P_{ij}$。

> 旋转关节（力矩/角速度）与平移关节（力/线速度）**形式完全一致**，方法通用。

---

## 7. 扩展后的完整问题（连续 + SOCP）

**连续形式**（式 25–32）：在 Verscheure 问题上，把原恒转矩盒式约束 (10) 替换/补充为新的速度相关约束 (32)：

$$
\min_{a,b,\tau}\ \int_0^1\frac{\mathrm ds}{\sqrt{b(s)}}
\quad\text{s.t.}\quad
\begin{cases}
\boldsymbol\tau=\mathbf m a+\mathbf c b+\mathbf g &(26)\\
b(0)=\dot s_0^2,\ b(1)=\dot s_T^2 &(27,28)\\
b'(s)=2a(s) &(29)\\
b(s)\ge0 &(30)\\
\underline{\mathbf f}\le \mathbf f a+\mathbf h b\le\bar{\mathbf f} &(31)\\
\mathbf T\boldsymbol\tau(s)+\mathbf U(s)b(s)\le\mathbf P &(32)\ \leftarrow\text{新 t–n 约束}
\end{cases}
$$

**离散 SOCP 形式**（式 33–43）：目标 $1/\sqrt b$ 用辅助变量 $c_k,d_k$ 经**旋转二阶锥**表达。

$$
\min\ \sum_{k=0}^{K-1}2\,\Delta s_k\,d_k\tag{33}
$$

主要约束：

$$
\begin{aligned}
&\boldsymbol\tau^k=\mathbf m(s_{k+1/2})a_k+\mathbf c(s_{k+1/2})b_{k+1/2}+\mathbf g(s_{k+1/2}) &&(34)\\
&b_0=\dot s_0^2,\quad b_K=\dot s_T^2 &&(35,36)\\
&b_{k+1}-b_k=2a_k\Delta s_k &&(37)\ [\,b'=2a\,\text{的离散}\,]\\
&b_k\ge0 &&(38)\\
&\underline{\mathbf f}\le \mathbf f a_k+\mathbf h b_{k+1/2}\le\bar{\mathbf f} &&(39,40)\\
&\big\|\,(2,\ c_k+c_{k+1}-d_k)\,\big\|_2\le c_k+c_{k+1}+d_k &&(41)\\
&\big\|\,(2c_k,\ b_k-1)\,\big\|_2\le b_k+1 &&(42)\\
&\mathbf T\boldsymbol\tau^k+\mathbf U^{k+1/2}b_{k+1/2}\le\mathbf P^{k+1/2} &&(43)\ \leftarrow\text{新 t–n 约束}
\end{aligned}
$$

**两条锥约束的含义**（便于理解目标是怎么"变时间"的）：

- 式(42)：$\|(2c_k,\ b_k-1)\|\le b_k+1 \iff 4c_k^2+(b_k-1)^2\le(b_k+1)^2 \iff c_k^2\le b_k$。即 $c_k\le\sqrt{b_k}=\dot s_k$，$c_k$ 是**路径速度** $\dot s_k$ 的（下）代理。
- 式(41)：$\|(2,\ c_k+c_{k+1}-d_k)\|\le c_k+c_{k+1}+d_k \iff 4\le4(c_k+c_{k+1})d_k \iff d_k\ge\dfrac{1}{c_k+c_{k+1}}$。于是目标 $\sum 2\Delta s_k d_k\approx\sum\dfrac{2\Delta s_k}{\dot s_k+\dot s_{k+1}}$ 正是**区间时长的梯形积分** $\approx\int 1/\dot s\,\mathrm ds$ = 总时间。

新增的式(43) 与原 SOCP 结构相容，**不破坏**任何锥约束，因此整体仍是可高效求解的 SOCP。论文亦指出：该扩展与 Verscheure 的**多目标标量化**框架（在时间最优、热耗能量、力矩变化率之间折中）**不冲突**。

---

## 8. 扩展：粘滞与库仑摩擦并入 t–n 约束

> 论文第 3 节把摩擦**定性**地画进了 Fig. 2 右的非对称梯形；本节把它**定量**落到式(14) 的系数上，说明为什么带摩擦的力矩模型能**原样套用**第 4–5 节的凸化，不需要任何新机制。

**能套用的前提**：破坏凸性的项**只依赖 $\dot q$（速度）、不依赖 $a$（加速度）**——这样它才能整体并入速度相关边界 $\tau_\text{avail}(\dot q)$，再用 $\dot q^2$ 凸化。粘滞摩擦 $F_{v,i}\dot q_i$ 与库仑摩擦 $F_{c,i}\,\mathrm{sgn}(\dot q_i)$ **恰好只依赖 $\dot q$**，与反电动势（counter-EMF）同类，因此都能这么搬。

### 8.1 摩擦系数直接并入斜率 $W_i$ 与截距 $P_i$

保持式(4) 的"无摩擦动力学力矩"仿射于 $(a,b)$：

$$
\tau_i^{\text{dyn}}(s)=m_i(s)\,a+c_i(s)\,b+g_i(s).
$$

电机电磁力矩还需**克服摩擦**：$\tau_i^{\text{dyn}}+F_{v,i}\dot q_i+F_{c,i}\,\mathrm{sgn}(\dot q_i)$；而可用力矩因反电动势随速度线性收窄为 $\tau_{0,i}-\kappa_i|\dot q_i|$。对 $\dot q_i>0$ 的**驱动侧上界**：

$$
\tau_i^{\text{dyn}}+F_{v,i}\dot q_i+F_{c,i}\ \le\ \tau_{0,i}-\kappa_i\dot q_i,
$$

移项即得论文式(14) 的形态 $\tilde T\,\tau_i+W_i\dot q_i\le P_i$：

$$
\underbrace{1}_{\tilde T}\cdot\tau_i^{\text{dyn}}+\underbrace{(\kappa_i+F_{v,i})}_{W_i}\,\dot q_i\ \le\ \underbrace{\tau_{0,i}-F_{c,i}}_{P_i}.
$$

即：**粘滞摩擦并进斜率、库仑摩擦并进截距**——完全落在论文式(14) 的框架内。

| 力矩成分 | 依赖 | 归宿（式 14 系数） |
|----------|------|--------------------|
| 惯量 / 科氏 / 重力 $m_i a+c_i b+g_i$ | $a,b$ | 留在 $\tau_i^{\text{dyn}}$（仿射，进式 4） |
| 反电动势 $\kappa_i\dot q_i$ | $\dot q$ | 斜率 $W_i$ |
| **粘滞摩擦 $F_{v,i}\dot q_i$** | $\dot q$ | **斜率 $W_i=\kappa_i+F_{v,i}$** |
| **库仑摩擦 $F_{c,i}\,\mathrm{sgn}(\dot q_i)$** | $\mathrm{sgn}(\dot q)$ | **截距 $P_i=\tau_{0,i}-F_{c,i}$** |

随后照搬第 5 节：这条仿射于 $\dot q_i$ 的约束在 $b$ 下非凸 → 改写成仿射于 $(\tau_i^{\text{dyn}},\dot q_i^2)$ 并做保守内逼近 → 代入 $\dot q_i^2=\mathbf q_i'^2(s)\,b$ 得仿射约束 $T_{ij}\tau_i+U_{ij}(s)b\le P_{ij}$，$U_{ij}(s)=\bar U_{ij}\mathbf q_i'^2(s)$。

### 8.2 摩擦造成的加速 / 减速不对称（正好对上 Fig. 2 右）

对同一 $\dot q_i>0$，把驱动侧与制动侧写成**两条仿射边**：

$$
\begin{aligned}
\text{驱动侧：}\quad &\tau_i^{\text{dyn}}+(\kappa_i+F_{v,i})\,\dot q_i\le\tau_{0,i}-F_{c,i},\\
\text{制动侧：}\quad &-\tau_i^{\text{dyn}}+(\kappa_i-F_{v,i})\,\dot q_i\le\tau_{0,i}+F_{c,i}.
\end{aligned}
$$

两侧斜率 $\kappa_i\pm F_{v,i}$ **不相等**：因 $F_{v,i}>0$，制动侧随速度衰减更慢（$\kappa_i-F_{v,i}$ 更小），故**高速时可用制动力矩更大**——这正是论文所说"decelerating 时摩擦增大可用力矩"，也正是 Fig. 2 右那块**非对称梯形**。截距上 $\tau_{0,i}\mp F_{c,i}$ 也体现库仑摩擦"低速即让驱动变难、制动变易"。论文对非对称约束只说一句"straightforward to extend"，这里的两条边就是它的具体形态。

### 8.3 库仑摩擦的 $\mathrm{sgn}$ 是**已知**的

$F_{c,i}\,\mathrm{sgn}(\dot q_i)$ 看似含符号非线性，但时间最优下 $\dot s\ge0$，故

$$
\mathrm{sgn}(\dot q_i)=\mathrm{sgn}\!\big(\mathbf q_i'(s)\,\dot s\big)=\mathrm{sgn}\!\big(\mathbf q_i'(s)\big)
$$

**完全由路径几何决定、与优化变量无关**。于是它只是一个**随 $s$ 变的已知常数**并入 $P_i(s)$，不破坏仿射性。仅在 $\mathbf q_i'(s)=0$（关节方向翻转）的孤立点处符号切换，按区间分段取号即可。

### 8.4 一个必须避开的陷阱，以及代价

- **陷阱**：**不要**把 $F_{v,i}\dot q_i$ 留在力矩**等式** (4) 里当作 $\tau$ 的定义。因 $\dot q_i=\mathbf q_i'\sqrt b$，那会让 $\sqrt b$ 进入一条**等式约束**；而凸问题的等式必须仿射，$\sqrt b$ 在等式里**无法用内 / 外逼近救回**（只有不等式能逼近）。正确做法就是 §8.1 的移项：$\tau^{\text{dyn}}$ 保持无摩擦（仿射于 $a,b$），摩擦全部并入**约束边界**。
- **适用边界**：能这么搬的**唯一理由**是摩擦只依赖 $\dot q$。若某力矩项**同时耦合 $a$ 与 $\sqrt b$**（例如某些与加速度相关的非线性损耗），则不在此列，需另行处理。
- **代价**：与论文一致——$\dot q^2$ 坐标下真实边界弯曲，需做**保守仿射内逼近**，切掉非凸尖角，损失一小块能力。

**小结**：粘滞 / 库仑摩擦不但可以、而且是这套方法的**典型适用对象**——粘滞项进斜率 $W_i$、库仑项进截距 $P_i$，非对称由 $\kappa_i\pm F_{v,i}$ 自然给出，其余一切照搬第 4–7 节。

---

## 9. 数值实验结论

基于 Verscheure (2009) 的 6 自由度算例（MATLAB + YALMIP + SeDuMi），关节 4–6 的约束不活跃、只展示关节 1–3。对比两套约束：

| 场景 | 约束 | 说明 |
|------|------|------|
| **Scenario 1** | 速度盒 + 恒转矩盒 | 速度上限设为 Verscheure 例中最大速度的 **90%** |
| **Scenario 2** | 仿射 t–n 约束 | 真实梯形（Fig. 2 右）在 $(\tau_i,\dot q_i^2)$ 平面做仿射保守逼近（Fig. 3） |

结果：

- **循环时间 9.83 s → 9.38 s，减少 4.6%**（Scenario 2 优于 Scenario 1）。
- **顶速提高 20%**（Scenario 2 相对 Scenario 1）；在**最大速度约束活跃**的轨迹上，循环时间差距可接近 **20%**。
- 凸逼近切掉了非凸尖角，**部分额外力矩无法利用**，但净收益仍为正——算法确实吃到了低速多出的转矩与高速多出的速度（Fig. 5、Fig. 6）。

**结论**：速度相关（t–n）约束是刻画执行器真实能力的**通用工具**，只要保证约束**对 $\tau_i$ 与 $\dot q_i^2$ 凸**，就能在充分利用力矩/速度的同时保持整体问题凸、可高效求解。

---

## 10. 方法要点总结

处理 t–n（转矩–转速）约束的**完整思路链**：

1. **物理**：执行器可用力矩随速度下降（热平台 + 反电动势线性收窄）再叠加粘滞摩擦 → 一条梯形 t–n 特性；恒转矩盒式约束太保守。
2. **难点**：把它写成"仿射于 $\dot q$"的约束后，因 $\dot q=\mathbf q'\sqrt b$ 里的 $\sqrt b$ 与 $\tau$ 耦合，在 $b=\dot s^2$ 参数化下**下水平集非凸**（且随 $\mathbf q'(s)$ 符号变化）。
3. **关键变换**：改用**关节速度平方 $\dot q^2$** 作坐标。因 $\dot q_i^2=\mathbf q_i'^2 b$ **线性于 $b$**，线性映射 $F$ 保凸 → 只要可行域**对 $(\tau_i,\dot q_i^2)$ 凸**，就在 $(b,\tau)$ 下保凸。
4. **逼近**：真实梯形的直边在 $\dot q^2$ 坐标下弯曲、局部非凸 → 用一组仿射不等式做**保守内逼近**，切掉非凸尖角（损失一小块能力换凸性）。
5. **落地**：每条边化为 $T_{ij}\tau_i+U_{ij}(s)b\le P_{ij}$，其中 $U_{ij}(s)=\bar U_{ij}\mathbf q_i'^2(s)$ **离线预计算**——速度相关性全部吸收进几何系数，优化变量里只剩仿射。
6. **求解**：并入 Verscheure 的 SOCP（式 43），结构相容、仍高效；对旋转/平移关节通用，且与多目标框架不冲突。

**一句话**：**"仿射于 $\dot q^2$（而非 $\dot q$）"** 是让速度相关力矩约束保凸的钥匙；代价是对真实非凸可行域做一次保守的仿射内逼近。

---

## 参考

- Ardeshiri, Norrlöf, Löfberg, Hansson (2010). *Convex Optimization approach for Time-Optimal Path Tracking of Robots with Speed Dependent Constraints*. LiTH-ISY-R-2970.
- Verscheure, Demeulenaere, Swevers, De Schutter, Diehl (2009). *Time-optimal path tracking for robots: A convex optimization approach*. IEEE TAC, 54(10), 2318–2327. —— 本文的凸形式基座。
- Boyd, Vandenberghe (2004). *Convex Optimization*. —— 仿射映射保凸（$F(S)$ 凸）的依据。
