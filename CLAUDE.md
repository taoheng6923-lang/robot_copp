# CLAUDE.md

## 交流语言

请始终用中文回复。

## 运行测试

本项目依赖（cvxpy / osqp / qdldl 等）的最新 Windows 轮子会段错误崩溃，已固定到稳定版本并安装在独立虚拟环境 `.venv` 中。请用该环境运行（`pyproject.toml` 已配置 `testpaths`/`pythonpath`，不用带路径参数）：

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

测试代码本体在 `trajectory-planning/copp/self-test/`（内核）与
`trajectory-planning/path/self-test/`（指令+降维），可视化产物落各自
`self-test/output/`（gitignored）。

注意：不要用 Anaconda base 环境安装本项目依赖（会把 numpy 升到 2.x，破坏 base 的 scipy/matplotlib）。

## 项目结构速览

- 顶层 `robot/`：机器人本体（`ur5.py` 真实 UR5 DH 运动学 + 解析 IK），独立于 `trajectory-planning/`
- `trajectory-planning/copp/`：TOTP-SPLP 数值求解核心（已实现，M1+M4）
- `trajectory-planning/path/`：路径构造（M2 已实现 `commands/`+`lowering/`；`blending/` 属 M3 占位）
- `trajectory-planning/planner/`：调度门面（M2+/M5，目前只有占位 `__init__.py`）

详细的模块状态、里程碑进度见 [`docs/README_M1.md`](docs/README_M1.md)（内核）与
[`docs/README_M2.md`](docs/README_M2.md)（指令+降维）——这两份是本仓库"如实反映当前代码"
的权威文档；`docs/` 下其余几份是设计/算法参考，可能与当前代码结构有出入，各自开头都有说明。
