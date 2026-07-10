# CLAUDE.md

## 交流语言

请始终用中文回复。

## 运行测试

本项目依赖（cvxpy / osqp / qdldl 等）的最新 Windows 轮子会段错误崩溃，已固定到稳定版本并安装在独立虚拟环境 `.venv` 中。请用该环境运行：

```powershell
.\.venv\Scripts\python.exe -m pytest tests/ -q
```

注意：不要用 Anaconda base 环境安装本项目依赖（会把 numpy 升到 2.x，破坏 base 的 scipy/matplotlib）。
