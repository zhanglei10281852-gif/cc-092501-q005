# 地下水同位素与污染迁移计算服务

这是一个面向水文地质研究团队和环境监管人员的模块化后端，集中管理地下水井、同位素观测、补给端元、混合源反演、污染物迁移、计算任务、参数版本、结果置信区间、用户权限、会话和审计。项目使用 FastAPI 与 SQLite，所有运行数据保存在单个本地数据库文件中，不依赖另行部署的数据库、缓存或消息队列。

## 主要模块

- 井点与样本：登记井点坐标、含水层、采样批次和实验室测量结果。
- 同位素计算：处理稳定同位素、溶质浓度、检测限和质量守恒约束，反演多个补给端元比例。
- 污染迁移：按从源区到监测井的有序含水层区段（长度、孔隙流速、弥散系数、一阶衰减）计算一维平流-弥散-衰减迁移；区段接口质量通量连续，输出目标井突破曲线、峰值、首达时间与累计质量；区段顺序错误、单位不一致或质量误差超限会被拒绝并留痕；相同配置幂等复用，可按参数版本指纹追溯。
- 任务与审计：保存参数版本、计算输入摘要、置信区间、失败重试和结果差异。
- 身份与权限：用户、角色、细粒度权限、会话令牌、账号停用和会话撤销。
- 审计记录：关键身份操作留痕，并对口令和令牌等敏感字段做过滤。
- 后台任务：使用 SQLite 保存待执行任务，支持去重、租约、重试和完成回执。

## 分段污染迁移计算

`POST /api/hydro/wells/{well_id}/transport` 接受两种负载：

- **分段模式（推荐）**：`segments` 为从源区到监测井的有序区段列表，每段含 `length`、`pore_velocity`、`dispersion`、`decay_rate` 与可选 `parameter_version`、`sequence`、`code`；另需 `source_mass`、`duration`、`step`。区段接口上质量通量连续（各段通量型逆高斯核卷积），返回突破曲线 `points`（含逐点 `cumulative_mass`）、`peak`、`first_arrival_days`、`cumulative_mass`、`expected_mass`、`mass_error` 与逐接口 `interface_flux`。
- **单一区段（兼容旧版）**：`source_concentration`、`distance_m`、`velocity_m_day`、`dispersion_m2_day`、`decay_per_day`、`duration_days`、`step_days`，内部视为一个区段。

单位可用 `length_unit`（m/km）与 `time_unit`（day/hour/second）声明，内部统一换算为米·天；区段级单位与运行级不一致会被拒绝。以下情况返回 422 并以 `status=rejected` 留痕（不产出曲线）：区段顺序错误（`segments_out_of_order`）、单位不一致（`inconsistent_units`）、质量误差超门限（`mass_error_exceeded`，响应中附建议的模拟时长 `suggested_duration`）。相同配置重复提交返回同一记录（幂等）；`GET /api/hydro/transport/runs/{id}` 与 `GET /api/hydro/transport/runs?well_id=` 可追溯历史运行采用的输入、求解器版本与参数版本指纹。

## 运行环境

- Python 3.11
- SQLite 3，由 Python 标准库提供
- Linux、macOS 或 Windows

## 安装

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

配置项均以 `TOWNSHIP_` 开头。可以复制 `.env.example` 后按需设置，默认数据库位于 `./data/township.db`。

## 初始化与检查

```bash
python -m app.cli init-db
python -m app.cli check-db
```

## 启动服务

```bash
uvicorn app.main:app --host 0.0.0.0 --port 8432
```

健康检查：

```bash
curl -sS http://127.0.0.1:8432/api/system/health
```

首次部署可创建唯一的初始管理员：

```bash
curl -sS -X POST http://127.0.0.1:8432/api/auth/bootstrap   -H 'Content-Type: application/json'   -d '{"username":"admin","password":"Admin!23456","client_label":"initial-setup"}'
```

之后通过 `/api/auth/login` 获取会话令牌，并在管理接口请求头中使用 `Authorization: Bearer <token>`。

## 测试

```bash
python -m pytest
```

测试覆盖身份初始化、登录、用户与角色维护、权限计算、账号停用后的会话撤销、审计脱敏、井点样本、同位素约束、迁移计算、任务恢复和数据库时间格式。

## 编译检查

```bash
python -m compileall -q app tests
```

## API 冒烟

```bash
python -m app.cli smoke
```

该命令在进程内启动应用并检查服务根路径与健康接口，适合部署前快速确认路由和数据库初始化是否正常。

## 目录结构

```text
app/
  api/             用户、角色、审计、认证和系统接口
  core/            时钟、安全、异常和分页能力
  repositories/    SQLite 查询与持久化读取
  routers/         灾情、事件、公告、部门和信访业务接口
  hydro/            地下水、同位素反演和污染迁移服务
  schemas/         管理接口输入模型
  services/        身份、审计和后台任务领域服务
  cli.py           初始化、检查和冒烟入口
  database.py      SQLite 连接、事务、表结构与基础权限
tests/             核心、管理接口和原有业务回归测试
tools/             本地维护脚本
```

## 数据一致性

SQLite 连接默认启用外键、WAL、busy timeout 与同步写入策略。需要跨多张表更新的管理操作在即时事务中执行，失败会整体回滚。会话令牌只保存摘要；用户停用会撤销仍有效的会话。审计事件保存操作者、动作、资源、结果和前后状态，但不会保存明文密码或令牌。
