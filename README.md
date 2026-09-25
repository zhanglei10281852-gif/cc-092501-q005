# 地下水同位素与污染迁移计算服务

这是一个面向水文地质研究团队和环境监管人员的模块化后端，集中管理地下水井、同位素观测、补给端元、混合源反演、污染物迁移、计算任务、参数版本、结果置信区间、用户权限、会话和审计。项目使用 FastAPI 与 SQLite，所有运行数据保存在单个本地数据库文件中，不依赖另行部署的数据库、缓存或消息队列。

## 主要模块

- 井点与样本：登记井点坐标、含水层、采样批次和实验室测量结果。
- 同位素计算：处理稳定同位素、溶质浓度、检测限和质量守恒约束，反演多个补给端元比例。
- 污染迁移：一维平流、弥散和一阶衰减。单段计算给出到达时间和浓度曲线；分段计算接受从源区到监测井的**有序含水层区段**(每段含长度、孔隙流速、弥散系数和一阶衰减),以区段脉冲响应串联卷积保证**区段接口质量通量连续**,输出目标井突破曲线、峰值、首达时间与累计质量,并对区段顺序、单位一致性和质量守恒误差进行校验,不合格即拒绝结果。
- 任务与审计：保存模型与参数版本、计算输入摘要、突破曲线、置信/质量守恒指标、失败重试和结果差异。相同配置重复运行复用同一条记录,可追溯采用的参数版本。
- 身份与权限：用户、角色、细粒度权限、会话令牌、账号停用和会话撤销。
- 审计记录：关键身份操作留痕，并对口令和令牌等敏感字段做过滤。
- 后台任务：使用 SQLite 保存待执行任务，支持去重、租约、重试和完成回执。

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

## 分段污染迁移

污染羽从源区到监测井会穿过渗透系数与衰减条件不同的含水层区段。分段迁移接口
`POST /api/hydro/wells/{well_id}/segment-transport` 接受按源区 → 监测井方向排序的
区段序列,每段给出长度、孔隙流速、弥散系数和一阶衰减;可选用 `start_m`/`end_m`
里程坐标显式锁定区段接口衔接与方向。

```json
{
  "source_mass_kg": 1000.0,
  "distance_m": 100.0,
  "segments": [
    {"length_m": 50, "velocity_m_day": 1.0, "dispersion_m2_day": 2.0, "decay_per_day": 0.0,
     "start_m": 0, "end_m": 50},
    {"length_m": 30, "velocity_m_day": 0.6, "dispersion_m2_day": 5.0, "decay_per_day": 0.01,
     "start_m": 50, "end_m": 80},
    {"length_m": 20, "velocity_m_day": 2.0, "dispersion_m2_day": 1.0, "decay_per_day": 0.0,
     "start_m": 80, "end_m": 100}
  ],
  "duration_days": 1200.0,
  "step_days": 1.0,
  "first_arrival_quantile": 0.01,
  "mass_tolerance": 0.02,
  "model_version": "ade-segment-1",
  "parameter_version": "param-2026q3",
  "units": {"length": "m", "time": "day", "mass": "kg",
            "velocity": "m/day", "dispersion": "m2/day", "decay": "1/day"}
}
```

返回记录中 `result_json` 含突破曲线 `points`、`peak`、`first_arrival_time_days`、
`cumulative_mass_kg` 与解析期望质量、质量误差、各段存活率和参数版本。求解采用区段
脉冲响应(逆高斯首达密度 × 一阶衰减存活因子)逐级卷积,上段出流质量通量即下段入流,
接口质量通量连续。以下情况返回 422 拒绝:

- 区段长度合计与 `distance_m` 不闭合,或里程坐标未衔接/方向倒退(顺序错误);
- `units` 与规范单位制(m、day、kg、m/day、m²/day、1/day)不一致或缺漏;
- 时间窗口未覆盖完整突破过程,数值累计质量与解析期望质量偏差超过 `mass_tolerance`。

相同配置(含 `parameter_version`)哈希为同一 `task_key`,重复运行返回同一记录;
`GET /api/hydro/transport-runs/{id}` 可追溯完整输入、模型版本与参数版本。

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
