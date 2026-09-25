from __future__ import annotations

import json
import math


def create_well(client, code="W-100"):
    response = client.post("/api/hydro/wells", json={
        "code": code, "name": "污染羽监测井", "latitude": 35.1, "longitude": 116.2,
        "aquifer": "多层孔隙含水层", "screen_depth_m": 42,
    })
    assert response.status_code == 201, response.text
    return response.json()


SEGMENTS = [
    {"sequence": 1, "code": "A-源区强透水层", "length": 50.0, "pore_velocity": 1.0,
     "dispersion": 2.0, "decay_rate": 0.01, "parameter_version": "k-zone-3"},
    {"sequence": 2, "code": "B-弱透水层", "length": 80.0, "pore_velocity": 0.5,
     "dispersion": 8.0, "decay_rate": 0.0, "parameter_version": "k-zone-3"},
    {"sequence": 3, "code": "C-井旁层", "length": 30.0, "pore_velocity": 3.0,
     "dispersion": 1.0, "decay_rate": 0.005, "parameter_version": "k-zone-4"},
]


def _run(client, well_id, body):
    return client.post(f"/api/hydro/wells/{well_id}/transport", json=body)


def test_segmented_transport_breakthrough_curve_and_indicators(client):
    well = create_well(client)
    response = _run(client, well["id"], {
        "source_mass": 10.0, "duration": 900.0, "step": 5.0,
        "model_version": "ade-seg-test", "segments": SEGMENTS,
    })
    assert response.status_code == 201, response.text
    run = response.json()
    assert run["status"] == "done"
    assert run["segment_count"] == 3
    assert run["solver_version"] == "ade-seg-1"
    result = json.loads(run["result_json"])
    points = result["points"]

    # 突破曲线单调上升后下降：存在唯一上升沿与完整回落
    concentrations = [p["concentration"] for p in points]
    peak_index = max(range(len(concentrations)), key=concentrations.__getitem__)
    assert 0 < peak_index < len(concentrations) - 2
    assert concentrations[0] < concentrations[peak_index]
    assert concentrations[-1] < concentrations[peak_index]

    # 累计质量单调不增
    cumulative = [p["cumulative_mass"] for p in points]
    assert all(b >= a - 1e-9 for a, b in zip(cumulative, cumulative[1:]))

    # 峰值位于报告点之间插值，且与报告网格峰值一致量级
    assert 100.0 < result["peak"]["time_days"] < 220.0
    assert abs(result["peak"]["concentration"] - max(concentrations)) < max(concentrations) * 0.02

    # 首达时间早于峰值
    assert result["first_arrival_days"] is not None
    assert result["first_arrival_days"] < result["peak"]["time_days"]
    # 平均到达时间 = Σ L/v = 50 + 160 + 10 = 220 天
    assert math.isclose(result["arrival_time_days"], 220.0, rel_tol=1e-12)

    # 累计质量 = 源质量 × 各段存活因子，数值积分误差在门限内
    assert math.isclose(result["cumulative_mass"], result["expected_mass"], rel_tol=1e-2)
    assert result["mass_error"] < result["mass_error_tolerance"]
    assert result["cumulative_mass"] < 10.0  # 衰减使到达质量小于源质量

    # 区段与参数版本完整回传，可追溯
    assert [s["parameter_version"] for s in result["segments"]] == ["k-zone-3", "k-zone-3", "k-zone-4"]
    assert result["parameter_fingerprint"]
    # 接口质量通量表：源区 → 每段出口
    interfaces = result["interface_flux"]
    assert [row["interface"] for row in interfaces] == [0, 1, 2, 3]
    assert math.isclose(interfaces[0]["cumulative_mass"], 10.0)
    assert interfaces[1]["cumulative_mass"] < interfaces[0]["cumulative_mass"]  # 第一段有衰减
    assert math.isclose(interfaces[2]["cumulative_mass"], interfaces[1]["cumulative_mass"])  # 第二段无衰减
    assert interfaces[3]["cumulative_mass"] < interfaces[2]["cumulative_mass"]  # 第三段有衰减
    assert math.isclose(interfaces[-1]["cumulative_mass"], result["expected_mass"], rel_tol=1e-12)


def test_mass_flux_continuity_matches_analytic_single_segment(client):
    """单段无衰减时，分段卷积必须与通量型解析解逐点一致（接口质量通量连续的特例）。"""
    well = create_well(client, "W-101")
    response = _run(client, well["id"], {
        "source_mass": 1.0, "duration": 400.0, "step": 2.0,
        "segments": [{"length": 100.0, "pore_velocity": 2.0, "dispersion": 5.0}],
    })
    assert response.status_code == 201, response.text
    result = json.loads(response.json()["result_json"])

    def analytic(t):
        return 100.0 / math.sqrt(4 * math.pi * 5.0 * t ** 3) * math.exp(-(100.0 - 2.0 * t) ** 2 / (4 * 5.0 * t))

    for point in result["points"]:
        assert math.isclose(point["concentration"], analytic(point["time_days"]), rel_tol=1e-8, abs_tol=1e-12)
    assert math.isclose(result["cumulative_mass"], 1.0, rel_tol=1e-9)
    assert math.isclose(result["expected_mass"], 1.0, rel_tol=1e-12)


def test_segment_order_errors_are_rejected(client):
    well = create_well(client, "W-102")
    bad = [
        {"sequence": 2, "length": 10.0, "pore_velocity": 1.0, "dispersion": 1.0},
        {"sequence": 1, "length": 10.0, "pore_velocity": 1.0, "dispersion": 1.0},
    ]
    response = _run(client, well["id"], {"source_mass": 1.0, "duration": 100.0, "segments": bad})
    assert response.status_code == 422
    assert response.json()["detail"]["code"] == "segments_out_of_order"

    # 拒绝结果被留痕，且同一坏配置重复提交只保留一条 rejected 记录
    again = _run(client, well["id"], {"source_mass": 1.0, "duration": 100.0, "segments": bad})
    assert again.status_code == 422
    runs = client.get(f"/api/hydro/transport/runs?well_id={well['id']}").json()["items"]
    rejected = [r for r in runs if r["status"] == "rejected"]
    assert len(rejected) == 1
    stored = json.loads(rejected[0]["result_json"])
    assert stored["rejection"] == "segments_out_of_order"


def test_inconsistent_units_are_rejected(client):
    well = create_well(client, "W-103")
    response = _run(client, well["id"], {
        "source_mass": 1.0, "duration": 100.0, "length_unit": "m",
        "segments": [{"length": 10.0, "pore_velocity": 1.0, "dispersion": 1.0, "length_unit": "km"}],
    })
    assert response.status_code == 422
    assert response.json()["detail"]["code"] == "inconsistent_units"


def test_unit_conversion_is_equivalent(client):
    well = create_well(client, "W-104")
    body_m = {"source_mass": 1.0, "duration": 100.0, "step": 1.0,
              "segments": [{"length": 100.0, "pore_velocity": 100.0, "dispersion": 10.0}]}
    body_km = {"source_mass": 1.0, "duration": 100.0, "step": 1.0, "length_unit": "km",
               "segments": [{"length": 0.1, "pore_velocity": 0.1, "dispersion": 1e-5}]}
    run_m = _run(client, well["id"], body_m).json()
    run_km = _run(client, well["id"], body_km).json()
    # 换算后物理参数完全一致 → 归并为同一次运行（单位无关的幂等键）
    assert run_m["id"] == run_km["id"]
    result = json.loads(run_m["result_json"])
    # 内部统一为 SI（米、天）
    assert result["unit_system"]["length"] == "m"
    assert result["unit_system"]["time"] == "day"

    # 先提交 km 的井上，requested_length 如实记录首次采用的单位
    well2 = create_well(client, "W-104B")
    run_km_first = _run(client, well2["id"], body_km).json()
    assert json.loads(run_km_first["result_json"])["unit_system"]["requested_length"] == "km"


def test_mass_error_exceeded_rejects_and_suggests_duration(client):
    well = create_well(client, "W-105")
    response = _run(client, well["id"], {
        "source_mass": 10.0, "duration": 120.0, "step": 5.0, "segments": SEGMENTS,
    })
    assert response.status_code == 422
    detail = response.json()["detail"]
    assert detail["code"] == "mass_error_exceeded"
    assert detail["context"]["mass_error"] > detail["context"]["tolerance"]
    assert detail["context"]["suggested_duration"] > 120.0

    # 按建议延长时长后结果被接受
    fixed = _run(client, well["id"], {
        "source_mass": 10.0, "duration": 1040.0, "step": 5.0, "segments": SEGMENTS,
    })
    assert fixed.status_code == 201, fixed.text
    result = json.loads(fixed.json()["result_json"])
    assert result["mass_error"] < result["mass_error_tolerance"]


def test_identical_configuration_is_idempotent_and_traceable(client):
    well = create_well(client, "W-106")
    body = {"source_mass": 10.0, "duration": 900.0, "step": 5.0,
            "model_version": "ade-seg-test", "segments": SEGMENTS}
    first = _run(client, well["id"], body)
    second = _run(client, well["id"], body)
    assert first.status_code == second.status_code == 201
    assert first.json()["id"] == second.json()["id"]
    assert first.json()["task_key"] == second.json()["task_key"]

    # 修改某一段的参数版本 → 视为新配置，产生新的运行记录与不同指纹
    changed = json.loads(json.dumps(SEGMENTS))
    changed[2]["parameter_version"] = "k-zone-5"
    third = _run(client, well["id"], {**body, "segments": changed})
    assert third.status_code == 201
    assert third.json()["id"] != first.json()["id"]
    r_old = json.loads(first.json()["result_json"])
    r_new = json.loads(third.json()["result_json"])
    assert r_old["parameter_fingerprint"] != r_new["parameter_fingerprint"]

    # GET 可追溯任一历史运行采用的输入与参数版本
    fetched = client.get(f"/api/hydro/transport/runs/{first.json()['id']}")
    assert fetched.status_code == 200
    stored_input = json.loads(fetched.json()["input_json"])
    assert stored_input["segments"][0]["parameter_version"] == "k-zone-3"
    assert fetched.json()["model_version"] == "ade-seg-test"


def test_legacy_single_segment_payload_still_supported(client):
    well = create_well(client, "W-107")
    response = _run(client, well["id"], {
        "source_concentration": 100, "distance_m": 100, "velocity_m_day": 2,
        "dispersion_m2_day": 5, "decay_per_day": 0.01, "duration_days": 100,
        "step_days": 5, "model_version": "ade-test",
    })
    assert response.status_code == 201, response.text
    run = response.json()
    assert run["status"] == "done"
    assert run["segment_count"] == 1
    result = json.loads(run["result_json"])
    assert result["peak"]["concentration"] > 0
    assert result["first_arrival_days"] is not None
    assert result["cumulative_mass"] > 0
