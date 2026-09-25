from __future__ import annotations

import json
import math

UNITS = {
    "length": "m",
    "time": "day",
    "mass": "kg",
    "velocity": "m/day",
    "dispersion": "m2/day",
    "decay": "1/day",
}


def create_well(client, code="W-SEG"):
    response = client.post(
        "/api/hydro/wells",
        json={
            "code": code,
            "name": "下游监测井",
            "latitude": 35.1,
            "longitude": 116.2,
            "aquifer": "多层孔隙含水层",
            "screen_depth_m": 60,
        },
    )
    assert response.status_code == 201, response.text
    return response.json()


def segment_payload(**overrides):
    payload = {
        "source_mass_kg": 1000.0,
        "distance_m": 100.0,
        "segments": [
            {"length_m": 50, "velocity_m_day": 1.0, "dispersion_m2_day": 2.0, "decay_per_day": 0.0},
            {"length_m": 30, "velocity_m_day": 0.6, "dispersion_m2_day": 5.0, "decay_per_day": 0.01},
            {"length_m": 20, "velocity_m_day": 2.0, "dispersion_m2_day": 1.0, "decay_per_day": 0.0},
        ],
        "duration_days": 1200.0,
        "step_days": 1.0,
        "first_arrival_quantile": 0.01,
        "mass_tolerance": 0.02,
        "model_version": "ade-segment-1",
        "parameter_version": "param-2026q3",
        "units": dict(UNITS),
    }
    payload.update(overrides)
    return payload


def post_segments(client, well_id, payload):
    return client.post(f"/api/hydro/wells/{well_id}/segment-transport", json=payload)


def test_segment_transport_breakthrough_curve(client):
    well = create_well(client)
    response = post_segments(client, well["id"], segment_payload())
    assert response.status_code == 201, response.text
    record = response.json()
    assert record["status"] == "done"
    assert record["run_type"] == "segment"
    assert record["parameter_version"] == "param-2026q3"

    result = json.loads(record["result_json"])
    # 突破曲线:时间单调递增、通量非负、存在明显峰
    points = result["points"]
    assert len(points) == 1200
    times = [p["time_days"] for p in points]
    assert times == sorted(times)
    assert all(p["mass_flux_kg_day"] >= 0 for p in points)

    peak = result["peak"]
    assert peak["time_days"] is not None and peak["mass_flux_kg_day"] > 0
    # 峰值应出现在平流到达时间(110 天)附近而非窗口边缘
    assert 30 < peak["time_days"] < 300
    # 首达时间早于峰值时间
    assert 0 < result["first_arrival_time_days"] < peak["time_days"]

    # 累计质量 = 源质量 × 各段解析存活率之积(仅中段有衰减)
    lam, length, v, d = 0.01, 30.0, 0.6, 5.0
    survival = math.exp(-2 * lam * length / (v + math.sqrt(v * v + 4 * d * lam)))
    expected = 1000.0 * survival
    assert result["expected_mass_kg"] == expected
    assert abs(result["cumulative_mass_kg"] - expected) < 1e-6 * expected
    assert result["mass_error"] <= 0.02
    assert result["n_segments"] == 3
    assert len(result["segment_summary"]) == 3
    assert result["model_version"] == "ade-segment-1"
    assert result["parameter_version"] == "param-2026q3"


def test_segment_transport_conserves_mass_without_decay(client):
    well = create_well(client, "W-SEG-CONS")
    payload = segment_payload(
        segments=[
            {"length_m": 50, "velocity_m_day": 1.0, "dispersion_m2_day": 2.0, "decay_per_day": 0.0},
            {"length_m": 30, "velocity_m_day": 0.6, "dispersion_m2_day": 5.0, "decay_per_day": 0.0},
            {"length_m": 20, "velocity_m_day": 2.0, "dispersion_m2_day": 1.0, "decay_per_day": 0.0},
        ]
    )
    response = post_segments(client, well["id"], payload)
    assert response.status_code == 201, response.text
    result = json.loads(response.json()["result_json"])
    # 区段接口质量通量连续:无衰减时全部源质量最终到达监测井
    assert abs(result["cumulative_mass_kg"] - 1000.0) < 1e-3
    assert result["survival_fraction"] == 1.0


def test_segment_transport_same_config_is_idempotent_and_traceable(client):
    well = create_well(client, "W-SEG-IDEM")
    payload = segment_payload()
    first = post_segments(client, well["id"], payload)
    assert first.status_code == 201, first.text
    second = post_segments(client, well["id"], payload)
    assert second.status_code == 201, second.text
    # 相同配置复用同一条计算记录,结果逐字节一致(可重复运行)
    assert first.json()["id"] == second.json()["id"]
    assert first.json()["result_json"] == second.json()["result_json"]
    assert first.json()["task_key"] == second.json()["task_key"]

    # 参数版本可追溯:按记录 id 查询,模型版本、参数版本与完整输入均在案
    trace = client.get(f"/api/hydro/transport-runs/{first.json()['id']}")
    assert trace.status_code == 200, trace.text
    record = trace.json()
    assert record["model_version"] == "ade-segment-1"
    assert record["parameter_version"] == "param-2026q3"
    stored_input = json.loads(record["input_json"])
    for stored, sent in zip(stored_input["segments"], payload["segments"]):
        for key in ("length_m", "velocity_m_day", "dispersion_m2_day", "decay_per_day"):
            assert stored[key] == sent[key]
    assert stored_input["source_mass_kg"] == 1000.0


def test_segment_transport_different_parameter_version_reruns(client):
    well = create_well(client, "W-SEG-VER")
    first = post_segments(client, well["id"], segment_payload())
    assert first.status_code == 201
    # 参数版本不同视为新配置,重新计算而非复用
    second = post_segments(client, well["id"], segment_payload(parameter_version="param-2026q4"))
    assert second.status_code == 201
    assert second.json()["id"] != first.json()["id"]
    assert second.json()["parameter_version"] == "param-2026q4"


def test_segment_transport_rejects_unclosed_segment_chain(client):
    well = create_well(client, "W-SEG-ORDER")
    # 区段长度合计 120 m,与源-井距离 100 m 不闭合:顺序/长度错误
    payload = segment_payload(
        segments=[
            {"length_m": 60, "velocity_m_day": 1.0, "dispersion_m2_day": 2.0, "decay_per_day": 0.0},
            {"length_m": 60, "velocity_m_day": 1.0, "dispersion_m2_day": 2.0, "decay_per_day": 0.0},
        ]
    )
    response = post_segments(client, well["id"], payload)
    assert response.status_code == 422, response.text


def test_segment_transport_rejects_unit_mismatch(client):
    well = create_well(client, "W-SEG-UNIT")
    # 时间单位声明为 hour,与规范单位 day 不一致
    bad_units = dict(UNITS)
    bad_units["time"] = "hour"
    response = post_segments(client, well["id"], segment_payload(units=bad_units))
    assert response.status_code == 422, response.text
    # 未显式声明单位同样拒绝
    missing = segment_payload()
    del missing["units"]
    response = post_segments(client, well["id"], missing)
    assert response.status_code == 422, response.text


def test_segment_transport_rejects_excess_mass_error(client):
    well = create_well(client, "W-SEG-MASS")
    # 模拟窗口只有 30 天,远小于 110 天平流到达时间,大部分质量未到达
    payload = segment_payload(duration_days=30.0, step_days=0.5)
    response = post_segments(client, well["id"], payload)
    assert response.status_code == 422, response.text
    assert "质量误差" in response.json()["detail"]

    # 相同配置重复提交:仍是同一条拒绝记录,不会写成成功结果
    again = post_segments(client, well["id"], payload)
    assert again.status_code == 422, again.text

    # 延长窗口覆盖完整突破过程后,同一源参数计算成功
    fixed = segment_payload(duration_days=1200.0, step_days=1.0)
    ok = post_segments(client, well["id"], fixed)
    assert ok.status_code == 201, ok.text
    assert ok.json()["status"] == "done"


def test_segment_transport_rejects_oversized_grid(client):
    well = create_well(client, "W-SEG-GRID")
    payload = segment_payload(duration_days=100000.0, step_days=1.0)
    response = post_segments(client, well["id"], payload)
    assert response.status_code == 422, response.text


def test_segment_transport_rejects_wrong_segment_order_by_chainage(client):
    well = create_well(client, "W-SEG-CHAIN")
    # 两段等长,仅靠长度合计无法发现交换;里程坐标暴露出接口未衔接
    segments = [
        {"length_m": 50, "velocity_m_day": 1.0, "dispersion_m2_day": 2.0,
         "decay_per_day": 0.0, "start_m": 50, "end_m": 100},
        {"length_m": 50, "velocity_m_day": 0.6, "dispersion_m2_day": 5.0,
         "decay_per_day": 0.0, "start_m": 0, "end_m": 50},
    ]
    response = post_segments(client, well["id"], segment_payload(segments=segments))
    assert response.status_code == 422, response.text
    assert "顺序" in response.text

    # 里程倒退(终点 <= 起点)同样拒绝
    reversed_one = [dict(segments[1]), dict(segments[0])]
    reversed_one[0]["start_m"], reversed_one[0]["end_m"] = 50, 0
    response = post_segments(client, well["id"], segment_payload(segments=reversed_one))
    assert response.status_code == 422, response.text


def test_segment_transport_accepts_consistent_chainage(client):
    well = create_well(client, "W-SEG-CHAIN-OK")
    base = segment_payload()
    chainage = [0, 50, 80, 100]
    for index, segment in enumerate(base["segments"]):
        segment["start_m"] = chainage[index]
        segment["end_m"] = chainage[index + 1]
    response = post_segments(client, well["id"], base)
    assert response.status_code == 201, response.text
    assert response.json()["status"] == "done"


def test_segment_transport_rejects_partial_chainage(client):
    well = create_well(client, "W-SEG-CHAIN-PART")
    payload = segment_payload()
    payload["segments"][0]["start_m"] = 0
    payload["segments"][0]["end_m"] = 50
    response = post_segments(client, well["id"], payload)
    assert response.status_code == 422, response.text


def test_segment_transport_rejects_grid_too_coarse(client):
    well = create_well(client, "W-SEG-COARSE")
    # 平流到达仅 1 天,但时间步长 100 天,脉冲响应采样不到任何质量
    payload = segment_payload(
        segments=[
            {"length_m": 100, "velocity_m_day": 100.0, "dispersion_m2_day": 5.0, "decay_per_day": 0.0}
        ],
        duration_days=1000.0,
        step_days=100.0,
    )
    response = post_segments(client, well["id"], payload)
    assert response.status_code == 422, response.text


def test_segment_transport_well_not_found(client):
    response = post_segments(client, 9999, segment_payload())
    assert response.status_code == 404, response.text


def test_segment_transport_single_segment_matches_analytical(client):
    well = create_well(client, "W-SEG-SINGLE")
    payload = segment_payload(
        distance_m=100.0,
        segments=[
            {"length_m": 100, "velocity_m_day": 1.0, "dispersion_m2_day": 5.0, "decay_per_day": 0.01}
        ],
        duration_days=400.0,
        step_days=0.5,
    )
    response = post_segments(client, well["id"], payload)
    assert response.status_code == 201, response.text
    result = json.loads(response.json()["result_json"])
    # 单段解析存活率 exp(-2λL/(v+sqrt(v^2+4Dλ)))
    survival = math.exp(-2 * 0.01 * 100 / (1 + math.sqrt(1 + 4 * 5 * 0.01)))
    assert abs(result["cumulative_mass_kg"] - 1000.0 * survival) < 1e-6 * 1000.0
    # 平流到达时间 100 天;首达(1% 分位)应显著早于它
    assert result["advective_travel_time_days"] == 100.0
    assert result["first_arrival_time_days"] < 100.0
