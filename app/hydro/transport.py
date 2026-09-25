"""一维分段溶质迁移计算引擎。

污染羽从源区到目标井依次穿过若干含水层区段，每段具有独立的长度 L、
孔隙流速 v、弥散系数 D 与一阶衰减 λ。区段之间以质量通量连续为约束：
第 i+1 段入口的通量浓度等于第 i 段出口的通量浓度，因此整条路径的
脉冲响应核是各段通量核（带一阶衰减的逆高斯分布）的卷积：

    g(t) = g_1 * g_2 * ... * g_N
    g_i(t; L_i, v_i, D_i, λ_i)
        = L_i / sqrt(4 π D_i t^3)
          * exp(-(L_i - v_i t)^2 / (4 D_i t) - λ_i t)

每段质量存活因子（解析值，来自拉氏变换 s → s+λ）：

    S_i = exp[v_i L_i/(2 D_i) * (1 - sqrt(1 + 4 D_i λ_i / v_i^2))]

无衰减时 S_i = 1，∫g_i dt = 1，接口处质量通量严格连续。
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any

# 支持的输入单位 → 内部规范单位（米、天）
LENGTH_TO_M: dict[str, float] = {"m": 1.0, "km": 1000.0}
TIME_TO_DAY: dict[str, float] = {"day": 1.0, "hour": 1.0 / 24.0, "second": 1.0 / 86400.0}

SOLVER_VERSION = "ade-seg-1"
_MAX_REPORT_STEPS = 12_000
_MAX_INTERNAL_STEPS = 24_000
_KERNEL_CUTOFF = 1e-13


class TransportError(ValueError):
    """迁移参数或质量核算未通过，结果必须拒绝。"""

    def __init__(self, code: str, message: str, context: dict[str, Any] | None = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.context = context or {}


@dataclass(frozen=True)
class Segment:
    sequence: int
    code: str | None
    length_m: float
    velocity_m_day: float
    dispersion_m2_day: float
    decay_per_day: float
    parameter_version: str


@dataclass(frozen=True)
class TransportConfig:
    source_mass: float
    duration_days: float
    step_days: float
    segments: tuple[Segment, ...]
    detection_limit: float
    relative_threshold: float
    mass_error_tolerance: float
    model_version: str
    length_unit: str
    time_unit: str
    time_factor: float  # 请求时间单位 → 天，用于把建议时长换算回用户单位


def _finite_positive(name: str, value: Any) -> float:
    if value is None:
        raise TransportError("invalid_transport_config", f"缺少参数 {name}")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise TransportError("invalid_transport_config", f"参数 {name} 不是数值") from exc
    if not math.isfinite(number) or number <= 0:
        raise TransportError("invalid_transport_config", f"参数 {name} 必须为正数")
    return number


def _build_segment(raw: dict[str, Any], index: int, lf: float, tf: float,
                   length_unit: str, time_unit: str,
                   used_codes: set[str]) -> Segment:
    sequence = raw.get("sequence", index)
    if not isinstance(sequence, int) or sequence < 1:
        raise TransportError("segments_out_of_order", f"第 {index} 段序号无效，必须从 1 开始递增")
    seg_length_unit = raw.get("length_unit")
    seg_time_unit = raw.get("time_unit")
    if seg_length_unit is not None and seg_length_unit != length_unit:
        raise TransportError(
            "inconsistent_units",
            f"第 {sequence} 段长度单位 {seg_length_unit!r} 与运行级单位 {length_unit!r} 不一致",
        )
    if seg_time_unit is not None and seg_time_unit != time_unit:
        raise TransportError(
            "inconsistent_units",
            f"第 {sequence} 段时间单位 {seg_time_unit!r} 与运行级单位 {time_unit!r} 不一致",
        )
    length = _finite_positive("length", raw.get("length")) * lf
    velocity = _finite_positive("pore_velocity", raw.get("pore_velocity")) * lf / tf
    dispersion = _finite_positive("dispersion", raw.get("dispersion")) * lf * lf / tf
    decay_raw = raw.get("decay_rate", 0.0) or 0.0
    try:
        decay = float(decay_raw) / tf
    except (TypeError, ValueError) as exc:
        raise TransportError("invalid_transport_config", "decay_rate 不是数值") from exc
    if not math.isfinite(decay) or decay < 0:
        raise TransportError("invalid_transport_config", "decay_rate 不能为负")
    code = raw.get("code")
    if code is not None:
        if code in used_codes:
            raise TransportError("duplicate_segment_code", f"区段编码 {code!r} 重复")
        used_codes.add(code)
    parameter_version = str(raw.get("parameter_version") or "unspecified")
    return Segment(sequence, code, length, velocity, dispersion, decay, parameter_version)


def normalize_config(payload: dict[str, Any]) -> TransportConfig:
    """校验并把请求负载换算为内部规范单位配置。"""
    length_unit = payload.get("length_unit") or "m"
    time_unit = payload.get("time_unit") or "day"
    if length_unit not in LENGTH_TO_M:
        raise TransportError("unsupported_unit", f"不支持的长度单位 {length_unit!r}，仅支持 m / km")
    if time_unit not in TIME_TO_DAY:
        raise TransportError("unsupported_unit", f"不支持的时间单位 {time_unit!r}，仅支持 day / hour / second")
    lf = LENGTH_TO_M[length_unit]
    tf = TIME_TO_DAY[time_unit]

    source = payload.get("source_mass")
    if source is None:
        source = payload.get("source_concentration")
    source_mass = _finite_positive("source_mass", source)

    raw_segments = payload.get("segments")
    used_codes: set[str] = set()
    if raw_segments:
        if not isinstance(raw_segments, list) or not raw_segments:
            raise TransportError("segments_required", "segments 必须为非空有序列表")
        segments = tuple(
            _build_segment(raw, i, lf, tf, length_unit, time_unit, used_codes)
            for i, raw in enumerate(raw_segments, start=1)
        )
        sequences = [segment.sequence for segment in segments]
        if sequences != list(range(1, len(segments) + 1)):
            raise TransportError(
                "segments_out_of_order",
                "区段必须按从源区到监测井的顺序给出，sequence 从 1 开始连续递增、不得缺号或重号",
                {"received": sequences},
            )
        duration_days = _finite_positive("duration", payload.get("duration")) * tf
        step_days = _finite_positive("step", payload.get("step")) * tf
    else:
        # 兼容单一参数（旧版接口）：整条路径视为一个区段
        missing = [
            name for name in ("distance_m", "velocity_m_day", "dispersion_m2_day", "duration_days")
            if payload.get(name) is None
        ]
        if missing:
            raise TransportError(
                "segments_required",
                f"未提供 segments 时必须给出单一区段参数，缺少：{', '.join(missing)}",
            )
        decay = float(payload.get("decay_per_day") or 0.0)
        if not math.isfinite(decay) or decay < 0:
            raise TransportError("invalid_transport_config", "decay_per_day 不能为负")
        segments = (Segment(
            sequence=1,
            code=None,
            length_m=_finite_positive("distance_m", payload["distance_m"]),
            velocity_m_day=_finite_positive("velocity_m_day", payload["velocity_m_day"]),
            dispersion_m2_day=_finite_positive("dispersion_m2_day", payload["dispersion_m2_day"]),
            decay_per_day=decay,
            parameter_version="legacy-single-segment",
        ),)
        duration_days = _finite_positive("duration_days", payload["duration_days"])
        step_value = payload.get("step_days")
        step_days = _finite_positive("step_days", 1.0 if step_value is None else step_value)

    if step_days > duration_days:
        raise TransportError("invalid_transport_config", "时间步长不得大于模拟时长")

    detection_limit = float(payload.get("detection_limit") or 0.0)
    relative_threshold = float(payload.get("relative_threshold") if payload.get("relative_threshold") is not None else 1e-3)
    tolerance = float(payload.get("mass_error_tolerance") if payload.get("mass_error_tolerance") is not None else 0.01)
    if not math.isfinite(detection_limit) or detection_limit < 0:
        raise TransportError("invalid_transport_config", "detection_limit 不能为负")
    if not 0 < relative_threshold < 1:
        raise TransportError("invalid_transport_config", "relative_threshold 必须位于 (0, 1)")
    if not 1e-6 <= tolerance <= 0.25:
        raise TransportError("invalid_transport_config", "mass_error_tolerance 必须位于 [1e-6, 0.25]")

    return TransportConfig(
        source_mass=source_mass,
        duration_days=duration_days,
        step_days=step_days,
        segments=segments,
        detection_limit=detection_limit,
        relative_threshold=relative_threshold,
        mass_error_tolerance=tolerance,
        model_version=str(payload.get("model_version") or SOLVER_VERSION),
        length_unit=length_unit,
        time_unit=time_unit,
        time_factor=tf,
    )


def _kernel_value(t: float, segment: Segment) -> float:
    """带一阶衰减的通量型逆高斯核 g_i(t)，t>0。"""
    if t <= 0:
        return 0.0
    ldv = segment.length_m - segment.velocity_m_day * t
    exponent = -(ldv * ldv) / (4.0 * segment.dispersion_m2_day * t) - segment.decay_per_day * t
    if exponent < -745.0:
        return 0.0
    return (
        segment.length_m
        / math.sqrt(4.0 * math.pi * segment.dispersion_m2_day * t ** 3)
        * math.exp(exponent)
    )


def _survival_factor(segment: Segment) -> float:
    ratio = 4.0 * segment.dispersion_m2_day * segment.decay_per_day / segment.velocity_m_day ** 2
    beta = segment.velocity_m_day * segment.length_m / (2.0 * segment.dispersion_m2_day)
    return math.exp(beta * (1.0 - math.sqrt(1.0 + ratio)))


def _segment_lags(segment: Segment, h: float) -> list[tuple[int, float]]:
    """在精细时间网格上离散化 g_i：返回 (滞后步数 k, 权重 g_i(kh)·h)，端点用梯形半权。"""
    mean = segment.length_m / segment.velocity_m_day
    variance = 2.0 * segment.dispersion_m2_day * segment.length_m / segment.velocity_m_day ** 3
    spread = math.sqrt(max(variance, 0.0))
    peak = _kernel_value(mean, segment)
    threshold = _KERNEL_CUTOFF * max(peak, 1e-300)

    low = mean
    walk = max(h, spread * 0.5)
    for _ in range(200):
        nxt = max(h, low - walk)
        low = nxt
        if _kernel_value(nxt, segment) <= threshold or nxt <= h:
            break
    high = mean
    for _ in range(200):
        nxt = high + walk
        high = nxt
        if _kernel_value(nxt, segment) <= threshold:
            break

    k0 = max(1, int(math.floor(low / h)))
    k1 = max(k0 + 1, int(math.ceil(high / h)))
    lags: list[tuple[int, float]] = []
    for k in range(k0, k1 + 1):
        weight = _kernel_value(k * h, segment) * h
        if k in (k0, k1):
            weight *= 0.5
        if weight > 0.0:
            lags.append((k, weight))
    return lags


def _choose_grid(config: TransportConfig) -> tuple[float, int, int]:
    """返回 (精细步长 h, 窗口步数 J, 含尾部核算的扩展步数 Jext)。"""
    means = [s.length_m / s.velocity_m_day for s in config.segments]
    variances = [
        2.0 * s.dispersion_m2_day * s.length_m / s.velocity_m_day ** 3
        for s in config.segments
    ]
    target = min(
        min(math.sqrt(v) / 8.0 if v > 0 else math.inf for v in variances),
        min(means) / 12.0,
    )
    steps_per_report = max(1, int(math.ceil(config.step_days / target)))
    h = config.step_days / steps_per_report
    j_window = int(math.ceil(config.duration_days / h - 1e-12))
    if j_window > _MAX_REPORT_STEPS:  # 极端高 Peclet 数下放宽网格，交由质量误差门限拒绝
        steps_per_report = max(1, int(config.step_days / (config.duration_days / _MAX_REPORT_STEPS)))
        h = config.step_days / steps_per_report
        j_window = int(round(config.duration_days / h))
    tail_cover = sum(means) + 8.0 * math.sqrt(sum(variances))
    j_extended = max(j_window, int(math.ceil(tail_cover / h)))
    if j_extended > _MAX_INTERNAL_STEPS:
        j_extended = j_window
    return h, j_window, j_extended


def _trapezoid(curve: list[float], h: float, upto: int) -> float:
    if upto <= 0:
        return 0.0
    total = 0.5 * curve[0] + 0.5 * curve[upto]
    total += sum(curve[1:upto])
    return total * h


def _cumulative(curve: list[float], h: float) -> list[float]:
    cumulative = [0.0] * len(curve)
    for j in range(1, len(curve)):
        cumulative[j] = cumulative[j - 1] + 0.5 * (curve[j - 1] + curve[j]) * h
    return cumulative


def parameter_fingerprint(config: TransportConfig) -> str:
    """参数版本指纹：唯一标识本次运行采用的物理参数版本组合。"""
    canonical = {
        "segments": [
            {
                "length_m": s.length_m,
                "velocity_m_day": s.velocity_m_day,
                "dispersion_m2_day": s.dispersion_m2_day,
                "decay_per_day": s.decay_per_day,
                "parameter_version": s.parameter_version,
            }
            for s in config.segments
        ],
        "source_mass": config.source_mass,
    }
    body = json.dumps(canonical, sort_keys=True, ensure_ascii=False).encode()
    return hashlib.sha256(body).hexdigest()[:16]


def simulate(config: TransportConfig) -> dict[str, Any]:
    """执行分段卷积模拟，返回突破曲线与指标；质量误差超限时抛出 TransportError。"""
    h, j_window, j_extended = _choose_grid(config)
    lag_tables = [_segment_lags(segment, h) for segment in config.segments]

    # 源项：t=0 处的单位面积脉冲 M0，离散为 M0/h
    curve = [0.0] * (j_extended + 1)
    curve[0] = config.source_mass / h
    for lags in lag_tables:
        outgoing = [0.0] * (j_extended + 1)
        for p, value in enumerate(curve):
            if value == 0.0:
                continue
            for k, weight in lags:
                j = p + k
                if j > j_extended:
                    break
                outgoing[j] += value * weight
        curve = outgoing

    survivals = [_survival_factor(s) for s in config.segments]
    expected_mass = config.source_mass
    interface_rows: list[dict[str, Any]] = []
    travel_time = 0.0
    interface_rows.append({
        "interface": 0,
        "travel_time_days": 0.0,
        "cumulative_mass": config.source_mass,
    })
    segment_rows = []
    for segment, survival in zip(config.segments, survivals):
        expected_mass *= survival
        travel_time += segment.length_m / segment.velocity_m_day
        interface_rows.append({
            "interface": segment.sequence,
            "travel_time_days": travel_time,
            "cumulative_mass": expected_mass,
        })
        segment_rows.append({
            "sequence": segment.sequence,
            "code": segment.code,
            "length_m": segment.length_m,
            "velocity_m_day": segment.velocity_m_day,
            "dispersion_m2_day": segment.dispersion_m2_day,
            "decay_per_day": segment.decay_per_day,
            "parameter_version": segment.parameter_version,
            "travel_time_days": segment.length_m / segment.velocity_m_day,
            "survival_factor": survival,
        })

    denominator = max(expected_mass, 1e-12 * config.source_mass)
    if j_extended > j_window:
        extended_mass = _trapezoid(curve, h, j_extended)
        solver_error = abs(extended_mass - expected_mass) / denominator
        if solver_error > config.mass_error_tolerance:
            raise TransportError(
                "mass_error_exceeded",
                f"分段卷积质量误差 {solver_error:.4%} 超过门限 {config.mass_error_tolerance:.4%}，结果被拒绝",
                {"mass_error": solver_error, "tolerance": config.mass_error_tolerance},
            )

    window_mass = _trapezoid(curve, h, j_window)
    mass_error = abs(window_mass - expected_mass) / denominator
    if mass_error > config.mass_error_tolerance:
        means = sum(s.length_m / s.velocity_m_day for s in config.segments)
        variance = sum(
            2.0 * s.dispersion_m2_day * s.length_m / s.velocity_m_day ** 3
            for s in config.segments
        )
        suggested_days = means + 8.0 * math.sqrt(variance)
        raise TransportError(
            "mass_error_exceeded",
            f"模拟时段内累计质量 {window_mass:.6g} 与应到达质量 {expected_mass:.6g} "
            f"相差 {mass_error:.4%}，超过门限 {config.mass_error_tolerance:.4%}；"
            "污染羽尾部尚未完全通过目标井，请延长模拟时长后重算",
            {
                "mass_error": mass_error,
                "tolerance": config.mass_error_tolerance,
                "window_mass": window_mass,
                "expected_mass": expected_mass,
                "suggested_duration": suggested_days / config.time_factor,
                "time_unit": config.time_unit,
            },
        )

    cumulative = _cumulative(curve, h)

    # 峰值（抛物线插值修正）
    peak_j = max(range(1, j_window + 1), key=lambda j: curve[j])
    y0, y1, y2 = curve[peak_j - 1], curve[peak_j], curve[min(peak_j + 1, j_window)]
    bend = y0 - 2.0 * y1 + y2
    if bend < 0.0:
        offset = 0.5 * (y0 - y2) / bend
        offset = max(-1.0, min(1.0, offset))
    else:
        offset = 0.0
    peak_time = (peak_j + offset) * h
    peak_concentration = y1 - 0.25 * (y0 - y2) * offset
    peak = {
        "time_days": peak_time,
        "concentration": peak_concentration,
        "cumulative_mass": cumulative[peak_j],
    }

    # 首达时间：阈值 = max(绝对检出限, 峰值的相对比例)，线性插值穿越时刻
    threshold = max(config.detection_limit, config.relative_threshold * peak_concentration)
    first_arrival: float | None = None
    for j in range(1, j_window + 1):
        if curve[j - 1] < threshold <= curve[j] and curve[j] > curve[j - 1]:
            first_arrival = (j - 1 + (threshold - curve[j - 1]) / (curve[j] - curve[j - 1])) * h
            break

    steps_per_report = max(1, int(round(config.step_days / h)))
    points = []
    j = steps_per_report
    while j <= j_window:
        points.append({
            "time_days": round(j * h, 10),
            "concentration": curve[j],
            "cumulative_mass": cumulative[j],
        })
        j += steps_per_report

    return {
        "model_version": config.model_version,
        "solver_version": SOLVER_VERSION,
        "parameter_fingerprint": parameter_fingerprint(config),
        "unit_system": {"length": "m", "time": "day",
                        "requested_length": config.length_unit, "requested_time": config.time_unit},
        "segment_count": len(config.segments),
        "segments": segment_rows,
        "total_length_m": sum(s.length_m for s in config.segments),
        "arrival_time_days": sum(s.length_m / s.velocity_m_day for s in config.segments),
        "source_mass": config.source_mass,
        "points": points,
        "peak": peak,
        "first_arrival_days": first_arrival,
        "detection_threshold": threshold,
        "cumulative_mass": window_mass,
        "expected_mass": expected_mass,
        "mass_error": mass_error,
        "mass_error_tolerance": config.mass_error_tolerance,
        "interface_flux": interface_rows,
    }
