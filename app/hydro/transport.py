"""分段一维平流-弥散-一阶衰减(ADR)迁移求解器。

物理模型
--------
污染羽从源区到监测井依次穿过若干渗透系数与衰减条件不同的含水层区段。对每个区段,
一维 ADR 方程 ``dc/dt = -v dc/dx + D d2c/dx2 - lambda c`` 在瞬时面源注入、上游
无限远边界下,出口处的单位脉冲响应(首达时间密度与衰减存活因子的乘积)为

    g_i(t) = L_i / sqrt(4*pi*D_i*t^3)
             * exp(-(L_i - v_i*t)^2 / (4*D_i*t))
             * exp(-lambda_i*t)

各区段线性串联时,系统总响应是各段脉冲响应的卷积。把区段 i 的出流质量通量作为
区段 i+1 的入流通量做卷积,即满足**区段接口上的质量通量连续**:无衰减时
``integral(g_out dt) == 1``,衰减只发生在区段内部,不会在接口处凭空产生或截留质量。

数值上采用等时间步 FFT 卷积(零填充到 2 的幂以避免循环卷绕)逐级传递通量序列,
并用各区段解析存活概率的乘积给出期望到达质量,对数值累计质量做守恒校验。

单位制:长度 m、时间 day、质量 kg,因此孔隙流速 m/day、弥散系数 m2/day、一阶
衰减 1/day、质量通量 kg/day。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

# FFT 内部最大网格,超出时改用直接(时域)卷积或报错前的保护上限
_MAX_FFT_SIZE = 16384


class TransportValidationError(ValueError):
    """输入配置无法形成物理上合法的分段迁移计算。"""


class MassBalanceError(RuntimeError):
    """数值累计质量与解析期望到达质量的偏差超过容差,结果被拒绝。"""


@dataclass(frozen=True)
class Segment:
    length_m: float
    velocity_m_day: float
    dispersion_m2_day: float
    decay_per_day: float


def segment_survival(segment: Segment) -> float:
    """单个区段一阶衰减下的解析质量存活比例。"""
    l, v, d, lam = (
        segment.length_m,
        segment.velocity_m_day,
        segment.dispersion_m2_day,
        segment.decay_per_day,
    )
    if lam <= 0.0:
        return 1.0
    return math.exp(-2.0 * lam * l / (v + math.sqrt(v * v + 4.0 * d * lam)))


def segment_mean_travel_time(segment: Segment) -> float:
    """单段平流参考时间 L/v,用于和无弥散情形对比。"""
    return segment.length_m / segment.velocity_m_day


def impulse_response(times: list[float], segment: Segment) -> list[float]:
    """区段在给定时间网格(天)上的单位脉冲响应 g_i(t)。"""
    l, v, d, lam = (
        segment.length_m,
        segment.velocity_m_day,
        segment.dispersion_m2_day,
        segment.decay_per_day,
    )
    prefactor = l / math.sqrt(4.0 * math.pi * d)
    out: list[float] = []
    for t in times:
        if t <= 0.0:
            out.append(0.0)
            continue
        value = prefactor * t ** -1.5 * math.exp(-((l - v * t) ** 2) / (4.0 * d * t) - lam * t)
        out.append(value if value >= 0.0 else 0.0)
    return out


def _convolve(a: list[float], b: list[float], limit: int) -> list[float]:
    """离散卷积 a*b 的前 limit 项;小核直接时域计算,大核走 FFT 避免 O(n^2)。"""
    if not a or not b:
        return [0.0] * limit
    if len(b) <= 64:
        # 直接卷积:out[n] = sum_j b[j]*a[n-j]
        out = [0.0] * limit
        for j, bj in enumerate(b):
            if bj == 0.0:
                continue
            upper = min(limit, len(a) + j)
            for n in range(j, upper):
                if n - j < len(a):
                    out[n] += bj * a[n - j]
        return out

    size = 1
    needed = len(a) + len(b) - 1
    while size < needed:
        size <<= 1
    if size > _MAX_FFT_SIZE:
        # 回退到分块直接卷积,保持接口稳定
        out = [0.0] * min(limit, needed)
        for i, ai in enumerate(a):
            if ai == 0.0:
                continue
            for j, bj in enumerate(b):
                n = i + j
                if n >= len(out):
                    break
                out[n] += ai * bj
        if len(out) < limit:
            out.extend([0.0] * (limit - len(out)))
        return out

    fa = _fft(a, size)
    fb = _fft(b, size)
    spectrum = [zr * z for zr, z in zip(fa, fb)]
    full = _ifft(spectrum, size)
    out = full[:limit]
    return out


# --------------------------------------------------------------------------- #
# 极简 radix-2 FFT(只依赖标准库);旋转因子用查找表一次性计算,避免迭代累乘漂移,
# 数值精度对本问题的质量守恒校验完全足够
# --------------------------------------------------------------------------- #
def _fft_transform(data: list[complex], inverse: bool) -> list[complex]:
    size = len(data)
    # 位反转置换
    j = 0
    for i in range(1, size):
        bit = size >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            data[i], data[j] = data[j], data[i]
    sign = 1.0 if inverse else -1.0
    roots = [
        complex(math.cos(sign * math.pi * k / (size >> 1)), math.sin(sign * math.pi * k / (size >> 1)))
        for k in range(size >> 1)
    ]
    length = 2
    while length <= size:
        half = length >> 1
        stride = size // length
        for start in range(0, size, length):
            for k in range(half):
                w = roots[k * stride]
                u = data[start + k]
                v = data[start + k + half] * w
                data[start + k] = u + v
                data[start + k + half] = u - v
        length <<= 1
    if inverse:
        data = [z / size for z in data]
    return data


def _fft(values: list[float], size: int) -> list[complex]:
    data = [complex(values[i]) if i < len(values) else 0j for i in range(size)]
    return _fft_transform(data, inverse=False)


def _ifft(spectrum: list[complex], size: int) -> list[float]:
    return [z.real for z in _fft_transform(list(spectrum), inverse=True)]


def validate_segments(raw_segments: list[dict[str, Any]], distance_m: float) -> list[Segment]:
    """校验区段有序、长度为正且总长度与源-井距离一致,返回有序 Segment 列表。

    顺序错误无法直接从单条记录判别,但本函数强制:至少一段、各段长度严格为正、
    总长度闭合到源-井距离;调用方必须按源区 -> 监测井的方向组装列表。
    """
    if not raw_segments:
        raise TransportValidationError("segments_required:至少需要一个含水层区段")
    segments: list[Segment] = []
    total = 0.0
    for index, raw in enumerate(raw_segments):
        try:
            segment = Segment(
                length_m=float(raw["length_m"]),
                velocity_m_day=float(raw["velocity_m_day"]),
                dispersion_m2_day=float(raw["dispersion_m2_day"]),
                decay_per_day=float(raw.get("decay_per_day", 0.0)),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise TransportValidationError(
                f"segment_{index}_invalid:区段 {index + 1} 参数缺失或不是数值"
            ) from exc
        if not math.isfinite(segment.length_m) or segment.length_m <= 0:
            raise TransportValidationError(
                f"segment_{index}_length_error:区段 {index + 1} 长度必须为正"
            )
        if not math.isfinite(segment.velocity_m_day) or segment.velocity_m_day <= 0:
            raise TransportValidationError(
                f"segment_{index}_velocity_error:区段 {index + 1} 孔隙流速必须为正"
            )
        if not math.isfinite(segment.dispersion_m2_day) or segment.dispersion_m2_day <= 0:
            raise TransportValidationError(
                f"segment_{index}_dispersion_error:区段 {index + 1} 弥散系数必须为正"
            )
        if not math.isfinite(segment.decay_per_day) or segment.decay_per_day < 0:
            raise TransportValidationError(
                f"segment_{index}_decay_error:区段 {index + 1} 一阶衰减不能为负"
            )
        segments.append(segment)
        total += segment.length_m
    if abs(total - distance_m) > max(1e-9, 0.005 * distance_m):
        raise TransportValidationError(
            f"segment_chain_unclosed:区段长度合计 {total:g} m 与源区至监测井距离 "
            f"{distance_m:g} m 不一致,区段顺序或长度有误"
        )
    # 若调用方提供了里程坐标,复核区段接口衔接与单调方向
    with_coords = [
        raw for raw in raw_segments
        if raw.get("start_m") is not None or raw.get("end_m") is not None
    ]
    if with_coords:
        if len(with_coords) != len(raw_segments):
            raise TransportValidationError("segment_order_error:里程坐标必须全部区段同时提供")
        cursor = 0.0
        for index, raw in enumerate(raw_segments):
            if raw.get("start_m") is None or raw.get("end_m") is None:
                raise TransportValidationError(
                    f"segment_order_error:第 {index + 1} 段里程坐标不完整"
                )
            start_m, end_m = float(raw["start_m"]), float(raw["end_m"])
            if abs(start_m - cursor) > 1e-6 or end_m <= start_m:
                raise TransportValidationError(
                    f"segment_order_error:第 {index + 1} 段未在里程 {cursor:g} m 处衔接或方向倒退"
                )
            if abs(end_m - start_m - segments[index].length_m) > 1e-6:
                raise TransportValidationError(
                    f"segment_order_error:第 {index + 1} 段里程跨度与长度不一致"
                )
            cursor = end_m
        if abs(cursor - distance_m) > max(1e-6, 0.005 * distance_m):
            raise TransportValidationError(
                f"segment_order_error:末段终点里程 {cursor:g} m 未到达目标井"
            )
    return segments


def solve_segment_transport(payload: dict[str, Any]) -> dict[str, Any]:
    """执行分段迁移计算并返回突破曲线与汇总指标。

    输入 payload 字段(经 SegmentTransportRequest 校验):
      source_mass_kg, distance_m, segments[], duration_days, step_days,
      first_arrival_quantile, mass_tolerance, model_version, parameter_version
    """
    distance_m = float(payload["distance_m"])
    step_days = float(payload["step_days"])
    duration_days = float(payload["duration_days"])
    source_mass_kg = float(payload["source_mass_kg"])
    quantile = float(payload["first_arrival_quantile"])
    tolerance = float(payload["mass_tolerance"])

    if step_days <= 0 or duration_days <= 0 or step_days > duration_days:
        raise TransportValidationError(
            "time_grid_invalid:step_days 必须为正且不大于 duration_days"
        )

    segments = validate_segments(payload["segments"], distance_m)

    n_steps = int(math.floor(duration_days / step_days + 1e-9))
    if n_steps < 2:
        raise TransportValidationError("time_grid_invalid:时间窗口内至少要有两个采样点")
    # 网格从 t=0 起:脉冲在源区 t=0 时刻注入,通量网格 f[n] 对应 t=n*h
    times = [step_days * i for i in range(n_steps + 1)]

    # 逐级卷积:第一段的入流是 t=0 的瞬时单位质量脉冲(通量序列 1/h),
    # 每经过一段就与该段脉冲响应卷积,接口处通量序列直接传递,保证质量通量连续。
    # 先用单位质量传递,末端再乘源质量,避免大质量初值放大 FFT 浮点噪声。
    segment_survivals = [segment_survival(segment) for segment in segments]
    flux = [1.0 / step_days] + [0.0] * n_steps
    for segment, survival in zip(segments, segment_survivals):
        # 离散核:梯形权(t=0 处半权,且 g(0)=0),再按该段解析存活率归一化,
        # 消除离散求积偏差在多段串联中的累积;接口传递的总质量即精确等于 S_i
        raw = impulse_response(times, segment)
        kernel = [step_days * value for value in raw]
        kernel[0] *= 0.5
        kernel_sum = sum(kernel)
        if kernel_sum <= 0.0:
            raise TransportValidationError(
                "time_grid_too_coarse:时间步长相对区段运移时间过大,脉冲响应采样不到质量,"
                "请减小 step_days"
            )
        scale = survival / kernel_sum
        kernel = [value * scale for value in kernel]
        flux = [value if value > 0.0 else 0.0 for value in _convolve(flux, kernel, n_steps + 1)]

    # 期望到达质量:各段解析存活概率之积
    survival = 1.0
    for segment in segments:
        survival *= segment_survival(segment)
    expected_mass = source_mass_kg * survival
    # 梯形积分累计到达质量(t=0 处通量为 0)
    cumulative = 0.0
    points: list[dict[str, Any]] = []
    peak_flux = 0.0
    peak_time: float | None = None
    first_arrival_time: float | None = None
    first_arrival_target = expected_mass * quantile
    for n in range(1, n_steps + 1):
        t = times[n]
        value = flux[n] * source_mass_kg
        cumulative += 0.5 * (flux[n - 1] + flux[n]) * source_mass_kg * step_days
        if peak_time is None or value > peak_flux:
            peak_flux = value
            peak_time = t
        if first_arrival_time is None and cumulative >= first_arrival_target:
            first_arrival_time = t
        points.append({"time_days": round(t, 8), "mass_flux_kg_day": value})

    mass_error = abs(cumulative - expected_mass) / expected_mass if expected_mass > 0 else 0.0
    if expected_mass > 0 and mass_error > tolerance:
        raise MassBalanceError(
            f"mass_balance_exceeded:窗口内累计质量 {cumulative:.6g} kg 与解析期望 "
            f"{expected_mass:.6g} kg 的相对偏差 {mass_error:.4f} 超过容差 {tolerance:.4f};"
            "请增大 duration_days 以覆盖完整突破过程或减小 step_days"
        )

    # 首达时间:累计到达质量首次达到 expected_mass*quantile 的网格时刻;
    # 衰减情形下即“存活质量分位数”的到达时间
    advective_arrival = sum(segment_mean_travel_time(segment) for segment in segments)

    return {
        "points": points,
        "peak": {
            "time_days": round(peak_time, 8) if peak_time is not None else None,
            "mass_flux_kg_day": peak_flux,
        },
        "first_arrival_time_days": round(first_arrival_time, 8)
        if first_arrival_time is not None
        else None,
        "advective_travel_time_days": advective_arrival,
        "cumulative_mass_kg": cumulative,
        "expected_mass_kg": expected_mass,
        "mass_error": mass_error,
        "survival_fraction": survival,
        "n_segments": len(segments),
        "total_distance_m": sum(s.length_m for s in segments),
        "model_version": payload["model_version"],
        "parameter_version": payload["parameter_version"],
        "units": payload.get("units", {}),
        "segment_summary": [
            {
                "index": i + 1,
                "length_m": s.length_m,
                "velocity_m_day": s.velocity_m_day,
                "dispersion_m2_day": s.dispersion_m2_day,
                "decay_per_day": s.decay_per_day,
                "travel_time_days": segment_mean_travel_time(s),
                "survival": segment_survival(s),
            }
            for i, s in enumerate(segments)
        ],
    }
