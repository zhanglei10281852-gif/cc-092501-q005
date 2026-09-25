from __future__ import annotations

from pydantic import BaseModel, Field, field_validator, model_validator

# 分段迁移求解器允许的最大时间步数,超出时要求调用方增大步长,避免纯 Python FFT 过慢
MAX_TRANSPORT_GRID_STEPS = 2000

# 分段迁移计算的规范单位制;请求中声明的单位必须与此完全一致,否则视为单位不一致
CANONICAL_TRANSPORT_UNITS = {
    "length": "m",
    "time": "day",
    "mass": "kg",
    "velocity": "m/day",
    "dispersion": "m2/day",
    "decay": "1/day",
}


class WellCreate(BaseModel):
    code: str = Field(..., min_length=2, max_length=50)
    name: str = Field(..., min_length=1, max_length=120)
    latitude: float = Field(..., ge=-90, le=90)
    longitude: float = Field(..., ge=-180, le=180)
    aquifer: str = Field(..., min_length=1, max_length=120)
    screen_depth_m: float = Field(..., gt=0, le=5000)

    @field_validator("code")
    @classmethod
    def normalize_code(cls, value: str) -> str:
        return value.strip().upper()


class EndmemberCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=100)
    isotope_d18o: float = Field(..., ge=-100, le=100)
    isotope_d2h: float = Field(..., ge=-800, le=800)
    solute_mg_l: float = Field(..., ge=0, le=100000)
    uncertainty: float = Field(default=0.1, gt=0, le=100)
    version: str = Field(default="v1", min_length=1, max_length=40)


class SampleCreate(BaseModel):
    sample_code: str = Field(..., min_length=3, max_length=64)
    sampled_at: str = Field(..., min_length=20, max_length=40)
    isotope_d18o: float | None = Field(default=None, ge=-100, le=100)
    isotope_d2h: float | None = Field(default=None, ge=-800, le=800)
    solute_mg_l: float | None = Field(default=None, ge=0, le=100000)
    detection_limit: float = Field(default=0, ge=0, le=100000)
    measurement_error: float = Field(default=0.05, ge=0, le=100)


class InversionRequest(BaseModel):
    endmember_ids: list[int] = Field(..., min_length=2, max_length=8)
    method: str = Field(default="weighted-least-squares", pattern="^(weighted-least-squares|projected-gradient)$")
    max_iterations: int = Field(default=500, ge=10, le=10000)
    tolerance: float = Field(default=1e-8, gt=0, le=0.1)
    model_version: str = Field(default="mix-1", min_length=1, max_length=40)


class TransportRequest(BaseModel):
    source_concentration: float = Field(..., ge=0, le=1000000)
    distance_m: float = Field(..., gt=0, le=1000000)
    velocity_m_day: float = Field(..., gt=0, le=10000)
    dispersion_m2_day: float = Field(..., gt=0, le=100000)
    decay_per_day: float = Field(default=0, ge=0, le=100)
    duration_days: float = Field(..., gt=0, le=100000)
    step_days: float = Field(default=1, gt=0, le=1000)
    model_version: str = Field(default="ade-1", min_length=1, max_length=40)


class TransportSegment(BaseModel):
    """从源区向监测井方向的一个含水层区段,字段顺序即迁移方向。

    可选的 start_m/end_m 是相对源区的里程坐标,用于校验区段接口衔接与方向;
    一旦任一区段给出坐标,所有区段都必须给出且必须从 0 单调衔接到 distance_m。
    """

    length_m: float = Field(..., gt=0, le=1000000)
    velocity_m_day: float = Field(..., gt=0, le=10000)
    dispersion_m2_day: float = Field(..., gt=0, le=100000)
    decay_per_day: float = Field(default=0, ge=0, le=100)
    start_m: float | None = Field(default=None, ge=0, le=1000000)
    end_m: float | None = Field(default=None, ge=0, le=1000000)


class SegmentTransportRequest(BaseModel):
    """有序区段的一维平流-弥散-一阶衰减迁移计算请求。

    区段必须按污染源 -> 监测井的顺序给出;总路径长度必须与 distance_m 一致(允许
    0.5% 的坐标误差);所有量纲固定为 SI 风格的 m/day/kg 单位制,由 units 显式确认。
    """

    source_mass_kg: float = Field(..., gt=0, le=1000000)
    distance_m: float = Field(..., gt=0, le=1000000)
    segments: list[TransportSegment] = Field(..., min_length=1, max_length=20)
    duration_days: float = Field(..., gt=0, le=100000)
    step_days: float = Field(default=1, gt=0, le=1000)
    first_arrival_quantile: float = Field(default=0.01, gt=0, lt=0.5)
    mass_tolerance: float = Field(default=0.02, gt=0, le=0.2)
    model_version: str = Field(default="ade-segment-1", min_length=1, max_length=40)
    parameter_version: str = Field(default="param-1", min_length=1, max_length=40)
    units: dict[str, str] = Field(default_factory=dict, validate_default=True)

    @field_validator("units")
    @classmethod
    def validate_units(cls, value: dict[str, str]) -> dict[str, str]:
        if value != CANONICAL_TRANSPORT_UNITS:
            raise ValueError(
                "单位不一致:长度 m、时间 day、质量 kg、速度 m/day、弥散系数 m2/day、衰减 1/day,"
                "需在 units 中显式声明全部规范单位"
            )
        return value

    @model_validator(mode="after")
    def validate_segments(self) -> "SegmentTransportRequest":
        total = sum(segment.length_m for segment in self.segments)
        if abs(total - self.distance_m) > max(1e-9, 0.005 * self.distance_m):
            raise ValueError(
                f"区段顺序/长度错误:各区段长度之和 {total:g} m 与源区至目标井距离 "
                f"{self.distance_m:g} m 不一致(允许 0.5% 误差)"
            )
        with_coords = [s for s in self.segments if s.start_m is not None or s.end_m is not None]
        if with_coords:
            if len(with_coords) != len(self.segments):
                raise ValueError("区段顺序错误:里程坐标必须全部区段同时提供或同时省略")
            cursor = 0.0
            for index, segment in enumerate(self.segments):
                assert segment.start_m is not None and segment.end_m is not None
                if abs(segment.start_m - cursor) > 1e-6:
                    raise ValueError(
                        f"区段顺序错误:第 {index + 1} 段起点里程 {segment.start_m:g} m "
                        f"未与上一段终点 {cursor:g} m 衔接"
                    )
                if segment.end_m <= segment.start_m:
                    raise ValueError(
                        f"区段顺序错误:第 {index + 1} 段终点里程不大于起点,方向与源区至监测井不符"
                    )
                if abs(segment.end_m - segment.start_m - segment.length_m) > max(
                    1e-6, 0.005 * segment.length_m
                ):
                    raise ValueError(
                        f"区段顺序错误:第 {index + 1} 段里程跨度与长度 {segment.length_m:g} m 不一致"
                    )
                cursor = segment.end_m
            if abs(cursor - self.distance_m) > max(1e-6, 0.005 * self.distance_m):
                raise ValueError(
                    f"区段顺序错误:末段终点里程 {cursor:g} m 未到达目标井距离 {self.distance_m:g} m"
                )
        if self.duration_days / self.step_days > MAX_TRANSPORT_GRID_STEPS:
            raise ValueError(
                f"时间网格过大:{self.duration_days / self.step_days:.0f} 步超过上限 "
                f"{MAX_TRANSPORT_GRID_STEPS},请增大 step_days 或缩短 duration_days"
            )
        return self

