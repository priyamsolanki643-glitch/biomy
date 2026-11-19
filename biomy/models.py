from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class AnalyteNode:
    name: str
    value: float
    unit: str
    reference_low: Optional[float] = None
    reference_high: Optional[float] = None
    source_text: str = ""


@dataclass(frozen=True)
class ExposureNode:
    name: str
    kind: str
    strength: Optional[float] = None
    unit: str = ""
    route: str = ""
    source_text: str = ""
