import math
from typing import Iterable


class ContinuousScalarProjector:
    """Map continuous clinical values into a bounded Fourier feature vector."""

    def __init__(self, harmonics: int = 8) -> None:
        if harmonics < 1:
            raise ValueError("harmonics must be at least 1")
        self.harmonics = harmonics

    def project(self, value: float, scale: float = 1.0) -> list[float]:
        if not math.isfinite(value) or not math.isfinite(scale) or scale <= 0:
            raise ValueError("value and scale must be finite; scale must be positive")
        normalized = value / scale
        result: list[float] = []
        for k in range(self.harmonics):
            angle = (2.0 ** k) * math.pi * normalized
            result.extend((math.sin(angle), math.cos(angle)))
        return result

    def distances(self, value: float, others: Iterable[float], scale: float = 1.0) -> list[float]:
        vector = self.project(value, scale)
        scores: list[float] = []
        for other in others:
            candidate = self.project(other, scale)
            scores.append(math.sqrt(sum((a - b) ** 2 for a, b in zip(vector, candidate))))
        return scores
