from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class VerificationResult:
    passed: bool
    violations: tuple[str, ...]


class BiochemicalInvariantGate:
    def verify(self, claims: Iterable[str], allowed_claims: Iterable[str]) -> VerificationResult:
        allowed = {claim.strip().lower() for claim in allowed_claims}
        violations = tuple(
            claim for claim in claims
            if claim.strip().lower() not in allowed
        )
        return VerificationResult(passed=not violations, violations=violations)


class ClearanceParityChecker:
    def verify(self, drug: str, prescribed_dose: float, egfr: float, max_safe_dose: float) -> VerificationResult:
        violations: list[str] = []
        if prescribed_dose < 0:
            violations.append(f"{drug}: dose cannot be negative")
        if egfr < 30 and prescribed_dose > max_safe_dose:
            violations.append(f"{drug}: dose exceeds renal safety limit")
        return VerificationResult(passed=not violations, violations=tuple(violations))
