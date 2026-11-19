from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class ClinicalReport:
    title: str
    findings: tuple[str, ...]
    evidence: tuple[str, ...]
    disclaimer: str = "Educational prototype only. Not a medical diagnosis or treatment recommendation."


class ClinicalHUD:
    def build_report(
        self,
        title: str,
        findings: Iterable[str],
        evidence: Iterable[str],
    ) -> ClinicalReport:
        return ClinicalReport(title, tuple(findings), tuple(evidence))
