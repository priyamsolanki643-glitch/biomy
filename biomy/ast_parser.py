import re
from typing import List, Tuple

from .models import AnalyteNode, ExposureNode


class BioASTParser:
    """Parse typed laboratory text and prescriptions into typed BioSymbols."""

    LAB_PATTERN = re.compile(
        r"(?P<name>[A-Za-z][A-Za-z0-9 _-]*?)"
        r"\s*[:=]?\s*(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+))"
        r"\s*(?P<unit>%|mg/dL|mg\/dL|mg/L|mg\/L|mmol/L|mmol\/L|mmol/mol|"
        r"umol/L|umol\/L|micromol/L|micromol\/L|mL/min|IU/L|g/dL|mcg/mL|ng/mL)",
        re.IGNORECASE,
    )

    EXPOSURE_PATTERN = re.compile(
        r"\b(?P<name>[A-Za-z][A-Za-z ]{2,30}?)\s+"
        r"(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*"
        r"(?P<unit>mg|mcg|µg|g|mL)\b",
        re.IGNORECASE,
    )

    def parse(self, text: str) -> Tuple[List[AnalyteNode], List[ExposureNode]]:
        analytes: List[AnalyteNode] = []
        exposures: List[ExposureNode] = []

        for match in self.LAB_PATTERN.finditer(text):
            name = re.sub(r"\s+", " ", match.group("name").strip())
            analytes.append(
                AnalyteNode(
                    name=name,
                    value=float(match.group("value")),
                    unit=match.group("unit"),
                    source_text=match.group(0),
                )
            )

        for match in self.EXPOSURE_PATTERN.finditer(text):
            exposures.append(
                ExposureNode(
                    name=re.sub(r"\s+", " ", match.group("name").strip()),
                    kind="exposure",
                    strength=float(match.group("value")),
                    unit=match.group("unit"),
                    source_text=match.group(0),
                )
            )

        return analytes, exposures
