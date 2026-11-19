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
        r"\b(?P<name>(?:[A-Za-z][A-Za-z0-9-]*(?:\s+[A-Za-z][A-Za-z0-9-]*){0,4})?)"
        r"\s*(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+))"
        r"\s*(?P<unit>mcg|µg|mg|g|mL)\b",
        re.IGNORECASE,
    )

    FREQUENCY_PATTERN = re.compile(
        r"\b(OD|BD|BID|TDS|TID|QDS|QID|QAM|QPM|HS|SOS|PRN|once daily|twice daily|"
        r"three times daily|every \d+\s*(?:hour|hours|h))\b",
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
            name = re.sub(r"\s+", " ", match.group("name").strip())
            name = re.sub(r"\b(?:tab(?:let)?|tab\.|caps(?:ule)?|cap\.|syrup|oral)\b\s*", "", name, flags=re.IGNORECASE).strip()
            if not name:
                continue
            nearby = text[max(0, match.start() - 50):min(len(text), match.end() + 50)]
            frequency_match = self.FREQUENCY_PATTERN.search(nearby)
            exposures.append(
                ExposureNode(
                    name=name,
                    kind="exposure",
                    strength=float(match.group("value")),
                    unit=match.group("unit"),
                    route=frequency_match.group(1).upper() if frequency_match else "",
                    source_text=match.group(0),
                )
            )

        return analytes, exposures
