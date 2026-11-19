import re
from typing import List, Tuple

from .models import AnalyteNode, ExposureNode


class BioASTParser:
    """Parse typed laboratory text and prescriptions into typed BioSymbols."""

    LAB_PATTERN = re.compile(
        r"(?P<name>[A-Za-z][A-Za-z0-9 _-]*?)"
        r"\s*[:=]?\s*(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+))"
        r"\s*(?P<unit>%|mg/dL|mg\/dL|mg/L|mg\/L|mmol/L|mmol\/L|mmol/mol|"
        r"umol/L|umol\/L|micromol/L|micromol\/L|mL/min|IU/L|g/dL|mcg/mL|ng/mL)\b",
        re.IGNORECASE,
    )

    EXPOSURE_PATTERN = re.compile(
        r"\b(?P<name>[A-Za-z][A-Za-z0-9-]*(?:\s+[A-Za-z][A-Za-z0-9-]*){0,3})"
        r"\s*(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*"
        r"(?P<unit>mcg|µg|mg|g|mL)\b",
        re.IGNORECASE,
    )

    FREQUENCY_PATTERN = re.compile(
        r"\b(OD|BD|BID|TDS|TID|QDS|QID|QAM|QPM|HS|SOS|PRN|once daily|twice daily|"
        r"three times daily|every \d+\s*(?:hour|hours|h))\b",
        re.IGNORECASE,
    )

    CLEAN_PREFIX_PATTERN = re.compile(
        r"^(?:tab(?:let)?|tab\.|caps(?:ule)?|cap\.|syrup|oral|inj(?:ection)?)\s*",
        re.IGNORECASE,
    )

    def parse(self, text: str) -> Tuple[List[AnalyteNode], List[ExposureNode]]:
        analytes: List[AnalyteNode] = []
        exposures: List[ExposureNode] = []

        # 1. Identify and extract lab analytes, recording their matched character spans
        lab_spans: List[Tuple[int, int]] = []
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
            lab_spans.append(match.span())

        # 2. Parse exposures line-by-line to prevent frequency/dose bleeding across entries
        # Also ensure exposure does not overlap with any lab analyte span
        lines = text.splitlines() if "\n" in text else text.split(";")
        current_offset = 0

        for raw_line in lines:
            line_str = raw_line.strip()
            line_start = text.find(raw_line, current_offset)
            if line_start == -1:
                line_start = current_offset
            current_offset = line_start + len(raw_line)

            for exp_match in self.EXPOSURE_PATTERN.finditer(line_str):
                global_start = line_start + exp_match.start()
                global_end = line_start + exp_match.end()

                # Check if this exposure overlaps with any lab span
                overlaps_lab = any(
                    max(global_start, l_start) < min(global_end, l_end)
                    for l_start, l_end in lab_spans
                )
                if overlaps_lab:
                    continue

                raw_name = exp_match.group("name").strip()
                cleaned_name = self.CLEAN_PREFIX_PATTERN.sub("", raw_name).strip()
                if not cleaned_name:
                    continue

                # Search frequency strictly within the same line
                freq_match = self.FREQUENCY_PATTERN.search(line_str)
                frequency = freq_match.group(1).upper() if freq_match else ""

                exposures.append(
                    ExposureNode(
                        name=cleaned_name,
                        kind="exposure",
                        strength=float(exp_match.group("value")),
                        unit=exp_match.group("unit"),
                        route=frequency,
                        source_text=exp_match.group(0),
                    )
                )

        return analytes, exposures
