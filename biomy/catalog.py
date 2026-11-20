import json
from pathlib import Path
from typing import Any

from .tpsm import TemporalDrug, TemporalHerb


CATALOG_PATH = Path(__file__).parent.parent / "data" / "temporal_catalog.json"


def load_temporal_catalog(path: Path = CATALOG_PATH) -> tuple[dict[str, TemporalDrug], dict[str, TemporalHerb], dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    drugs = {
        item["name"].lower(): TemporalDrug(
            name=item["name"],
            bioavailability=item["bioavailability"],
            absorption_rate_h=item["absorption_rate_h"],
            volume_l=item["volume_l"],
            clearance_l_h=item["clearance_l_h"],
            therapeutic_min_mg_l=item.get("therapeutic_min_mg_l"),
            toxicity_threshold_mg_l=item.get("toxicity_threshold_mg_l"),
            half_life_h=item.get("half_life_h"),
            age_restriction=item.get("age_restriction", ""),
            side_effects=tuple(item.get("side_effects", ())),
        )
        for item in payload.get("drugs", ())
    }
    herbs = {
        item["name"].lower(): TemporalHerb(
            name=item["name"],
            bioavailability=item["bioavailability"],
            absorption_rate_h=item["absorption_rate_h"],
            volume_l=item["volume_l"],
            clearance_l_h=item["clearance_l_h"],
            enzyme=item.get("enzyme"),
            inhibition_strength=item.get("inhibition_strength", 0.0),
            enzyme_recovery_half_life_h=item.get("enzyme_recovery_half_life_h"),
            onset_h=item.get("onset_h"),
            duration_h=item.get("duration_h"),
            age_restriction=item.get("age_restriction", ""),
            side_effects=tuple(item.get("side_effects", ())),
        )
        for item in payload.get("herbs", ())
    }
    return drugs, herbs, payload.get("metadata", {})
