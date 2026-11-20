"""Catalog loader, parameter provenance mapper, and registry manager."""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .models import (
    DrugParameters,
    EliminationMode,
    EnzymePathwayParam,
    EvidenceLevel,
    HerbParameters,
    ParameterStatus,
    PathwayType,
    ProvenanceRecord,
    TPSMValidationError,
    UncertaintyRange,
)


DEFAULT_CATALOG_PATH = Path(__file__).parent.parent.parent / "data" / "temporal_catalog.json"


class ParameterCatalog:
    """Manages typed drug and botanical parameters with source provenance and audit metadata."""

    def __init__(
        self,
        drugs: Dict[str, DrugParameters],
        herbs: Dict[str, HerbParameters],
        metadata: Optional[Dict[str, Any]] = None,
    ):
        self.drugs = drugs
        self.herbs = herbs
        self.metadata = metadata or {}

    @classmethod
    def from_json_file(cls, path: Any = DEFAULT_CATALOG_PATH) -> "ParameterCatalog":
        """Load and parse typed catalog from JSON file."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Temporal catalog not found at {path}")

        with path.open(encoding="utf-8") as handle:
            raw = json.load(handle)

        drugs: Dict[str, DrugParameters] = {}
        for item in raw.get("drugs", []):
            prov_raw = item.get("provenance", {})
            unc_raw = prov_raw.get("uncertainty_range")
            unc = UncertaintyRange(**unc_raw) if unc_raw else None
            prov = ProvenanceRecord(
                source_reference=prov_raw.get("source_reference", item.get("source_reference", "Clinical standard")),
                evidence_level=EvidenceLevel(prov_raw.get("evidence_level", "in_vivo_human")),
                parameter_status=ParameterStatus(prov_raw.get("parameter_status", "sourced")),
                uncertainty_range=unc,
                preparation=prov_raw.get("preparation", item.get("formulation", "oral")),
                population_assumptions=prov_raw.get("population_assumptions", "Adult 70kg reference"),
                limitations=prov_raw.get("limitations", ""),
            )

            elim_mode = EliminationMode.MICHAELIS_MENTEN if item.get("nonlinear") else EliminationMode.LINEAR

            drug = DrugParameters(
                name=item["name"],
                role=item.get("role", "victim"),
                formulation=item.get("formulation", "oral_solid"),
                bioavailability=float(item["bioavailability"]),
                ka_h=float(item.get("ka_h", item.get("absorption_rate_h", 1.0))),
                volume_central_l=float(item.get("volume_central_l", item.get("volume_l", 10.0))),
                volume_peripheral_l=float(item["volume_peripheral_l"]) if item.get("volume_peripheral_l") else None,
                inter_compartment_clearance_l_h=float(item["inter_compartment_clearance_l_h"]) if item.get("inter_compartment_clearance_l_h") else None,
                clearance_l_h=float(item.get("clearance_l_h", 1.0)),
                fraction_unbound=float(item.get("fraction_unbound", 1.0)),
                half_life_h=float(item["half_life_h"]) if item.get("half_life_h") else None,
                elimination_mode=elim_mode,
                vmax_mg_h=float(item["vmax_mg_h"]) if item.get("vmax_mg_h") else None,
                km_mg_l=float(item["km_mg_l"]) if item.get("km_mg_l") else None,
                fm_pathways=item.get("fm_pathways", {}),
                pathway_types=item.get("pathway_types", {}),
                therapeutic_min_mg_l=float(item["therapeutic_min_mg_l"]) if item.get("therapeutic_min_mg_l") is not None else None,
                toxicity_threshold_mg_l=float(item["toxicity_threshold_mg_l"]) if item.get("toxicity_threshold_mg_l") is not None else None,
                age_restrictions=item.get("age_restrictions", item.get("age_restriction", "Requires clinical review")),
                pregnancy_restrictions=item.get("pregnancy_restrictions", ""),
                renal_hepatic_cautions=item.get("renal_hepatic_cautions", ""),
                side_effects=tuple(item.get("side_effects", ())),
                provenance=prov,
            )
            drugs[drug.name.lower()] = drug

        herbs: Dict[str, HerbParameters] = {}
        for item in raw.get("herbs", []):
            prov_raw = item.get("provenance", {})
            unc_raw = prov_raw.get("uncertainty_range")
            unc = UncertaintyRange(**unc_raw) if unc_raw else None
            prov = ProvenanceRecord(
                source_reference=prov_raw.get("source_reference", item.get("source_reference", "Literature monograph")),
                evidence_level=EvidenceLevel(prov_raw.get("evidence_level", "in_vitro")),
                parameter_status=ParameterStatus(prov_raw.get("parameter_status", "estimated")),
                uncertainty_range=unc,
                preparation=prov_raw.get("preparation", item.get("preparation", "standardized_extract")),
                population_assumptions=prov_raw.get("population_assumptions", "Healthy adult"),
                limitations=prov_raw.get("limitations", "Concentration and formulation dependent"),
            )

            # Build pathway effects
            pathways: Dict[str, EnzymePathwayParam] = {}
            raw_pathways = item.get("pathway_effects", {})
            if not raw_pathways and item.get("enzyme"):
                enz = item["enzyme"]
                pathways[enz] = EnzymePathwayParam(
                    pathway_name=enz,
                    pathway_type=PathwayType.CYP_ENZYME,
                    reversible_ki_u_mg_l=item.get("reversible_ki_u_mg_l"),
                    kinact_h=item.get("kinact_h"),
                    ki_u_mbi_mg_l=item.get("ki_u_mbi_mg_l"),
                    recovery_half_life_h=item.get("enzyme_recovery_half_life_h", 36.0),
                )
            else:
                for enz_name, p_data in raw_pathways.items():
                    pathways[enz_name] = EnzymePathwayParam(
                        pathway_name=enz_name,
                        pathway_type=PathwayType(p_data.get("pathway_type", "cyp_enzyme")),
                        is_transporter=p_data.get("is_transporter", False),
                        reversible_ki_u_mg_l=p_data.get("reversible_ki_u_mg_l"),
                        kinact_h=p_data.get("kinact_h"),
                        ki_u_mbi_mg_l=p_data.get("ki_u_mbi_mg_l"),
                        recovery_half_life_h=p_data.get("recovery_half_life_h", 36.0),
                    )

            herb = HerbParameters(
                name=item["name"],
                preparation=item.get("preparation", "standardized_extract"),
                active_constituent=item.get("active_constituent", ""),
                dose_unit=item.get("dose_unit", "mg"),
                bioavailability=float(item.get("bioavailability", 0.5)),
                ka_h=float(item.get("ka_h", item.get("absorption_rate_h", 1.0))),
                volume_central_l=float(item.get("volume_central_l", item.get("volume_l", 50.0))),
                clearance_l_h=float(item.get("clearance_l_h", 10.0)),
                fraction_unbound=float(item.get("fraction_unbound", 0.1)),
                pathway_effects=pathways,
                onset_h=float(item["onset_h"]) if item.get("onset_h") else None,
                duration_h=float(item["duration_h"]) if item.get("duration_h") else None,
                age_restrictions=item.get("age_restrictions", item.get("age_restriction", "Requires clinical review")),
                pregnancy_restrictions=item.get("pregnancy_restrictions", ""),
                side_effects=tuple(item.get("side_effects", ())),
                provenance=prov,
            )
            herbs[herb.name.lower()] = herb

        return cls(drugs=drugs, herbs=herbs, metadata=raw.get("metadata", {}))

    def get_drug(self, name: str) -> DrugParameters:
        key = name.strip().lower()
        if key not in self.drugs:
            raise TPSMValidationError(f"Drug '{name}' not found in catalog")
        return self.drugs[key]

    def get_herb(self, name: str) -> HerbParameters:
        key = name.strip().lower()
        if key not in self.herbs:
            raise TPSMValidationError(f"Herb '{name}' not found in catalog")
        return self.herbs[key]
