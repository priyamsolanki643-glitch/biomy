"""Unit tests for TPSM models, parameter validation, and provenance records."""

import pytest
from biomy.tpsm.models import (
    CompartmentMode,
    DrugParameters,
    EliminationMode,
    EnzymePathwayParam,
    EvidenceLevel,
    HerbParameters,
    ParameterStatus,
    PathwayType,
    ProvenanceRecord,
    SimulationResult,
    TPSMValidationError,
    UncertaintyRange,
)
from biomy.tpsm.parameters import ParameterCatalog
from biomy.tpsm.validation import SimulationValidator


def test_drug_parameter_validation_catches_invalid_bioavailability():
    drug = DrugParameters(
        name="BadDrug",
        bioavailability=1.5,  # Invalid: must be in (0, 1]
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
    )
    with pytest.raises(TPSMValidationError, match="bioavailability"):
        SimulationValidator.validate_drug(drug, CompartmentMode.ONE_COMPARTMENT)


def test_drug_parameter_validation_catches_invalid_volume():
    drug = DrugParameters(
        name="BadDrug",
        bioavailability=0.9,
        ka_h=1.0,
        volume_central_l=-5.0,  # Invalid: must be positive
        clearance_l_h=1.0,
    )
    with pytest.raises(TPSMValidationError, match="volume"):
        SimulationValidator.validate_drug(drug, CompartmentMode.ONE_COMPARTMENT)


def test_drug_parameter_validation_catches_michaelis_menten_missing_km():
    drug = DrugParameters(
        name="BadDrug",
        bioavailability=0.9,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        elimination_mode=EliminationMode.MICHAELIS_MENTEN,
        vmax_mg_h=10.0,
        km_mg_l=None,  # Missing Km
    )
    with pytest.raises(TPSMValidationError, match="Km must be strictly positive"):
        SimulationValidator.validate_drug(drug, CompartmentMode.ONE_COMPARTMENT)


def test_herb_parameter_validation_catches_negative_kinact():
    herb = HerbParameters(
        name="BadHerb",
        bioavailability=0.5,
        ka_h=1.0,
        volume_central_l=20.0,
        clearance_l_h=5.0,
        pathway_effects={
            "CYP3A4": EnzymePathwayParam(
                pathway_name="CYP3A4",
                kinact_h=-0.5,  # Invalid: negative kinact
                ki_u_mbi_mg_l=1.0,
            )
        },
    )
    with pytest.raises(TPSMValidationError, match="kinact must be non-negative"):
        SimulationValidator.validate_herb(herb)


def test_provenance_record_instantiation_and_serialization():
    prov = ProvenanceRecord(
        source_reference="Clinical Study 2024",
        evidence_level=EvidenceLevel.IN_VIVO_HUMAN,
        parameter_status=ParameterStatus.SOURCED,
        uncertainty_range=UncertaintyRange(low=10.0, high=15.0, unit="L/h", confidence="high"),
        preparation="tablet",
        population_assumptions="Adult",
        limitations="None",
    )
    assert prov.source_reference == "Clinical Study 2024"
    assert prov.evidence_level == EvidenceLevel.IN_VIVO_HUMAN
    assert prov.uncertainty_range.low == 10.0
