"""Regression and dynamic physiological recovery tests."""

import numpy as np
import pytest
from biomy.tpsm.engine import TPSMEngine
from biomy.tpsm.events import DoseEvent
from biomy.tpsm.models import (
    DrugParameters,
    EnzymePathwayParam,
    EvidenceLevel,
    HerbParameters,
    ParameterStatus,
    PathwayType,
    ProvenanceRecord,
    TPSMValidationError,
)
from biomy.tpsm.parameters import ParameterCatalog
from biomy.tpsm.serialization import TPSMSerializer


@pytest.fixture
def mbi_system():
    drug = DrugParameters(
        name="SensitiveDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        fm_pathways={"CYP3A4": 0.8},
        provenance=ProvenanceRecord(
            source_reference="Clinical PK Study 2024",
            evidence_level=EvidenceLevel.IN_VIVO_HUMAN,
            parameter_status=ParameterStatus.SOURCED,
        ),
    )
    herb = HerbParameters(
        name="MBIHerb",
        bioavailability=1.0,
        ka_h=2.0,
        volume_central_l=20.0,
        clearance_l_h=10.0,  # Fast clearance (clears within 12h)
        pathway_effects={
            "CYP3A4": EnzymePathwayParam(
                pathway_name="CYP3A4",
                pathway_type=PathwayType.CYP_ENZYME,
                kinact_h=1.5,
                ki_u_mbi_mg_l=0.2,
                recovery_half_life_h=36.0,  # Slow enzyme resynthesis (36h half-life)
            )
        },
    )
    return drug, herb


def test_enzyme_does_not_instantly_recover_when_inhibitor_stops(mbi_system):
    drug, herb = mbi_system
    engine = TPSMEngine([drug], [herb])
    # Single dose of herb at t=0
    res = engine.simulate([DoseEvent("MBIHerb", 50.0, 0.0)], duration_h=72.0, step_h=1.0)

    times = res.time_points
    enz_curve = res.enzyme_activity_curves["CYP3A4"]
    herb_curve = res.concentration_curves[herb.name]

    # At t=12h, herb is mostly cleared
    idx_12h = times.index(12.0) if 12.0 in times else 12
    assert herb_curve[idx_12h] < 0.2

    # But enzyme at t=12h is still severely depleted (< 60%)
    assert enz_curve[idx_12h] < 0.6


def test_enzyme_gradually_recovers_after_inhibitor_discontinuation(mbi_system):
    drug, herb = mbi_system
    engine = TPSMEngine([drug], [herb])
    res = engine.simulate([DoseEvent("MBIHerb", 50.0, 0.0)], duration_h=120.0, step_h=1.0)

    times = res.time_points
    enz_curve = res.enzyme_activity_curves["CYP3A4"]

    # Nadir (minimum active enzyme)
    min_enz = min(enz_curve)
    nadir_idx = enz_curve.index(min_enz)

    # By t=120h, enzyme should have recovered substantially from nadir
    final_enz = enz_curve[-1]
    assert final_enz > min_enz
    assert final_enz > 0.80  # Approaches baseline after ~3 half-lives


def test_mbi_produces_delayed_victim_washout(mbi_system):
    drug, herb = mbi_system
    # Run baseline drug without herb
    engine_baseline = TPSMEngine([drug])
    res_base = engine_baseline.simulate([DoseEvent("SensitiveDrug", 20.0, 24.0)], duration_h=48.0)

    # Run drug with herb given at t=0 (herb clears by t=24h, but enzyme remains suppressed)
    engine_co = TPSMEngine([drug], [herb])
    res_co = engine_co.simulate(
        [
            DoseEvent("MBIHerb", 50.0, 0.0),
            DoseEvent("SensitiveDrug", 20.0, 24.0),
        ],
        duration_h=48.0,
    )

    # Drug AUC under MBI condition must be significantly higher due to delayed clearance
    auc_base = res_base.auc[drug.name]
    auc_co = res_co.auc[drug.name]
    assert auc_co > auc_base * 1.3


def test_parameter_provenance_is_included_in_output(mbi_system):
    drug, herb = mbi_system
    engine = TPSMEngine([drug], [herb])
    res = engine.simulate([], duration_h=10.0)

    prov = res.parameter_provenance
    assert drug.name in prov
    assert prov[drug.name]["source"] == "Clinical PK Study 2024"
    assert prov[drug.name]["evidence_level"] == "in_vivo_human"


def test_biomistral_serializer_structure(mbi_system):
    drug, herb = mbi_system
    engine = TPSMEngine([drug], [herb])
    res = engine.simulate([DoseEvent("SensitiveDrug", 10.0, 0.0)], duration_h=24.0)

    context = TPSMSerializer.to_biomistral_context(res)
    assert "system_instruction_for_llm" in context
    assert "DO NOT recalculate" in context["system_instruction_for_llm"]
    assert drug.name in context["exposure_metrics"]["auc_mg_h_l"]

    prompt = TPSMSerializer.to_formatted_prompt(res)
    assert "BIOMY DETERMINISTIC TPSM FACT MATRIX" in prompt


def test_dose_at_duration_h_is_applied_and_logged():
    """Verify doses scheduled at t = duration_h are correctly applied and logged."""
    drug = DrugParameters(
        name="BoundaryDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
    )
    engine = TPSMEngine([drug])
    res = engine.simulate([DoseEvent("BoundaryDrug", 100.0, 24.0)], duration_h=24.0)

    # Event log must record the dose at t=24h
    assert len(res.event_log) == 1
    assert res.event_log[0]["time_h"] == 24.0
    assert res.event_log[0]["dose_mg"] == 100.0
    # The final point in time_points must be 24.0
    assert res.time_points[-1] == 24.0


def test_genetic_baseline_remains_stable_without_inhibitor():
    """Verify that genetic poor metabolizer baselines do not spontaneously revert to wild-type."""
    drug = DrugParameters(
        name="PolymorphicDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        fm_pathways={"CYP2C9": 0.8},
    )
    # Set genetic baseline to 50% capacity (e.g. *2/*3 intermediate/poor metabolizer)
    engine = TPSMEngine([drug], enzyme_baseline={"CYP2C9": 0.5})
    res = engine.simulate([], duration_h=168.0, step_h=1.0)

    cyp2c9_curve = res.enzyme_activity_curves["CYP2C9"]
    # The enzyme activity must stay at 0.5 across the entire 168-hour horizon
    for val in cyp2c9_curve:
        assert np.isclose(val, 0.5, atol=1e-4)


def test_catalog_enzyme_recovery_halflife_used_from_herb():
    """Verify that catalog-defined enzyme recovery half-lives are prioritized over the default."""
    catalog = ParameterCatalog.from_json_file()
    drug = catalog.get_drug("Warfarin")
    herb = catalog.get_herb("Piperine")
    engine = TPSMEngine([drug], [herb])

    # Piperine in catalog defines CYP3A4 recovery_half_life_h = 28.9h
    expected_kdeg_cyp3a4 = np.log(2.0) / 28.9
    assert np.isclose(engine.enzyme_kdeg["CYP3A4"], expected_kdeg_cyp3a4, rtol=1e-4)

    # Piperine in catalog defines CYP2C9 recovery_half_life_h = 38.0h
    expected_kdeg_cyp2c9 = np.log(2.0) / 38.0
    assert np.isclose(engine.enzyme_kdeg["CYP2C9"], expected_kdeg_cyp2c9, rtol=1e-4)


def test_degenerate_interval_solver_resilience():
    """Verify that degenerate and closely spaced events do not cause IndexError or solver crashes."""
    drug = DrugParameters(
        name="StiffDrug",
        bioavailability=1.0,
        ka_h=2.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
    )
    engine = TPSMEngine([drug])
    # Doses spaced by zero-duration or tiny delta
    events = [
        DoseEvent("StiffDrug", 50.0, 0.0),
        DoseEvent("StiffDrug", 50.0, 12.0),
        DoseEvent("StiffDrug", 50.0, 12.0),  # Exact duplicate timestamp
        DoseEvent("StiffDrug", 50.0, 24.0),
    ]
    res = engine.simulate(events, duration_h=24.0, step_h=0.5)
    assert res.solver_success is True
    assert len(res.time_points) >= 2


def test_validation_rejects_inverted_therapeutic_thresholds():
    """Verify validation rejects therapeutic_min_mg_l >= toxicity_threshold_mg_l."""
    drug = DrugParameters(
        name="InvertedDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        therapeutic_min_mg_l=20.0,
        toxicity_threshold_mg_l=10.0,  # Inverted!
    )
    with pytest.raises(TPSMValidationError, match="strictly less than"):
        TPSMEngine([drug])


def test_validation_rejects_negative_fm_pathways():
    """Verify validation rejects negative individual fm pathways."""
    drug = DrugParameters(
        name="NegativeFmDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        fm_pathways={"CYP3A4": -0.2, "CYP2C9": 0.8},
    )
    with pytest.raises(TPSMValidationError, match="must be non-negative"):
        TPSMEngine([drug])


def test_validation_rejects_duplicate_compound_names():
    """Verify validation rejects duplicate compound names in drugs and herbs."""
    drug1 = DrugParameters(name="DuplicateName", volume_central_l=10.0, clearance_l_h=1.0)
    drug2 = DrugParameters(name="DuplicateName", volume_central_l=20.0, clearance_l_h=2.0)
    with pytest.raises(TPSMValidationError, match="Duplicate compound name"):
        TPSMEngine([drug1, drug2])


def test_simulation_rejects_invalid_step_h():
    """Verify engine.simulate rejects negative, zero, or excessive step_h."""
    drug = DrugParameters(name="StepDrug", volume_central_l=10.0, clearance_l_h=1.0)
    engine = TPSMEngine([drug])
    with pytest.raises(TPSMValidationError, match="step_h must be finite and positive"):
        engine.simulate([], duration_h=24.0, step_h=-0.5)
    with pytest.raises(TPSMValidationError, match="step_h must be finite and positive"):
        engine.simulate([], duration_h=24.0, step_h=0.0)
    with pytest.raises(TPSMValidationError, match="cannot exceed duration_h"):
        engine.simulate([], duration_h=24.0, step_h=48.0)
