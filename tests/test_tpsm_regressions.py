"""Regression and dynamic physiological recovery tests."""

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
)
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
