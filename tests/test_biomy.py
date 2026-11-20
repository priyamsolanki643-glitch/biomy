import math
import pytest

from biomy.ast_parser import BioASTParser
from biomy.calculations import BioCalculations
from biomy.ripple import ClinicalRippleEngine, Interaction
from biomy.scalar import ContinuousScalarProjector
from biomy.verification import BiochemicalInvariantGate, ClearanceParityChecker
from biomy.hud import ClinicalHUD
from biomy.catalog import load_temporal_catalog
from biomy.tpsm import DoseEvent, TemporalPKStateMachine


def test_parser_extracts_analytes_and_exposures():
    text = "Creatinine: 1.2 mg/dL; HbA1c: 48 mmol/mol; Warfarin 5 mg"
    analytes, exposures = BioASTParser().parse(text)
    assert [(node.name.lower(), node.value) for node in analytes] == [
        ("creatinine", 1.2),
        ("hba1c", 48.0),
    ]
    assert exposures[0].name.lower() == "warfarin"
    assert exposures[0].strength == 5.0


def test_parser_handles_no_space_and_frequency():
    analytes, exposures = BioASTParser().parse("Tab. Warfarin 5mg OD")
    assert len(analytes) == 0
    assert len(exposures) == 1
    assert exposures[0].name.lower() == "warfarin"
    assert exposures[0].strength == 5.0


def test_parser_keeps_frequency_with_its_own_line():
    _, exposures = BioASTParser().parse(
        "Warfarin 5 mg OD\nMidazolam 2 mg BD"
    )
    assert len(exposures) == 2
    assert exposures[0].name.lower() == "warfarin"
    assert exposures[0].route == "OD"
    assert exposures[1].name.lower() == "midazolam"
    assert exposures[1].route == "BD"


def test_parser_does_not_treat_lab_value_as_exposure():
    analytes, exposures = BioASTParser().parse("Creatinine: 1.2 mg/dL")
    assert len(analytes) == 1
    assert len(exposures) == 0


def test_ckd_epi_returns_finite_positive_result():
    result = BioCalculations.egfr_ckd_epi_2021(1.2, 50, False)
    assert math.isfinite(result)
    assert result > 0


@pytest.mark.parametrize("age", [1, 141])
def test_ckd_epi_rejects_invalid_adult_age(age):
    with pytest.raises(ValueError):
        BioCalculations.egfr_ckd_epi_2021(1.0, age, False)


def test_r_ratio_rejects_zero_alp():
    with pytest.raises(ValueError):
        BioCalculations.hepatocyte_r_ratio(40, 40, 0, 120)


def test_ripple_frontier_finds_transitive_interactions():
    engine = ClinicalRippleEngine()
    engine.add_interaction(Interaction("Piperine", "CYP3A4", "CYP3A4", "inhibition", "moderate", "REF-1"))
    engine.add_interaction(Interaction("CYP3A4", "Midazolam", "CYP3A4", "clearance decrease", "high", "REF-2"))
    frontier = engine.ripple_frontier(["Piperine"])
    assert any(record["node"] == "Midazolam" for record in frontier["Piperine"])


def test_pair_screen_does_not_attach_wrong_herb_evidence():
    engine = ClinicalRippleEngine()
    engine.add_interaction(
        Interaction(
            "St. John's wort",
            "Cyclosporine",
            "CYP3A4",
            "may decrease drug exposure",
            "high",
            "STJOHN-CYCLOSPORINE-REF",
        )
    )
    result = engine.screen_pair("Piperine", "Cyclosporine")
    assert result["status"] == "no_interaction_found"


def test_scalar_projector_has_equal_scale_properties():
    projector = ContinuousScalarProjector(harmonics=4)
    assert len(projector.project(1.8, 10)) == 8
    with_value_error = False
    try:
        projector.project(float("nan"))
    except ValueError:
        with_value_error = True
    assert with_value_error


def test_verification_rejects_unapproved_claim():
    result = BiochemicalInvariantGate().verify(
        ["verified CYP3A4 inhibition", "made-up enzyme"],
        ["verified CYP3A4 inhibition"],
    )
    assert not result.passed
    assert "made-up enzyme" in result.violations


def test_renal_dose_check():
    result = ClearanceParityChecker().verify("Metformin", 1000, 25, 500)
    assert not result.passed


def test_hud_has_safety_disclaimer():
    report = ClinicalHUD().build_report("Test", ["Finding"], ["Reference"])
    assert "Educational prototype" in report.disclaimer


def test_temporal_catalog_contains_common_starter_entries():
    drugs, herbs, metadata = load_temporal_catalog()
    assert "paracetamol" in drugs
    assert "warfarin" in drugs
    assert "piperine" in herbs
    assert metadata["status"] == "educational starter data"


def test_tpsm_tracks_concentration_and_washout():
    drugs, herbs, _ = load_temporal_catalog()
    result = TemporalPKStateMachine([drugs["warfarin"]], [herbs["piperine"]]).simulate(
        [
            DoseEvent("Warfarin", 5.0, 0.0),
            DoseEvent("Piperine", 20.0, 0.0),
        ],
        duration_h=48,
        step_h=1.0,
    )
    series = result.series("Warfarin")
    assert len(series) == 49
    assert series[1][1] > 0
    assert series[-1][1] < series[1][1]
    assert result.points[0].active_enzyme_fraction["CYP3A4"] == 1.0


def test_tpsm_rejects_unknown_event_and_invalid_window():
    drugs, _, _ = load_temporal_catalog()
    engine = TemporalPKStateMachine([drugs["warfarin"]])
    with pytest.raises(ValueError):
        engine.simulate([DoseEvent("Unknown", 1.0, 0.0)], 24)
    with pytest.raises(ValueError):
        engine.simulate([], 0)
