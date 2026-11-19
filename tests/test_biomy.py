import math

from biomy.ast_parser import BioASTParser
from biomy.calculations import BioCalculations
from biomy.ripple import ClinicalRippleEngine, Interaction
from biomy.scalar import ContinuousScalarProjector
from biomy.verification import BiochemicalInvariantGate, ClearanceParityChecker
from biomy.hud import ClinicalHUD


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


def test_ckd_epi_returns_finite_positive_result():
    result = BioCalculations.egfr_ckd_epi_2021(1.2, 50, False)
    assert math.isfinite(result)
    assert result > 0


def test_ripple_frontier_finds_transitive_interactions():
    engine = ClinicalRippleEngine()
    engine.add_interaction(Interaction("Piperine", "CYP3A4", "CYP3A4", "inhibition", "moderate", "REF-1"))
    engine.add_interaction(Interaction("CYP3A4", "Midazolam", "CYP3A4", "clearance decrease", "high", "REF-2"))
    frontier = engine.ripple_frontier(["Piperine"])
    assert any(record["node"] == "Midazolam" for record in frontier["Piperine"])


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
