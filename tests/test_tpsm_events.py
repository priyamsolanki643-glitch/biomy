"""Unit tests for discrete dose event scheduling, timestamp precision, and input validation."""

import pytest
from biomy.tpsm.engine import TPSMEngine
from biomy.tpsm.events import DoseEvent
from biomy.tpsm.models import DrugParameters, TPSMValidationError


@pytest.fixture
def test_drug():
    return DrugParameters(
        name="EventDrug",
        bioavailability=1.0,
        ka_h=2.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        toxicity_threshold_mg_l=25.0,
    )


def test_dose_at_t24_is_applied_exactly_at_t24(test_drug):
    engine = TPSMEngine([test_drug])
    res = engine.simulate([DoseEvent("EventDrug", 50.0, 24.0)], duration_h=48.0, step_h=0.5)

    # Before t=24, concentration is 0
    t_arr = res.time_points
    c_arr = res.concentration_curves[test_drug.name]

    for t, c in zip(t_arr, c_arr):
        if t < 24.0:
            assert c == 0.0
        elif t > 25.0:
            assert c > 0.0


def test_repeated_doses_accumulate(test_drug):
    engine = TPSMEngine([test_drug])
    # Single dose vs two doses 6 hours apart
    res_single = engine.simulate([DoseEvent("EventDrug", 50.0, 0.0)], duration_h=24.0)
    res_double = engine.simulate(
        [
            DoseEvent("EventDrug", 50.0, 0.0),
            DoseEvent("EventDrug", 50.0, 6.0),
        ],
        duration_h=24.0,
    )

    c_single_max = res_single.c_max[test_drug.name]
    c_double_max = res_double.c_max[test_drug.name]
    assert c_double_max > c_single_max


def test_no_dose_is_applied_twice(test_drug):
    engine = TPSMEngine([test_drug])
    res = engine.simulate([DoseEvent("EventDrug", 10.0, 0.0)], duration_h=10.0)
    # The event log should record exactly one dose event
    assert len(res.event_log) == 1
    assert res.event_log[0]["dose_mg"] == 10.0


def test_solver_rejects_negative_dose(test_drug):
    engine = TPSMEngine([test_drug])
    with pytest.raises(TPSMValidationError, match="positive"):
        engine.simulate([DoseEvent("EventDrug", -10.0, 0.0)], duration_h=24.0)


def test_solver_rejects_unknown_compound(test_drug):
    engine = TPSMEngine([test_drug])
    with pytest.raises(TPSMValidationError, match="Unknown compound"):
        engine.simulate([DoseEvent("NonExistentCompound", 10.0, 0.0)], duration_h=24.0)


def test_solver_rejects_event_outside_simulation_range(test_drug):
    engine = TPSMEngine([test_drug])
    with pytest.raises(TPSMValidationError, match="exceeds simulation duration"):
        engine.simulate([DoseEvent("EventDrug", 10.0, 30.0)], duration_h=24.0)
