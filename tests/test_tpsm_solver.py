"""Unit tests for SciPy ODE solver integration, stiffness handling, and numerical stability."""

import numpy as np
import pytest
from biomy.tpsm.engine import TPSMEngine
from biomy.tpsm.events import DoseEvent
from biomy.tpsm.models import (
    CompartmentMode,
    DrugParameters,
    TPSMSolverError,
)
from biomy.tpsm.solver import ODESolverWrapper


@pytest.fixture
def standard_drug():
    return DrugParameters(
        name="StandardDrug",
        bioavailability=1.0,
        ka_h=1.5,
        volume_central_l=10.0,
        volume_peripheral_l=5.0,
        inter_compartment_clearance_l_h=2.0,
        clearance_l_h=2.0,
        toxicity_threshold_mg_l=10.0,
    )


def test_zero_dose_state_remains_zero(standard_drug):
    engine = TPSMEngine([standard_drug])
    res = engine.simulate([], duration_h=24.0)
    curve = res.concentration_curves[standard_drug.name]
    assert all(c == 0.0 for c in curve)
    assert res.auc[standard_drug.name] == 0.0


def test_one_oral_dose_produces_nonzero_curve_and_decays(standard_drug):
    engine = TPSMEngine([standard_drug])
    res = engine.simulate([DoseEvent("StandardDrug", 100.0, 0.0)], duration_h=24.0)
    curve = res.concentration_curves[standard_drug.name]

    # Starts at 0, rises to peak, then decays
    assert curve[0] == 0.0
    peak_idx = int(np.argmax(curve))
    assert peak_idx > 0
    assert curve[peak_idx] > 0.0
    assert curve[-1] < curve[peak_idx]  # Decays after peak


def test_solver_failure_is_explicit():
    wrapper = ODESolverWrapper(method="Radau")

    # Define a pathological ODE that blows up instantly
    def exploding_fun(t, y):
        return 1e12 * y ** 2

    with pytest.raises(TPSMSolverError):
        wrapper.integrate_segment(exploding_fun, (0.0, 10.0), np.array([100.0]))


def test_non_finite_state_is_rejected():
    wrapper = ODESolverWrapper(method="Radau")

    def dummy_fun(t, y):
        return np.zeros_like(y)

    with pytest.raises(TPSMSolverError, match="non-finite"):
        wrapper.integrate_segment(dummy_fun, (0.0, 1.0), np.array([np.nan]))


def test_two_compartment_mass_remains_numerically_non_negative(standard_drug):
    engine = TPSMEngine([standard_drug], compartment_mode=CompartmentMode.TWO_COMPARTMENT)
    res = engine.simulate([DoseEvent("StandardDrug", 50.0, 0.0)], duration_h=48.0)
    curve = res.concentration_curves[standard_drug.name]
    assert all(c >= 0.0 for c in curve)
    assert min(curve) >= 0.0


@pytest.mark.parametrize("method", ["Radau", "BDF", "RK45"])
def test_solver_methods_execute_successfully(standard_drug, method):
    engine = TPSMEngine([standard_drug], solver_method=method)
    res = engine.simulate([DoseEvent("StandardDrug", 25.0, 0.0)], duration_h=12.0)
    assert res.solver_success
    assert res.c_max[standard_drug.name] > 0.0
