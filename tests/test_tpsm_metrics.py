"""Unit tests for exposure metrics, AUC calculation, TMR behavior, and toxicity root-finding."""

import numpy as np
import pytest
from biomy.tpsm.engine import TPSMEngine
from biomy.tpsm.events import DoseEvent
from biomy.tpsm.metrics import PKMetricsCalculator
from biomy.tpsm.models import DrugParameters


@pytest.fixture
def metric_drug():
    return DrugParameters(
        name="MetricDrug",
        bioavailability=1.0,
        ka_h=2.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        therapeutic_min_mg_l=1.0,
        toxicity_threshold_mg_l=5.0,
    )


def test_auc_is_finite_and_positive_after_exposure(metric_drug):
    engine = TPSMEngine([metric_drug])
    res = engine.simulate([DoseEvent("MetricDrug", 20.0, 0.0)], duration_h=24.0)
    auc_val = res.auc[metric_drug.name]
    assert np.isfinite(auc_val)
    assert auc_val > 0.0


def test_toxicity_event_reports_first_upward_crossing(metric_drug):
    engine = TPSMEngine([metric_drug])
    # Dose high enough to cross 5.0 mg/L (e.g. 100 mg in 10 L -> peak ~8 mg/L)
    res = engine.simulate([DoseEvent("MetricDrug", 100.0, 0.0)], duration_h=24.0, step_h=0.1)
    t_breach = res.first_toxicity_time(metric_drug.name)
    assert t_breach is not None
    assert 0.0 < t_breach < 2.0  # Rises rapidly in early hours


def test_no_toxic_crossing_returns_none(metric_drug):
    engine = TPSMEngine([metric_drug])
    # Very small dose (1 mg in 10 L -> peak ~0.08 mg/L << 5.0)
    res = engine.simulate([DoseEvent("MetricDrug", 1.0, 0.0)], duration_h=24.0)
    t_breach = res.first_toxicity_time(metric_drug.name)
    assert t_breach is None


def test_tmr_decreases_as_concentration_approaches_toxicity_threshold(metric_drug):
    engine = TPSMEngine([metric_drug])
    res = engine.simulate([DoseEvent("MetricDrug", 40.0, 0.0)], duration_h=12.0, step_h=0.2)
    tmr_series = res.therapeutic_margin_ratio_curves[metric_drug.name]
    conc_series = res.concentration_curves[metric_drug.name]

    # At t=0, C=0 -> TMR = (5 - 0) / (5 - 1) = 1.25
    assert tmr_series[0] == 1.25

    # Find point of peak concentration
    peak_idx = int(np.argmax(conc_series))
    assert peak_idx > 0
    # As concentration rises, TMR must decrease
    assert tmr_series[peak_idx] < tmr_series[0]


def test_missing_toxicity_threshold_yields_none_tmr():
    drug_no_tox = DrugParameters(
        name="BenignDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        toxicity_threshold_mg_l=None,
    )
    engine = TPSMEngine([drug_no_tox])
    res = engine.simulate([DoseEvent("BenignDrug", 10.0, 0.0)], duration_h=10.0)
    tmr_series = res.therapeutic_margin_ratio_curves[drug_no_tox.name]
    assert all(val is None for val in tmr_series)
