"""Tests for BioMy Temporal Scenarios (data/temporal_scenarios.json).

Validates complete benchmark runs against catalog data and ensures
mechanistic expectations hold across all defined scenarios.
"""

import json
from pathlib import Path
import pytest
import numpy as np

from biomy.tpsm.engine import TPSMEngine
from biomy.tpsm.events import DoseEvent
from biomy.tpsm.models import CompartmentMode
from biomy.tpsm.parameters import ParameterCatalog


@pytest.fixture
def catalog():
    cat_path = Path(__file__).resolve().parent.parent / "data" / "temporal_catalog.json"
    return ParameterCatalog.from_json_file(cat_path)


@pytest.fixture
def scenarios_data():
    scen_path = Path(__file__).resolve().parent.parent / "data" / "temporal_scenarios.json"
    with open(scen_path, "r", encoding="utf-8") as f:
        return json.load(f)["scenarios"]


def test_scenario_a_warfarin_piperine_mbi(catalog):
    """Scenario A: Warfarin + Piperine MBI Recovery Hysteresis."""
    drug = catalog.get_drug("Warfarin")
    herb = catalog.get_herb("Piperine")
    
    engine = TPSMEngine(
        drugs=[drug],
        herbs=[herb],
        compartment_mode=CompartmentMode.ONE_COMPARTMENT,
        solver_method="Radau",
    )
    
    events = [
        DoseEvent(name="Warfarin", dose_mg=5.0, time_h=0.0),
        DoseEvent(name="Warfarin", dose_mg=5.0, time_h=24.0),
        DoseEvent(name="Warfarin", dose_mg=5.0, time_h=48.0),
        DoseEvent(name="Piperine", dose_mg=20.0, time_h=48.0),
        DoseEvent(name="Warfarin", dose_mg=5.0, time_h=72.0),
        DoseEvent(name="Piperine", dose_mg=20.0, time_h=72.0),
        DoseEvent(name="Warfarin", dose_mg=5.0, time_h=96.0),
        DoseEvent(name="Piperine", dose_mg=20.0, time_h=96.0),
        DoseEvent(name="Warfarin", dose_mg=5.0, time_h=120.0),
        DoseEvent(name="Warfarin", dose_mg=5.0, time_h=144.0),
    ]
    
    result = engine.simulate(
        events=events,
        duration_h=168.0,
        step_h=0.5,
    )
    
    assert result.success is True
    assert result.solver_status == "converged"
    
    # Check enzyme activity suppression during piperine phase (t=48 to t=96)
    times = np.array(result.time_points)
    enzyme_cyp3a4 = np.array(result.enzyme_activity_curves.get("CYP3A4", []))
    assert len(enzyme_cyp3a4) > 0
    
    # At t=0, enzyme is 1.0
    assert np.isclose(enzyme_cyp3a4[0], 1.0, atol=1e-3)
    # At t=96 (peak inhibition), enzyme should be depressed
    idx_96 = np.argmin(np.abs(times - 96.0))
    assert enzyme_cyp3a4[idx_96] < 0.95
    
    # At t=168 (washout), enzyme recovery is occurring
    idx_168 = np.argmin(np.abs(times - 168.0))
    assert enzyme_cyp3a4[idx_168] > enzyme_cyp3a4[idx_96]


def test_scenario_b_phenytoin_curcumin_nonlinear(catalog):
    """Scenario B: Phenytoin Michaelis-Menten Nonlinear Elimination Saturation."""
    drug = catalog.get_drug("Phenytoin")
    herb = catalog.get_herb("Curcumin")
    
    engine = TPSMEngine(
        drugs=[drug],
        herbs=[herb],
        compartment_mode=CompartmentMode.ONE_COMPARTMENT,
        solver_method="Radau",
    )
    
    events = [
        DoseEvent(name="Phenytoin", dose_mg=300.0, time_h=0.0),
        DoseEvent(name="Phenytoin", dose_mg=300.0, time_h=24.0),
        DoseEvent(name="Curcumin", dose_mg=500.0, time_h=24.0),
        DoseEvent(name="Phenytoin", dose_mg=300.0, time_h=48.0),
        DoseEvent(name="Curcumin", dose_mg=500.0, time_h=48.0),
        DoseEvent(name="Phenytoin", dose_mg=300.0, time_h=72.0),
        DoseEvent(name="Curcumin", dose_mg=500.0, time_h=72.0),
        DoseEvent(name="Phenytoin", dose_mg=300.0, time_h=96.0),
    ]
    
    result = engine.simulate(
        events=events,
        duration_h=120.0,
        step_h=0.5,
    )
    
    assert result.success is True
    conc = np.array(result.concentration_curves["Phenytoin"])
    times = np.array(result.time_points)
    
    idx_24 = np.argmin(np.abs(times - 24.0))
    idx_96 = np.argmin(np.abs(times - 96.0))
    
    # Trough at 96h should exceed trough at 24h due to accumulation
    assert conc[idx_96] > conc[idx_24]
    assert result.c_max["Phenytoin"] > 5.0


def test_scenario_c_midazolam_piperine_reversible(catalog):
    """Scenario C: Midazolam single dose with vs without Piperine."""
    drug = catalog.get_drug("Midazolam")
    herb = catalog.get_herb("Piperine")
    
    # Run combo
    engine_combo = TPSMEngine(
        drugs=[drug],
        herbs=[herb],
        solver_method="Radau",
    )
    events_combo = [
        DoseEvent(name="Midazolam", dose_mg=7.5, time_h=0.0),
        DoseEvent(name="Piperine", dose_mg=20.0, time_h=0.0),
    ]
    res_combo = engine_combo.simulate(
        events=events_combo,
        duration_h=24.0,
        step_h=0.25,
    )
    
    # Run alone
    engine_alone = TPSMEngine(
        drugs=[drug],
        herbs=[],
        solver_method="Radau",
    )
    events_alone = [
        DoseEvent(name="Midazolam", dose_mg=7.5, time_h=0.0),
    ]
    res_alone = engine_alone.simulate(
        events=events_alone,
        duration_h=24.0,
        step_h=0.25,
    )
    
    assert res_combo.success is True
    assert res_alone.success is True
    
    auc_combo = res_combo.auc["Midazolam"]
    auc_alone = res_alone.auc["Midazolam"]
    assert auc_combo > auc_alone


def test_scenario_d_paracetamol_repeated_accumulation(catalog):
    """Scenario D: Paracetamol Q6H multi-dose accumulation and complete washout."""
    drug = catalog.get_drug("Paracetamol")
    
    engine = TPSMEngine(
        drugs=[drug],
        herbs=[],
        solver_method="Radau",
    )
    
    events = [
        DoseEvent(name="Paracetamol", dose_mg=500.0, time_h=0.0),
        DoseEvent(name="Paracetamol", dose_mg=500.0, time_h=6.0),
        DoseEvent(name="Paracetamol", dose_mg=500.0, time_h=12.0),
        DoseEvent(name="Paracetamol", dose_mg=500.0, time_h=18.0),
    ]
    
    result = engine.simulate(
        events=events,
        duration_h=48.0,
        step_h=0.25,
    )
    
    assert result.success is True
    conc = np.array(result.concentration_curves["Paracetamol"])
    times = np.array(result.time_points)
    
    # Peak occurs after multiple doses (at or after dose 3 at t=12h)
    peak_time = times[np.argmax(conc)]
    assert peak_time >= 12.0
    
    # By 48h (after 30h washout), paracetamol is cleared
    assert conc[-1] < 0.2


def test_scenario_e_digoxin_renal_caution(catalog):
    """Scenario E: Digoxin accumulation and narrow therapeutic index."""
    drug = catalog.get_drug("Digoxin")
    
    engine = TPSMEngine(
        drugs=[drug],
        herbs=[],
        solver_method="Radau",
    )
    
    events = [
        DoseEvent(name="Digoxin", dose_mg=0.25, time_h=0.0),
        DoseEvent(name="Digoxin", dose_mg=0.25, time_h=24.0),
        DoseEvent(name="Digoxin", dose_mg=0.25, time_h=48.0),
        DoseEvent(name="Digoxin", dose_mg=0.25, time_h=72.0),
    ]
    
    result = engine.simulate(
        events=events,
        duration_h=96.0,
        step_h=0.5,
    )
    
    assert result.success is True
    conc = np.array(result.concentration_curves["Digoxin"])
    times = np.array(result.time_points)
    
    trough_24 = conc[np.argmin(np.abs(times - 24.0))]
    trough_72 = conc[np.argmin(np.abs(times - 72.0))]
    assert trough_72 > trough_24


def test_all_scenarios_from_json(catalog, scenarios_data):
    """Test loading and running every scenario directly from the JSON specification."""
    for scen in scenarios_data:
        drugs = [catalog.get_drug(d) for d in scen["drugs"]]
        herbs = [catalog.get_herb(h) for h in scen["herbs"]]
        
        events = [
            DoseEvent(name=e["name"], dose_mg=float(e["dose_mg"]), time_h=float(e["time_h"]))
            for e in scen["events"]
        ]
        
        comp_mode = (
            CompartmentMode.TWO_COMPARTMENT
            if scen.get("compartment_mode") == "2-compartment"
            else CompartmentMode.ONE_COMPARTMENT
        )
        
        engine = TPSMEngine(
            drugs=drugs,
            herbs=herbs,
            compartment_mode=comp_mode,
            solver_method=scen.get("solver_method", "Radau"),
        )
        
        result = engine.simulate(
            events=events,
            duration_h=float(scen["duration_h"]),
            step_h=float(scen["step_h"]),
        )
        
        assert result.success is True, f"Failed on {scen['scenario_id']}: {result.solver_message}"
        assert len(result.time_points) > 0
        for d in scen["drugs"]:
            assert result.auc[d] >= 0.0
