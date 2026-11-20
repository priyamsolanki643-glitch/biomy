"""Unit tests for PK/PD differential equations, rate calculations, and vector field evaluation."""

import numpy as np
import pytest
from biomy.tpsm.equations import PKPDEquations, StateLayout
from biomy.tpsm.models import (
    CompartmentMode,
    DrugParameters,
    EliminationMode,
    EnzymePathwayParam,
    HerbParameters,
    PathwayType,
)


@pytest.fixture
def sample_system():
    drug = DrugParameters(
        name="TestDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=1.0,
        fm_pathways={"CYP3A4": 0.8},
        toxicity_threshold_mg_l=5.0,
    )
    herb = HerbParameters(
        name="TestHerb",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=20.0,
        clearance_l_h=2.0,
        fraction_unbound=1.0,
        pathway_effects={
            "CYP3A4": EnzymePathwayParam(
                pathway_name="CYP3A4",
                pathway_type=PathwayType.CYP_ENZYME,
                kinact_h=1.5,
                ki_u_mbi_mg_l=0.5,
                recovery_half_life_h=30.0,
            )
        },
    )
    layout = StateLayout([drug], [herb], ["CYP3A4"], CompartmentMode.ONE_COMPARTMENT)
    eq = PKPDEquations(layout)
    return drug, herb, layout, eq


def test_enzyme_starts_at_100_percent_baseline(sample_system):
    _, _, layout, _ = sample_system
    y0 = layout.initial_state()
    assert y0[layout.idx_e] == 1.0


def test_inhibitor_lowers_enzyme_activity(sample_system):
    _, _, layout, eq = sample_system
    y = layout.initial_state()
    # Place perpetrator in central compartment (high concentration)
    y[layout.idx_h_ac] = 50.0  # 50 mg in 20 L -> 2.5 mg/L

    dy = eq.rhs(0.0, y)
    # dE_rel/dt should be negative due to active MBI
    d_enz = dy[layout.idx_e]
    assert d_enz < 0.0


def test_linear_and_nonlinear_elimination_produce_different_trajectories():
    drug_linear = DrugParameters(
        name="LinearDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=2.0,
        elimination_mode=EliminationMode.LINEAR,
    )
    drug_nonlinear = DrugParameters(
        name="NonlinearDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=2.0,
        elimination_mode=EliminationMode.MICHAELIS_MENTEN,
        vmax_mg_h=20.0,
        km_mg_l=2.0,
    )
    layout_lin = StateLayout([drug_linear], [], [], CompartmentMode.ONE_COMPARTMENT)
    layout_nonlin = StateLayout([drug_nonlinear], [], [], CompartmentMode.ONE_COMPARTMENT)

    eq_lin = PKPDEquations(layout_lin)
    eq_nonlin = PKPDEquations(layout_nonlin)

    # State with drug in central compartment: 100 mg -> 10 mg/L (well above Km)
    y_lin = layout_lin.initial_state()
    y_lin[layout_lin.idx_v_ac] = 100.0

    y_nonlin = layout_nonlin.initial_state()
    y_nonlin[layout_nonlin.idx_v_ac] = 100.0

    dy_lin = eq_lin.rhs(0.0, y_lin)
    dy_nonlin = eq_nonlin.rhs(0.0, y_nonlin)

    elim_rate_lin = -dy_lin[layout_lin.idx_v_ac]
    elim_rate_nonlin = -dy_nonlin[layout_nonlin.idx_v_ac]

    # Linear: CL * C = 2.0 * 10 = 20 mg/h
    # Nonlinear: Vmax * C / (Km + C) = 20 * 10 / (2 + 10) = 16.67 mg/h
    assert abs(elim_rate_lin - elim_rate_nonlin) > 1.0


def test_michaelis_menten_clearance_saturates_near_km():
    drug_mm = DrugParameters(
        name="MMDrug",
        bioavailability=1.0,
        ka_h=1.0,
        volume_central_l=10.0,
        clearance_l_h=2.0,
        elimination_mode=EliminationMode.MICHAELIS_MENTEN,
        vmax_mg_h=20.0,
        km_mg_l=4.0,
    )
    layout = StateLayout([drug_mm], [], [], CompartmentMode.ONE_COMPARTMENT)
    eq = PKPDEquations(layout)

    # C = 40 mg/L (10x Km) -> elimination approaches Vmax (20 mg/h)
    y = layout.initial_state()
    y[layout.idx_v_ac] = 400.0  # 400 mg in 10 L = 40 mg/L

    dy = eq.rhs(0.0, y)
    elim_rate = -dy[layout.idx_v_ac]
    # Vmax * 40 / 44 = 18.18 mg/h (close to 20 mg/h saturation)
    assert 17.5 < elim_rate <= 20.0
