"""Coupled pharmacokinetic and pharmacodynamic differential equations.

Implements:
1. One- and two-compartment open oral disposition models.
2. Linear and Michaelis-Menten nonlinear saturable elimination.
3. Reversible inhibition and irreversible Mechanism-Based Inactivation (MBI).
4. De novo enzyme synthesis and post-washout recovery dynamics.
5. Exact toxicity event functions for ODE solver root-finding.
"""

from math import log
from typing import Callable, Dict, List, Optional, Tuple
import numpy as np

from .models import CompartmentMode, DrugParameters, EliminationMode, HerbParameters


class StateLayout:
    """Manages index mapping for the continuous ODE state vector."""

    def __init__(
        self,
        drugs: List[DrugParameters],
        herbs: List[HerbParameters],
        enzymes: List[str],
        compartment_mode: CompartmentMode = CompartmentMode.ONE_COMPARTMENT,
    ):
        self.drugs = drugs
        self.herbs = herbs
        self.enzymes = sorted(enzymes)
        self.compartment_mode = compartment_mode

        self.n_d = len(drugs)
        self.n_h = len(herbs)
        self.n_e = len(self.enzymes)

        # Offsets
        self.idx_v_ag = 0
        self.idx_v_ac = self.n_d
        self.idx_v_ap = 2 * self.n_d if compartment_mode == CompartmentMode.TWO_COMPARTMENT else None

        v_slots = 3 * self.n_d if compartment_mode == CompartmentMode.TWO_COMPARTMENT else 2 * self.n_d

        self.idx_h_ag = v_slots
        self.idx_h_ac = v_slots + self.n_h
        self.idx_e = v_slots + 2 * self.n_h

        self.total_size = v_slots + 2 * self.n_h + self.n_e

    def initial_state(self, baseline_enzymes: Optional[Dict[str, float]] = None) -> np.ndarray:
        """Create baseline zero-amount state with 100% active enzymes."""
        y0 = np.zeros(self.total_size, dtype=np.float64)
        # Initialize enzymes to 1.0 (or specified baseline fraction)
        for i, enz in enumerate(self.enzymes):
            val = (baseline_enzymes or {}).get(enz, 1.0)
            y0[self.idx_e + i] = np.clip(val, 0.0, 1.0)
        return y0

    def get_victim_central_amount(self, y: np.ndarray, drug_idx: int) -> float:
        return float(y[self.idx_v_ac + drug_idx])

    def get_victim_concentration(self, y: np.ndarray, drug_idx: int) -> float:
        amount = self.get_victim_central_amount(y, drug_idx)
        vc = self.drugs[drug_idx].volume_central_l
        return max(0.0, amount / vc)


class PKPDEquations:
    """Evaluates the vector field dy/dt for coupled multi-drug, multi-herb, multi-enzyme PK/PD."""

    def __init__(
        self,
        layout: StateLayout,
        enzyme_kdeg: Optional[Dict[str, float]] = None,
        enzyme_baseline: Optional[Dict[str, float]] = None,
    ):
        self.layout = layout
        self.drugs = layout.drugs
        self.herbs = layout.herbs
        self.enzymes = layout.enzymes
        self.mode = layout.compartment_mode
        self.enzyme_baseline = enzyme_baseline or {}

        # Enzyme degradation rate constants k_deg (1/h)
        # Default turnover half-life ~ 36h -> k_deg = ln(2)/36 ~ 0.01925 1/h
        self.k_deg = np.zeros(layout.n_e, dtype=np.float64)
        self.e_base = np.ones(layout.n_e, dtype=np.float64)
        for i, enz in enumerate(self.enzymes):
            if enzyme_kdeg and enz in enzyme_kdeg:
                self.k_deg[i] = enzyme_kdeg[enz]
            else:
                self.k_deg[i] = log(2.0) / 36.0
            if enzyme_baseline and enz in enzyme_baseline:
                self.e_base[i] = float(np.clip(enzyme_baseline[enz], 0.0, 1.0))
            else:
                self.e_base[i] = 1.0

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        """Compute system derivative dy/dt."""
        dy = np.zeros_like(y)
        layout = self.layout

        # Extract current state segments
        v_ag = np.maximum(y[layout.idx_v_ag : layout.idx_v_ac], 0.0)
        v_ac = np.maximum(y[layout.idx_v_ac : layout.idx_v_ac + layout.n_d], 0.0)
        if layout.compartment_mode == CompartmentMode.TWO_COMPARTMENT:
            v_ap = np.maximum(y[layout.idx_v_ap : layout.idx_v_ap + layout.n_d], 0.0)

        h_ag = np.maximum(y[layout.idx_h_ag : layout.idx_h_ac], 0.0)
        h_ac = np.maximum(y[layout.idx_h_ac : layout.idx_h_ac + layout.n_h], 0.0)
        e_rel = np.clip(y[layout.idx_e : layout.idx_e + layout.n_e], 0.0, 1.0)

        # 1. Perpetrator (Herb) PK & Unbound Concentration
        h_unbound = np.zeros(layout.n_h, dtype=np.float64)
        for m, herb in enumerate(self.herbs):
            ka = herb.ka_h
            cl = herb.clearance_l_h
            vc = herb.volume_central_l
            dy[layout.idx_h_ag + m] = -ka * h_ag[m]
            dy[layout.idx_h_ac + m] = ka * h_ag[m] - (cl / vc) * h_ac[m]
            h_unbound[m] = herb.fraction_unbound * (h_ac[m] / vc)

        # 2. Dynamic Enzyme Pool (MBI + Synthesis/Degradation)
        # dE_rel/dt = k_deg * (1 - E_rel) - sum_m(k_obs_m) * E_rel
        rev_denom = {enz: 1.0 for enz in self.enzymes}

        for j, enz in enumerate(self.enzymes):
            k_obs_sum = 0.0
            for m, herb in enumerate(self.herbs):
                pathway = herb.pathway_effects.get(enz)
                if not pathway:
                    continue
                iu = h_unbound[m]
                # Reversible inhibition denominator: 1 + Iu / Ki_rev
                if pathway.reversible_ki_u_mg_l and pathway.reversible_ki_u_mg_l > 0:
                    rev_denom[enz] += iu / pathway.reversible_ki_u_mg_l
                # MBI inactivation rate: kinact * Iu / (KIu + Iu)
                if pathway.kinact_h and pathway.ki_u_mbi_mg_l and pathway.ki_u_mbi_mg_l > 0:
                    k_obs_sum += (pathway.kinact_h * iu) / (pathway.ki_u_mbi_mg_l + iu)

            kdeg = self.k_deg[j]
            ebase = self.e_base[j]
            dy[layout.idx_e + j] = kdeg * (ebase - e_rel[j]) - k_obs_sum * e_rel[j]

        # Map enzyme name to relative activity
        enz_map = {enz: e_rel[j] for j, enz in enumerate(self.enzymes)}

        # 3. Victim Drug PK with Dynamic Clearance Feedback
        for k, drug in enumerate(self.drugs):
            ka = drug.ka_h
            vc = drug.volume_central_l
            conc = v_ac[k] / vc

            # Clearance scaling: non-metabolized fraction + sum(fm_j * E_rel_j / rev_denom_j)
            fm_sum = sum(drug.fm_pathways.values())
            non_metabolized = max(0.0, 1.0 - fm_sum)
            cl_scale = non_metabolized
            for enz, fm_val in drug.fm_pathways.items():
                e_val = enz_map.get(enz, 1.0)
                denom = rev_denom.get(enz, 1.0)
                cl_scale += fm_val * (e_val / denom)

            cl_scale = max(0.01, cl_scale)

            # Elimination rate (mg/h)
            if drug.elimination_mode == EliminationMode.MICHAELIS_MENTEN and drug.vmax_mg_h and drug.km_mg_l:
                vmax_eff = drug.vmax_mg_h * cl_scale
                elim_rate = (vmax_eff * conc) / (drug.km_mg_l + conc)
            else:
                cl_eff = drug.clearance_l_h * cl_scale
                elim_rate = (cl_eff / vc) * v_ac[k]

            # Gut differential
            dy[layout.idx_v_ag + k] = -ka * v_ag[k]

            # Disposition differential
            if layout.compartment_mode == CompartmentMode.TWO_COMPARTMENT:
                vp = drug.volume_peripheral_l or (vc * 0.5)
                q = drug.inter_compartment_clearance_l_h or (drug.clearance_l_h * 0.5)
                k12 = q / vc
                k21 = q / vp
                dy[layout.idx_v_ac + k] = ka * v_ag[k] - elim_rate - k12 * v_ac[k] + k21 * v_ap[k]
                dy[layout.idx_v_ap + k] = k12 * v_ac[k] - k21 * v_ap[k]
            else:
                dy[layout.idx_v_ac + k] = ka * v_ag[k] - elim_rate

        return dy

    def make_toxicity_events(self) -> List[Callable[[float, np.ndarray], float]]:
        """Construct SciPy root-finding events for exact toxicity threshold crossings."""
        events = []
        for k, drug in enumerate(self.drugs):
            if drug.toxicity_threshold_mg_l is None:
                continue

            threshold = float(drug.toxicity_threshold_mg_l)
            drug_idx = k

            def event_func(t: float, y: np.ndarray, d_idx=drug_idx, thresh=threshold) -> float:
                conc = self.layout.get_victim_concentration(y, d_idx)
                return conc - thresh

            event_func.direction = 1.0  # Upward threshold breach
            events.append(event_func)

        return events
