"""High-level simulation orchestrator for the Temporal Pharmacokinetic State Machine."""

from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple
import numpy as np

from .equations import PKPDEquations, StateLayout
from .events import DoseEvent, EventTimelineProcessor
from .metrics import PKMetricsCalculator
from .models import (
    CompartmentMode,
    DrugParameters,
    HerbParameters,
    SimulationResult,
    TemporalPoint,
    TPSMSolverError,
    TPSMValidationError,
)
from .solver import ODESolverWrapper
from .validation import SimulationValidator


class TPSMEngine:
    """Orchestrates exact event segmentation, stiff numerical integration, and metric derivation."""

    def __init__(
        self,
        drugs: Iterable[DrugParameters],
        herbs: Iterable[HerbParameters] = (),
        compartment_mode: CompartmentMode = CompartmentMode.ONE_COMPARTMENT,
        solver_method: str = "Radau",
        enzyme_baseline: Optional[Dict[str, float]] = None,
        enzyme_recovery_halflives: Optional[Dict[str, float]] = None,
    ):
        self.drugs = list(drugs)
        self.herbs = list(herbs)
        self.compartment_mode = compartment_mode
        self.solver_method = solver_method
        self.enzyme_baseline = enzyme_baseline or {}
        self.enzyme_kdeg: Dict[str, float] = {}

        # Validate all drug and herb parameters upfront
        for d in self.drugs:
            SimulationValidator.validate_drug(d, compartment_mode)
        for h in self.herbs:
            SimulationValidator.validate_herb(h)

        # Collect all active enzymes/pathways
        all_enzymes: Set[str] = set()
        for d in self.drugs:
            all_enzymes.update(d.fm_pathways.keys())
        for h in self.herbs:
            all_enzymes.update(h.pathway_effects.keys())

        self.enzymes = sorted(all_enzymes)

        # Set up enzyme turnover k_deg (1/h)
        for enz in self.enzymes:
            hl = (enzyme_recovery_halflives or {}).get(enz, 36.0)
            self.enzyme_kdeg[enz] = np.log(2.0) / max(1.0, hl)

        self.layout = StateLayout(
            drugs=self.drugs,
            herbs=self.herbs,
            enzymes=self.enzymes,
            compartment_mode=compartment_mode,
        )

        self.equations = PKPDEquations(
            layout=self.layout,
            enzyme_kdeg=self.enzyme_kdeg,
        )

        self.solver = ODESolverWrapper(method=solver_method)

    def simulate(
        self,
        events: Iterable[DoseEvent],
        duration_h: float,
        step_h: float = 0.25,
    ) -> SimulationResult:
        """
        Execute deterministic timeline simulation across discrete dose events.
        """
        # 1. Validate events
        known_names = set(d.name.lower() for d in self.drugs) | set(h.name.lower() for h in self.herbs)
        sorted_events = EventTimelineProcessor.validate_events(events, duration_h, known_names)

        # 2. Segment timeline
        segments = EventTimelineProcessor.segment_timeline(sorted_events, duration_h)

        # 3. Setup initial state
        y_current = self.layout.initial_state(self.enzyme_baseline)

        # Tracking arrays
        all_times: List[float] = []
        all_states: List[np.ndarray] = []
        toxicity_events_collected: Dict[str, List[float]] = {d.name: [] for d in self.drugs}
        event_log: List[dict] = []

        # Name mapping
        drug_map = {d.name.lower(): (i, d) for i, d in enumerate(self.drugs)}
        herb_map = {h.name.lower(): (m, h) for m, h in enumerate(self.herbs)}

        # Build toxicity root-finding events
        root_events = self.equations.make_toxicity_events()

        # If there are no segments (e.g. duration_h == 0 is rejected by validation), handle [0, duration_h]
        if not segments:
            segments = [EventTimelineProcessor.segment_timeline([], duration_h)[0]]

        last_successful_time = 0.0

        for seg_idx, seg in enumerate(segments):
            # Apply impulses at segment start
            for dose_ev in seg.doses_at_start:
                name_key = dose_ev.name.lower()
                if name_key in drug_map:
                    idx, drug = drug_map[name_key]
                    absorbed_mass = dose_ev.dose_mg * drug.bioavailability
                    y_current[self.layout.idx_v_ag + idx] += absorbed_mass
                elif name_key in herb_map:
                    idx, herb = herb_map[name_key]
                    absorbed_mass = dose_ev.dose_mg * herb.bioavailability
                    y_current[self.layout.idx_h_ag + idx] += absorbed_mass

                event_log.append({
                    "time_h": seg.t_start,
                    "compound": dose_ev.name,
                    "dose_mg": dose_ev.dose_mg,
                    "bioavailability": drug.bioavailability if name_key in drug_map else herb.bioavailability,
                })

            # Check state before integration
            if not np.all(np.isfinite(y_current)):
                raise TPSMSolverError(
                    f"Non-finite state detected at t={seg.t_start}h prior to ODE integration step"
                )

            # Determine number of output evaluation points for this segment
            duration_seg = seg.t_end - seg.t_start
            n_eval = max(2, int(round(duration_seg / max(0.0001, step_h))) + 1)

            # Integrate segment
            try:
                t_eval, y_eval, t_ev_list, success, msg = self.solver.integrate_segment(
                    fun=self.equations.rhs,
                    t_span=(seg.t_start, seg.t_end),
                    y0=y_current,
                    events=root_events,
                    num_points=n_eval,
                )
            except Exception as e:
                # Do not advance simulation clock on failure!
                raise TPSMSolverError(
                    f"TPSM solver failure at t={last_successful_time:.2f}h on interval [{seg.t_start:.2f}, {seg.t_end:.2f}]h: {e}"
                ) from e

            # Collect toxicity crossing times
            ev_idx = 0
            for k, drug in enumerate(self.drugs):
                if drug.toxicity_threshold_mg_l is not None:
                    if ev_idx < len(t_ev_list) and len(t_ev_list[ev_idx]) > 0:
                        for ev_time in t_ev_list[ev_idx]:
                            toxicity_events_collected[drug.name].append(float(ev_time))
                    ev_idx += 1

            # Append trajectory (avoid duplicating segment boundary point)
            if all_times and len(t_eval) > 0 and abs(t_eval[0] - all_times[-1]) < 1e-9:
                t_eval = t_eval[1:]
                y_eval = y_eval[:, 1:]

            for t_pt, y_pt in zip(t_eval, y_eval.T):
                all_times.append(float(t_pt))
                all_states.append(y_pt.copy())

            # Advance current state to end of segment
            y_current = y_eval[:, -1].copy()
            last_successful_time = seg.t_end

        # If zero events were provided or t_eval was empty, ensure at least [0, duration] points
        if not all_times:
            all_times = [0.0, duration_h]
            all_states = [y_current, y_current]

        time_array = np.array(all_times, dtype=np.float64)
        state_matrix = np.array(all_states, dtype=np.float64)  # Shape: [N_points, Total_vars]

        # 4. Extract Concentration and Enzyme Curves
        conc_curves: Dict[str, List[float]] = {}
        for k, drug in enumerate(self.drugs):
            amounts = state_matrix[:, self.layout.idx_v_ac + k]
            concs = [float(max(0.0, a / drug.volume_central_l)) for a in amounts]
            conc_curves[drug.name] = concs

        for m, herb in enumerate(self.herbs):
            amounts = state_matrix[:, self.layout.idx_h_ac + m]
            concs = [float(max(0.0, a / herb.volume_central_l)) for a in amounts]
            conc_curves[herb.name] = concs

        enzyme_curves: Dict[str, List[float]] = {}
        for j, enz in enumerate(self.enzymes):
            e_vals = state_matrix[:, self.layout.idx_e + j]
            enzyme_curves[enz] = [float(np.clip(val, 0.0, 1.0)) for val in e_vals]

        # 5. Compute PK Metrics
        auc_results: Dict[str, float] = {}
        c_max_results: Dict[str, float] = {}
        c_min_results: Dict[str, float] = {}
        tmr_curves: Dict[str, List[Optional[float]]] = {}
        first_breach_times: Dict[str, Optional[float]] = {}

        for k, drug in enumerate(self.drugs):
            c_array = np.array(conc_curves[drug.name], dtype=np.float64)
            auc_results[drug.name] = PKMetricsCalculator.calculate_auc(time_array, c_array)
            cmax, cmin = PKMetricsCalculator.calculate_cmax_cmin(c_array)
            c_max_results[drug.name] = cmax
            c_min_results[drug.name] = cmin
            tmr_curves[drug.name] = PKMetricsCalculator.calculate_tmr_curve(c_array, drug)

            # Earliest breach from solver root-finding or grid
            t_events_for_drug = toxicity_events_collected.get(drug.name)
            breach_t = PKMetricsCalculator.detect_first_toxicity_breach(
                time_array, c_array, drug.toxicity_threshold_mg_l,
                [np.array(t_events_for_drug)] if t_events_for_drug else None
            )
            first_breach_times[drug.name] = breach_t

        for m, herb in enumerate(self.herbs):
            c_array = np.array(conc_curves[herb.name], dtype=np.float64)
            auc_results[herb.name] = PKMetricsCalculator.calculate_auc(time_array, c_array)
            cmax, cmin = PKMetricsCalculator.calculate_cmax_cmin(c_array)
            c_max_results[herb.name] = cmax
            c_min_results[herb.name] = cmin

        # 6. Provenance, Assumptions, Warnings
        provenance_dict: Dict[str, dict] = {}
        warnings: List[str] = []
        for d in self.drugs:
            if d.provenance:
                provenance_dict[d.name] = {
                    "source": d.provenance.source_reference,
                    "evidence_level": d.provenance.evidence_level.value,
                    "status": d.provenance.parameter_status.value,
                    "uncertainty": {
                        "low": d.provenance.uncertainty_range.low,
                        "high": d.provenance.uncertainty_range.high,
                    } if d.provenance.uncertainty_range else None,
                    "population": d.provenance.population_assumptions,
                    "limitations": d.provenance.limitations,
                }
            if d.age_restrictions:
                warnings.append(f"{d.name}: {d.age_restrictions}")
            if d.renal_hepatic_cautions:
                warnings.append(f"{d.name}: {d.renal_hepatic_cautions}")

        for h in self.herbs:
            if h.provenance:
                provenance_dict[h.name] = {
                    "source": h.provenance.source_reference,
                    "evidence_level": h.provenance.evidence_level.value,
                    "status": h.provenance.parameter_status.value,
                    "preparation": h.provenance.preparation,
                    "limitations": h.provenance.limitations,
                }
            if h.age_restrictions:
                warnings.append(f"{h.name}: {h.age_restrictions}")

        assumptions = [
            f"{self.compartment_mode.value.capitalize()} open disposition model with first-order absorption.",
            "Dynamic metabolic clearance scaling with Mechanism-Based Inactivation (MBI) and de novo enzyme synthesis.",
            "Unbound perpetrator concentration drives target enzyme inhibition kinetics.",
            "All parameters reflect healthy reference populations unless renal/hepatic adjustments are explicitly applied.",
            "Research prototype simulation: does not replace professional clinical therapeutic drug monitoring.",
        ]

        uncertainty_metadata = {
            "model_type": self.compartment_mode.value,
            "solver": self.solver_method,
            "discretization_step_h": step_h,
            "enzymes_modeled": self.enzymes,
            "enzyme_turnover_kdeg": {k: round(v, 6) for k, v in self.enzyme_kdeg.items()},
        }

        # Final state dictionary
        final_state: Dict[str, float] = {}
        for k, drug in enumerate(self.drugs):
            final_state[f"{drug.name}_conc_mg_l"] = conc_curves[drug.name][-1]
        for m, herb in enumerate(self.herbs):
            final_state[f"{herb.name}_conc_mg_l"] = conc_curves[herb.name][-1]
        for enz in self.enzymes:
            final_state[f"{enz}_activity_fraction"] = enzyme_curves[enz][-1]

        return SimulationResult(
            solver_success=True,
            solver_method=self.solver_method,
            solver_message="Exact event-segmented integration completed successfully.",
            time_points=[round(t, 4) for t in all_times],
            concentration_curves={k: [round(v, 6) for v in vals] for k, vals in conc_curves.items()},
            enzyme_activity_curves={k: [round(v, 6) for v in vals] for k, vals in enzyme_curves.items()},
            auc={k: round(v, 4) for k, v in auc_results.items()},
            c_max={k: round(v, 4) for k, v in c_max_results.items()},
            c_min={k: round(v, 4) for k, v in c_min_results.items()},
            first_toxicity_breach_time=first_breach_times,
            therapeutic_margin_ratio_curves=tmr_curves,
            final_state=final_state,
            event_log=event_log,
            parameter_provenance=provenance_dict,
            warnings=sorted(set(warnings)),
            assumptions=assumptions,
            uncertainty_metadata=uncertainty_metadata,
        )


class TemporalPKStateMachine:
    """
    Backward-compatible facade maintaining the public interface of previous prototype
    while delegating to the research-grade TPSMEngine.
    """

    def __init__(
        self,
        drugs: Iterable[Any],
        herbs: Iterable[Any] = (),
        enzyme_baseline: Optional[Dict[str, float]] = None,
        compartment_mode: CompartmentMode = CompartmentMode.ONE_COMPARTMENT,
        solver_method: str = "Radau",
    ):
        converted_drugs = [self._convert_drug(d) for d in drugs]
        converted_herbs = [self._convert_herb(h) for h in herbs]
        self._engine = TPSMEngine(
            drugs=converted_drugs,
            herbs=converted_herbs,
            compartment_mode=compartment_mode,
            solver_method=solver_method,
            enzyme_baseline=enzyme_baseline,
        )
        self.drugs = {d.name.lower(): d for d in converted_drugs}
        self.herbs = {h.name.lower(): h for h in converted_herbs}

    def simulate(
        self,
        events: Iterable[DoseEvent],
        duration_h: float,
        step_h: float = 0.25,
    ) -> SimulationResult:
        return self._engine.simulate(events, duration_h, step_h)

    @staticmethod
    def _convert_drug(d: Any) -> DrugParameters:
        if isinstance(d, DrugParameters):
            return d
        # Convert from legacy TemporalDrug if needed
        return DrugParameters(
            name=getattr(d, "name"),
            bioavailability=getattr(d, "bioavailability", 1.0),
            ka_h=getattr(d, "absorption_rate_h", getattr(d, "ka_h", 1.0)),
            volume_central_l=getattr(d, "volume_l", getattr(d, "volume_central_l", 10.0)),
            clearance_l_h=getattr(d, "clearance_l_h", 1.0),
            therapeutic_min_mg_l=getattr(d, "therapeutic_min_mg_l", None),
            toxicity_threshold_mg_l=getattr(d, "toxicity_threshold_mg_l", None),
            half_life_h=getattr(d, "half_life_h", None),
            age_restrictions=getattr(d, "age_restriction", getattr(d, "age_restrictions", "")),
            side_effects=tuple(getattr(d, "side_effects", ())),
        )

    @staticmethod
    def _convert_herb(h: Any) -> HerbParameters:
        if isinstance(h, HerbParameters):
            return h
        # Convert from legacy TemporalHerb if needed
        enz = getattr(h, "enzyme", None)
        pathways = {}
        if enz:
            from .models import EnzymePathwayParam, PathwayType
            pathways[enz] = EnzymePathwayParam(
                pathway_name=enz,
                pathway_type=PathwayType.CYP_ENZYME,
                kinact_h=getattr(h, "inhibition_strength", 0.5),
                ki_u_mbi_mg_l=0.5,
                recovery_half_life_h=getattr(h, "enzyme_recovery_half_life_h", 36.0),
            )
        return HerbParameters(
            name=getattr(h, "name"),
            bioavailability=getattr(h, "bioavailability", 0.5),
            ka_h=getattr(h, "absorption_rate_h", getattr(h, "ka_h", 1.0)),
            volume_central_l=getattr(h, "volume_l", getattr(h, "volume_central_l", 50.0)),
            clearance_l_h=getattr(h, "clearance_l_h", 10.0),
            pathway_effects=pathways,
            onset_h=getattr(h, "onset_h", None),
            duration_h=getattr(h, "duration_h", None),
            age_restrictions=getattr(h, "age_restriction", getattr(h, "age_restrictions", "")),
            side_effects=tuple(getattr(h, "side_effects", ())),
        )
