"""Deterministic temporal pharmacokinetic state machine.

This is an educational one-compartment simulator. It is not a PBPK model and
must not be used for dosing or patient-specific treatment decisions.
"""

from dataclasses import dataclass, field
from math import exp, isfinite
from typing import Dict, Iterable, List, Optional


@dataclass(frozen=True)
class TemporalDrug:
    name: str
    bioavailability: float
    absorption_rate_h: float
    volume_l: float
    clearance_l_h: float
    therapeutic_min_mg_l: Optional[float] = None
    toxicity_threshold_mg_l: Optional[float] = None
    half_life_h: Optional[float] = None
    age_restriction: str = "Requires clinical review"
    side_effects: tuple[str, ...] = ()


@dataclass(frozen=True)
class TemporalHerb:
    name: str
    bioavailability: float
    absorption_rate_h: float
    volume_l: float
    clearance_l_h: float
    enzyme: Optional[str] = None
    inhibition_strength: float = 0.0
    enzyme_recovery_half_life_h: Optional[float] = None
    onset_h: Optional[float] = None
    duration_h: Optional[float] = None
    age_restriction: str = "Requires clinical review"
    side_effects: tuple[str, ...] = ()


@dataclass(frozen=True)
class DoseEvent:
    name: str
    dose_mg: float
    time_h: float


@dataclass(frozen=True)
class TemporalPoint:
    time_h: float
    concentrations_mg_l: Dict[str, float]
    active_enzyme_fraction: Dict[str, float]
    toxicity_flags: Dict[str, bool]


@dataclass(frozen=True)
class TPSMResult:
    points: tuple[TemporalPoint, ...]
    assumptions: tuple[str, ...]
    warnings: tuple[str, ...]

    def series(self, name: str) -> list[tuple[float, float]]:
        key = name.lower()
        return [
            (point.time_h, point.concentrations_mg_l.get(key, 0.0))
            for point in self.points
        ]

    def first_toxicity_time(self, name: str) -> Optional[float]:
        target = name.lower()
        for point in self.points:
            if any(
                drug_name.lower() == target and flag
                for drug_name, flag in point.toxicity_flags.items()
            ):
                return point.time_h
        return None


class TemporalPKStateMachine:
    """Simulate approximate exposure and enzyme recovery over a timeline."""

    def __init__(
        self,
        drugs: Iterable[TemporalDrug],
        herbs: Iterable[TemporalHerb] = (),
        enzyme_baseline: Optional[Dict[str, float]] = None,
    ) -> None:
        self.drugs = {item.name.lower(): item for item in drugs}
        self.herbs = {item.name.lower(): item for item in herbs}
        self.enzyme_baseline = {
            key: value for key, value in (enzyme_baseline or {}).items()
        }
        for value in self.enzyme_baseline.values():
            if not 0.0 <= value <= 1.0:
                raise ValueError("enzyme baseline fractions must be between 0 and 1")

    def simulate(
        self,
        events: Iterable[DoseEvent],
        duration_h: float,
        step_h: float = 0.25,
    ) -> TPSMResult:
        if not isfinite(duration_h) or duration_h <= 0:
            raise ValueError("duration_h must be finite and positive")
        if not isfinite(step_h) or step_h <= 0 or step_h > 24:
            raise ValueError("step_h must be finite, positive, and at most 24 hours")

        event_list = sorted(events, key=lambda event: event.time_h)
        if any(event.time_h < 0 or event.time_h > duration_h for event in event_list):
            raise ValueError("dose event times must be within the simulation window")
        if any(event.dose_mg <= 0 or not isfinite(event.dose_mg) for event in event_list):
            raise ValueError("dose amounts must be finite and positive")

        known = {**self.drugs, **self.herbs}
        if any(event.name.lower() not in known for event in event_list):
            raise ValueError("all dose events must reference a registered drug or herb")

        amounts = {name: 0.0 for name in known}
        enzyme_state = dict(self.enzyme_baseline)
        for herb in self.herbs.values():
            if herb.enzyme and herb.enzyme not in enzyme_state:
                enzyme_state[herb.enzyme] = 1.0
        points: List[TemporalPoint] = []
        assumptions = (
            "One-compartment oral model with first-order absorption and elimination.",
            "Herb parameters and enzyme inhibition strengths are illustrative until sourced.",
            "This simulator is not validated for dosing, diagnosis, or individual treatment.",
        )
        warnings = tuple(
            sorted(
                {
                    f"{item.name}: {item.age_restriction}"
                    for item in (*self.drugs.values(), *self.herbs.values())
                    if item.age_restriction
                }
            )
        )
        event_index = 0
        time_h = 0.0
        while time_h <= duration_h + 1e-9:
            while event_index < len(event_list) and event_list[event_index].time_h <= time_h + 1e-9:
                event = event_list[event_index]
                item = known[event.name.lower()]
                amounts[item.name.lower()] += event.dose_mg * item.bioavailability
                event_index += 1

            concentrations = {
                name: amount / known[name].volume_l for name, amount in amounts.items()
            }
            points.append(
                TemporalPoint(
                    time_h=round(time_h, 6),
                    concentrations_mg_l=dict(concentrations),
                    active_enzyme_fraction=dict(enzyme_state),
                    toxicity_flags={
                        item.name: (
                            item.toxicity_threshold_mg_l is not None
                            and concentrations[item.name.lower()] >= item.toxicity_threshold_mg_l
                        )
                        for item in self.drugs.values()
                    },
                )
            )
            if time_h >= duration_h:
                break

            dt = min(step_h, duration_h - time_h)
            herb_inhibition: Dict[str, float] = {}
            for herb in self.herbs.values():
                if herb.enzyme and herb_inhibition.get(herb.enzyme, 0.0) is not None:
                    concentration = concentrations[herb.name.lower()]
                    herb_inhibition[herb.enzyme] = max(
                        herb_inhibition.get(herb.enzyme, 0.0),
                        herb.inhibition_strength * concentration,
                    )

            for herb in self.herbs.values():
                if herb.enzyme:
                    recovery = 1.0
                    if herb.enzyme_recovery_half_life_h:
                        recovery = 1.0 - 2 ** (-dt / herb.enzyme_recovery_half_life_h)
                    inhibition = min(1.0, herb_inhibition.get(herb.enzyme, 0.0))
                    enzyme_state[herb.enzyme] = min(
                        1.0,
                        max(0.0, enzyme_state.get(herb.enzyme, 1.0) + recovery * (1.0 - enzyme_state.get(herb.enzyme, 1.0)) - inhibition * dt),
                    )

            for name, item in known.items():
                rate = item.clearance_l_h / item.volume_l
                if isinstance(item, TemporalDrug) and item.toxicity_threshold_mg_l is not None:
                    enzyme_factor = min(enzyme_state.values()) if enzyme_state else 1.0
                    rate *= max(0.05, enzyme_factor)
                amounts[name] *= exp(-rate * dt)
            time_h += dt

        return TPSMResult(tuple(points), assumptions, warnings)
