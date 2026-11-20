"""Dose event processing, validation, and timeline segmentation."""

from dataclasses import dataclass
from math import isfinite
from typing import Dict, Iterable, List, Sequence, Set

from .models import DoseEvent, TPSMValidationError


@dataclass(frozen=True)
class EventSegment:
    """Simulation window bounded by discrete event impulses at t_start."""
    t_start: float
    t_end: float
    doses_at_start: List[DoseEvent]


class EventTimelineProcessor:
    """Validates, sorts, deduplicates, and segments dose events for exact numerical integration."""

    @staticmethod
    def validate_events(
        events: Iterable[DoseEvent],
        simulation_duration_h: float,
        known_compounds: Set[str],
    ) -> List[DoseEvent]:
        """Validate input events against simulation horizon and known compound catalog."""
        if not isfinite(simulation_duration_h) or simulation_duration_h <= 0:
            raise TPSMValidationError(f"duration_h must be finite and positive: got {simulation_duration_h}")

        validated: List[DoseEvent] = []
        for i, ev in enumerate(events):
            if not isinstance(ev, DoseEvent):
                raise TPSMValidationError(f"Item {i} is not a DoseEvent: got {type(ev)}")

            if not isfinite(ev.dose_mg) or ev.dose_mg <= 0:
                raise TPSMValidationError(f"Dose amount must be finite and strictly positive: got {ev.dose_mg} for {ev.name}")

            if not isfinite(ev.time_h) or ev.time_h < 0:
                raise TPSMValidationError(f"Dose timestamp must be finite and non-negative: got {ev.time_h} for {ev.name}")

            if ev.time_h > simulation_duration_h:
                raise TPSMValidationError(
                    f"Dose event at t={ev.time_h}h for {ev.name} exceeds simulation duration {simulation_duration_h}h"
                )

            name_lower = ev.name.strip().lower()
            if known_compounds and name_lower not in known_compounds:
                raise TPSMValidationError(f"Unknown compound '{ev.name}' in dose events; not present in catalog")

            validated.append(
                DoseEvent(
                    name=name_lower,
                    dose_mg=float(ev.dose_mg),
                    time_h=float(ev.time_h),
                    route=ev.route,
                    notes=ev.notes,
                )
            )

        # Sort strictly by timestamp
        return sorted(validated, key=lambda ev: ev.time_h)

    @staticmethod
    def segment_timeline(
        sorted_events: Sequence[DoseEvent],
        duration_h: float,
    ) -> List[EventSegment]:
        """
        Segment the total timeline [0, duration_h] by discrete event timestamps.
        Each segment integrates over [t_start, t_end] with impulses applied at t_start.
        """
        # Collect unique timestamps within [0, duration_h]
        event_times = sorted(set(ev.time_h for ev in sorted_events))
        if not event_times or event_times[0] > 0.0:
            timestamps = [0.0] + [t for t in event_times if t > 0.0]
        else:
            timestamps = list(event_times)

        if timestamps[-1] < duration_h:
            timestamps.append(duration_h)

        # Group doses by their exact timestamp
        doses_by_time: Dict[float, List[DoseEvent]] = {}
        for ev in sorted_events:
            doses_by_time.setdefault(ev.time_h, []).append(ev)

        segments: List[EventSegment] = []
        for i in range(len(timestamps) - 1):
            t_start = timestamps[i]
            t_end = timestamps[i + 1]
            if t_end <= t_start:
                continue
            doses = doses_by_time.get(t_start, [])
            segments.append(EventSegment(t_start=t_start, t_end=t_end, doses_at_start=doses))

        # If an event exists exactly at duration_h
        if duration_h in doses_by_time:
            segments.append(
                EventSegment(
                    t_start=duration_h,
                    t_end=duration_h,
                    doses_at_start=doses_by_time[duration_h],
                )
            )

        return segments
