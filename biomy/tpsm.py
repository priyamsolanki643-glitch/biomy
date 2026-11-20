"""Compatibility facade for biomy.tpsm module.

Exposes backward-compatible interfaces mapping to the research-grade biomy.tpsm package.
"""

from .tpsm import (
    TPSMEngine,
    TemporalPKStateMachine,
    DoseEvent,
    TemporalDrug,
    TemporalHerb,
    SimulationResult as TPSMResult,
    TemporalPoint,
    DrugParameters,
    HerbParameters,
    CompartmentMode,
    EliminationMode,
    TPSMError,
    TPSMValidationError,
    TPSMSolverError,
)

__all__ = [
    "TPSMEngine",
    "TemporalPKStateMachine",
    "DoseEvent",
    "TemporalDrug",
    "TemporalHerb",
    "TPSMResult",
    "TemporalPoint",
    "DrugParameters",
    "HerbParameters",
    "CompartmentMode",
    "EliminationMode",
    "TPSMError",
    "TPSMValidationError",
    "TPSMSolverError",
]
