"""Research-grade Temporal Pharmacokinetic State Machine (TPSM) package.

Models coupled multi-compartment disposition, enzyme MBI/synthesis,
and therapeutic margin dynamics for clinical decision support research.
"""

from .engine import TPSMEngine, TemporalPKStateMachine
from .events import DoseEvent, EventTimelineProcessor
from .models import (
    CompartmentMode,
    DrugParameters,
    EliminationMode,
    EnzymePathwayParam,
    EvidenceLevel,
    HerbParameters,
    ParameterStatus,
    PathwayType,
    ProvenanceRecord,
    SimulationResult,
    TemporalPoint,
    TPSMError,
    TPSMSolverError,
    TPSMValidationError,
    UncertaintyRange,
)
from .parameters import ParameterCatalog
from .serialization import TPSMSerializer

# Backward-compatibility aliases
TemporalDrug = DrugParameters
TemporalHerb = HerbParameters

__all__ = [
    "TPSMEngine",
    "TemporalPKStateMachine",
    "DoseEvent",
    "EventTimelineProcessor",
    "CompartmentMode",
    "EliminationMode",
    "EnzymePathwayParam",
    "EvidenceLevel",
    "HerbParameters",
    "DrugParameters",
    "TemporalDrug",
    "TemporalHerb",
    "ParameterStatus",
    "PathwayType",
    "ProvenanceRecord",
    "SimulationResult",
    "TemporalPoint",
    "TPSMError",
    "TPSMSolverError",
    "TPSMValidationError",
    "UncertaintyRange",
    "ParameterCatalog",
    "TPSMSerializer",
]
