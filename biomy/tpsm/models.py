"""Data models and typed schemas for the Temporal Pharmacokinetic State Machine (TPSM).

Research-grade temporal PK prototype with illustrative and sourced parameter records.
Not validated for individual clinical dosing or treatment decisions.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple


class CompartmentMode(str, Enum):
    ONE_COMPARTMENT = "1-compartment"
    TWO_COMPARTMENT = "2-compartment"


class EliminationMode(str, Enum):
    LINEAR = "linear"
    MICHAELIS_MENTEN = "michaelis-menten"


class PathwayType(str, Enum):
    CYP_ENZYME = "cyp_enzyme"
    TRANSPORTER = "transporter"


class ParameterStatus(str, Enum):
    SOURCED = "sourced"
    ESTIMATED = "estimated"
    ILLUSTRATIVE = "illustrative"


class EvidenceLevel(str, Enum):
    CLINICAL_TRIAL = "clinical_trial"
    IN_VIVO_HUMAN = "in_vivo_human"
    IN_VITRO = "in_vitro"
    OBSERVATIONAL = "observational"
    THEORETICAL = "theoretical"
    UNKNOWN = "unknown"


class TPSMError(Exception):
    """Base exception for all TPSM errors."""
    pass


class TPSMValidationError(TPSMError, ValueError):
    """Raised when simulation parameters or inputs fail validation."""
    pass


class TPSMSolverError(TPSMError):
    """Raised when the numerical ODE solver encounters an error."""
    pass


@dataclass(frozen=True)
class UncertaintyRange:
    """Quantitative uncertainty bounds for kinetic parameters."""
    low: Optional[float] = None
    high: Optional[float] = None
    unit: str = ""
    confidence: Optional[str] = None
    notes: str = ""


@dataclass(frozen=True)
class ProvenanceRecord:
    """Provenance and audit tracking for empirical and modeled parameters."""
    source_reference: str
    evidence_level: EvidenceLevel = EvidenceLevel.UNKNOWN
    parameter_status: ParameterStatus = ParameterStatus.ILLUSTRATIVE
    uncertainty_range: Optional[UncertaintyRange] = None
    preparation: str = ""
    population_assumptions: str = "Healthy Adult (70 kg reference)"
    limitations: str = "In vitro to in vivo extrapolation; inter-individual variability applies."


@dataclass(frozen=True)
class EnzymePathwayParam:
    """Parameters for a single enzyme or transporter interaction pathway."""
    pathway_name: str  # e.g., "CYP3A4", "CYP2C9", "OATP1A2", "P-gp"
    pathway_type: PathwayType = PathwayType.CYP_ENZYME
    fm: float = 0.0  # Fraction metabolized/cleared by this pathway
    is_transporter: bool = False
    reversible_ki_u_mg_l: Optional[float] = None  # Unbound Ki for reversible inhibition (mg/L)
    kinact_h: Optional[float] = None  # Maximum inactivation rate constant (1/h)
    ki_u_mbi_mg_l: Optional[float] = None  # Unbound perpetrator conc producing half-max inactivation (mg/L)
    recovery_half_life_h: Optional[float] = None  # De novo enzyme synthesis half-life (h)
    provenance: Optional[ProvenanceRecord] = None


@dataclass(frozen=True)
class DrugParameters:
    """Kinetic parameters for a pharmaceutical agent (victim or substrate)."""
    name: str
    role: str = "victim"
    formulation: str = "oral_solid"
    bioavailability: float = 1.0  # F in (0, 1]
    ka_h: float = 1.0  # Absorption rate constant (1/h)
    volume_central_l: float = 10.0  # Central volume Vc / Vd (L)
    volume_peripheral_l: Optional[float] = None  # Peripheral volume Vp (L) for 2-comp
    inter_compartment_clearance_l_h: Optional[float] = None  # Inter-compartment clearance Q (L/h)
    clearance_l_h: float = 1.0  # Baseline systemic clearance CL (L/h)
    fraction_unbound: float = 1.0  # fu in (0, 1]
    half_life_h: Optional[float] = None
    elimination_mode: EliminationMode = EliminationMode.LINEAR
    vmax_mg_h: Optional[float] = None  # Vmax (mg/h) for Michaelis-Menten
    km_mg_l: Optional[float] = None  # Km (mg/L) for Michaelis-Menten
    fm_pathways: Dict[str, float] = field(default_factory=dict)
    pathway_types: Dict[str, str] = field(default_factory=dict)
    therapeutic_min_mg_l: Optional[float] = None  # Minimum effective concentration (MEC)
    toxicity_threshold_mg_l: Optional[float] = None  # Minimum toxic concentration (MTC)
    age_restrictions: str = "Requires clinical review"
    pregnancy_restrictions: str = ""
    renal_hepatic_cautions: str = ""
    side_effects: tuple[str, ...] = ()
    provenance: Optional[ProvenanceRecord] = None


@dataclass(frozen=True)
class HerbParameters:
    """Kinetic parameters for an herbal or botanical perpetrator."""
    name: str
    preparation: str = "standardized_extract"
    active_constituent: str = ""
    dose_unit: str = "mg"
    bioavailability: float = 0.5
    ka_h: float = 1.0
    volume_central_l: float = 50.0
    clearance_l_h: float = 10.0
    fraction_unbound: float = 0.1
    pathway_effects: Dict[str, EnzymePathwayParam] = field(default_factory=dict)
    onset_h: Optional[float] = None
    duration_h: Optional[float] = None
    age_restrictions: str = "Requires clinical review"
    pregnancy_restrictions: str = ""
    side_effects: tuple[str, ...] = ()
    provenance: Optional[ProvenanceRecord] = None


@dataclass(frozen=True)
class DoseEvent:
    """Discrete administration event of a drug or botanical compound."""
    name: str = ""
    dose_mg: float = 0.0
    time_h: float = 0.0
    route: str = "oral"
    notes: str = ""

    def __init__(
        self,
        name: str = "",
        dose_mg: float = 0.0,
        time_h: float = 0.0,
        route: str = "oral",
        notes: str = "",
        substance: Optional[str] = None,
    ):
        target_name = substance if substance is not None else name
        object.__setattr__(self, "name", target_name)
        object.__setattr__(self, "dose_mg", float(dose_mg))
        object.__setattr__(self, "time_h", float(time_h))
        object.__setattr__(self, "route", str(route))
        object.__setattr__(self, "notes", str(notes))

    @property
    def substance(self) -> str:
        return self.name


@dataclass(frozen=True)
class TemporalPoint:
    """Snapshot of continuous state variables at a specific point in time."""
    time_h: float
    concentrations_mg_l: Dict[str, float]
    amounts_mg: Dict[str, float]
    active_enzyme_fraction: Dict[str, float]
    toxicity_flags: Dict[str, bool]
    tmr_values: Dict[str, Optional[float]]


@dataclass(frozen=True)
class SimulationResult:
    """Comprehensive output record emitted by the TPSM engine."""
    solver_success: bool
    solver_method: str
    solver_message: str
    time_points: List[float]
    concentration_curves: Dict[str, List[float]]  # compound_name -> concentrations (mg/L)
    enzyme_activity_curves: Dict[str, List[float]]  # enzyme_name -> activity fraction (0.0 to 1.0)
    auc: Dict[str, float]  # compound_name -> AUC (mg*h/L)
    c_max: Dict[str, float]  # compound_name -> peak concentration (mg/L)
    c_min: Dict[str, float]  # compound_name -> trough concentration (mg/L)
    first_toxicity_breach_time: Dict[str, Optional[float]]  # compound_name -> first breach time (h) or None
    therapeutic_margin_ratio_curves: Dict[str, List[Optional[float]]]
    final_state: Dict[str, float]
    event_log: List[Dict[str, Any]]
    parameter_provenance: Dict[str, Any]
    warnings: List[str]
    assumptions: List[str]
    uncertainty_metadata: Dict[str, Any]

    @property
    def success(self) -> bool:
        """Convenience alias for solver_success."""
        return self.solver_success

    @property
    def solver_status(self) -> str:
        """Convenience status string."""
        return "converged" if self.solver_success else "failed"

    @property
    def points(self) -> List["TemporalPoint"]:
        """Backward-compatible access to time-point snapshots."""
        pts: List[TemporalPoint] = []
        for i, t in enumerate(self.time_points):
            concs = {}
            for c in self.concentration_curves:
                if i < len(self.concentration_curves[c]):
                    val = self.concentration_curves[c][i]
                    concs[c] = val
                    concs[c.lower()] = val
            enzs = {e: self.enzyme_activity_curves[e][i] for e in self.enzyme_activity_curves if i < len(self.enzyme_activity_curves[e])}
            tmrs = {c: self.therapeutic_margin_ratio_curves[c][i] for c in self.therapeutic_margin_ratio_curves if i < len(self.therapeutic_margin_ratio_curves[c])}
            pts.append(
                TemporalPoint(
                    time_h=t,
                    concentrations_mg_l=concs,
                    amounts_mg={},
                    active_enzyme_fraction=enzs,
                    toxicity_flags={},
                    tmr_values=tmrs,
                )
            )
        return pts

    def series(self, name: str) -> List[Tuple[float, float]]:
        """Return (time_h, concentration_mg_l) tuples for a named compound."""
        key = name.lower()
        for k, curve in self.concentration_curves.items():
            if k.lower() == key:
                return list(zip(self.time_points, curve))
        return [(t, 0.0) for t in self.time_points]

    def first_toxicity_time(self, name: str) -> Optional[float]:
        """Return first toxicity breach time for a named compound."""
        key = name.lower()
        for k, t_breach in self.first_toxicity_breach_time.items():
            if k.lower() == key:
                return t_breach
        return None
