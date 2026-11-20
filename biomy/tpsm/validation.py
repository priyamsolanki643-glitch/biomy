"""Parameter and state validation for clinical simulation safety."""

from math import isfinite
from typing import Dict, Iterable, List, Optional

from .models import CompartmentMode, DrugParameters, EliminationMode, HerbParameters, TPSMValidationError


class SimulationValidator:
    """Performs rigorous pre-simulation checks on kinetic parameters and boundary conditions."""

    @staticmethod
    def validate_drug(drug: DrugParameters, mode: CompartmentMode) -> None:
        """Validate kinetic parameters for a pharmaceutical substrate."""
        if not drug.name or not drug.name.strip():
            raise TPSMValidationError("Drug name cannot be empty")

        if not isfinite(drug.bioavailability) or not (0.0 < drug.bioavailability <= 1.0):
            raise TPSMValidationError(f"{drug.name}: bioavailability F must be in (0, 1], got {drug.bioavailability}")

        if not isfinite(drug.ka_h) or drug.ka_h <= 0:
            raise TPSMValidationError(f"{drug.name}: absorption rate ka must be positive, got {drug.ka_h}")

        if not isfinite(drug.volume_central_l) or drug.volume_central_l <= 0:
            raise TPSMValidationError(f"{drug.name}: central volume Vc must be positive, got {drug.volume_central_l}")

        if not isfinite(drug.clearance_l_h) or drug.clearance_l_h < 0:
            raise TPSMValidationError(f"{drug.name}: clearance CL must be non-negative, got {drug.clearance_l_h}")

        if not isfinite(drug.fraction_unbound) or not (0.0 < drug.fraction_unbound <= 1.0):
            raise TPSMValidationError(f"{drug.name}: fraction unbound fu must be in (0, 1], got {drug.fraction_unbound}")

        if drug.elimination_mode == EliminationMode.MICHAELIS_MENTEN:
            if drug.vmax_mg_h is None or drug.vmax_mg_h <= 0 or not isfinite(drug.vmax_mg_h):
                raise TPSMValidationError(f"{drug.name}: Vmax must be strictly positive for Michaelis-Menten elimination")
            if drug.km_mg_l is None or drug.km_mg_l <= 0 or not isfinite(drug.km_mg_l):
                raise TPSMValidationError(f"{drug.name}: Km must be strictly positive for Michaelis-Menten elimination")

        if mode == CompartmentMode.TWO_COMPARTMENT:
            if drug.volume_peripheral_l is not None and (drug.volume_peripheral_l <= 0 or not isfinite(drug.volume_peripheral_l)):
                raise TPSMValidationError(f"{drug.name}: peripheral volume Vp must be positive for 2-compartment mode")
            if drug.inter_compartment_clearance_l_h is not None and (
                drug.inter_compartment_clearance_l_h <= 0 or not isfinite(drug.inter_compartment_clearance_l_h)
            ):
                raise TPSMValidationError(f"{drug.name}: inter-compartment clearance Q must be positive")

        # Validate therapeutic and toxicity thresholds
        if drug.therapeutic_min_mg_l is not None:
            if not isfinite(drug.therapeutic_min_mg_l) or drug.therapeutic_min_mg_l <= 0:
                raise TPSMValidationError(f"{drug.name}: therapeutic_min_mg_l must be strictly positive, got {drug.therapeutic_min_mg_l}")
        if drug.toxicity_threshold_mg_l is not None:
            if not isfinite(drug.toxicity_threshold_mg_l) or drug.toxicity_threshold_mg_l <= 0:
                raise TPSMValidationError(f"{drug.name}: toxicity_threshold_mg_l must be strictly positive, got {drug.toxicity_threshold_mg_l}")
        if drug.therapeutic_min_mg_l is not None and drug.toxicity_threshold_mg_l is not None:
            if drug.therapeutic_min_mg_l >= drug.toxicity_threshold_mg_l:
                raise TPSMValidationError(
                    f"{drug.name}: therapeutic_min_mg_l ({drug.therapeutic_min_mg_l}) must be strictly less than toxicity_threshold_mg_l ({drug.toxicity_threshold_mg_l})"
                )

        # Validate fractions metabolized
        for enz, fm_val in drug.fm_pathways.items():
            if not isfinite(fm_val) or fm_val < 0.0:
                raise TPSMValidationError(f"{drug.name}: fm pathway '{enz}' must be non-negative, got {fm_val}")
        fm_sum = sum(drug.fm_pathways.values())
        if fm_sum > 1.0001:
            raise TPSMValidationError(f"{drug.name}: sum of fm pathways exceeds 1.0 (got {fm_sum})")

    @staticmethod
    def validate_compounds(drugs: Iterable[DrugParameters], herbs: Iterable[HerbParameters] = ()) -> None:
        """Validate uniqueness of compound names across drugs and herbs."""
        seen_names: set[str] = set()
        for d in drugs:
            d_name = d.name.strip().lower()
            if d_name in seen_names:
                raise TPSMValidationError(f"Duplicate compound name '{d.name}' detected in drugs")
            seen_names.add(d_name)
        for h in herbs:
            h_name = h.name.strip().lower()
            if h_name in seen_names:
                raise TPSMValidationError(f"Duplicate compound name '{h.name}' detected across drugs and herbs")
            seen_names.add(h_name)

    @staticmethod
    def validate_herb(herb: HerbParameters) -> None:
        """Validate kinetic parameters for an herbal perpetrator."""
        if not herb.name or not herb.name.strip():
            raise TPSMValidationError("Herb name cannot be empty")

        if not isfinite(herb.bioavailability) or not (0.0 < herb.bioavailability <= 1.0):
            raise TPSMValidationError(f"{herb.name}: bioavailability F must be in (0, 1], got {herb.bioavailability}")

        if not isfinite(herb.ka_h) or herb.ka_h <= 0:
            raise TPSMValidationError(f"{herb.name}: absorption rate ka must be positive, got {herb.ka_h}")

        if not isfinite(herb.volume_central_l) or herb.volume_central_l <= 0:
            raise TPSMValidationError(f"{herb.name}: volume Vc must be positive, got {herb.volume_central_l}")

        if not isfinite(herb.clearance_l_h) or herb.clearance_l_h < 0:
            raise TPSMValidationError(f"{herb.name}: clearance CL must be non-negative, got {herb.clearance_l_h}")

        for enz, pathway in herb.pathway_effects.items():
            if pathway.kinact_h is not None and (pathway.kinact_h < 0 or not isfinite(pathway.kinact_h)):
                raise TPSMValidationError(f"{herb.name} ({enz}): kinact must be non-negative, got {pathway.kinact_h}")
            if pathway.kinact_h is not None and pathway.kinact_h > 0:
                if pathway.ki_u_mbi_mg_l is None or pathway.ki_u_mbi_mg_l <= 0 or not isfinite(pathway.ki_u_mbi_mg_l):
                    raise TPSMValidationError(f"{herb.name} ({enz}): MBI KIu must be positive when kinact > 0, got {pathway.ki_u_mbi_mg_l}")
            if pathway.reversible_ki_u_mg_l is not None and (
                pathway.reversible_ki_u_mg_l <= 0 or not isfinite(pathway.reversible_ki_u_mg_l)
            ):
                raise TPSMValidationError(f"{herb.name} ({enz}): reversible Ki must be positive, got {pathway.reversible_ki_u_mg_l}")
