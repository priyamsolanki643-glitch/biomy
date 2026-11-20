"""Pharmacokinetic metrics, exposure accumulation, and toxicity threshold evaluation."""

from typing import Dict, List, Optional, Tuple
import numpy as np

from .models import DrugParameters


class PKMetricsCalculator:
    """Calculates exposure metrics (AUC, Cmax, Cmin, TMR) from concentration trajectories."""

    @staticmethod
    def calculate_auc(time_points: np.ndarray, concentrations: np.ndarray) -> float:
        """Calculate Area Under the Curve (AUC, mg*h/L) using trapezoidal rule."""
        if len(time_points) < 2:
            return 0.0
        # Use np.trapezoid for numpy >= 2.0 or np.trapz fallback
        if hasattr(np, "trapezoid"):
            val = float(np.trapezoid(concentrations, time_points))
        else:
            val = float(np.trapz(concentrations, time_points))
        return max(0.0, val)

    @staticmethod
    def calculate_cmax_cmin(concentrations: np.ndarray) -> Tuple[float, float]:
        """Calculate peak (Cmax) and trough (Cmin) concentration."""
        if len(concentrations) == 0:
            return 0.0, 0.0
        return float(np.max(concentrations)), float(np.min(concentrations))

    @staticmethod
    def calculate_tmr_curve(
        concentrations: np.ndarray,
        drug: DrugParameters,
    ) -> List[Optional[float]]:
        """
        Calculate Therapeutic Margin Ratio (TMR) time series:
        TMR(t) = (Ctox - C(t)) / (Ctox - CMEC)
        If CMEC is missing: TMR(t) = (Ctox - C(t)) / Ctox
        If Ctox is missing: returns [None, ...]
        """
        if drug.toxicity_threshold_mg_l is None:
            return [None] * len(concentrations)

        c_tox = float(drug.toxicity_threshold_mg_l)
        c_mec = float(drug.therapeutic_min_mg_l) if drug.therapeutic_min_mg_l is not None else 0.0

        denom = c_tox - c_mec
        if abs(denom) < 1e-9:
            denom = c_tox if c_tox > 0 else 1.0

        tmr_list: List[Optional[float]] = []
        for c in concentrations:
            tmr = (c_tox - float(c)) / denom
            tmr_list.append(round(tmr, 4))
        return tmr_list

    @staticmethod
    def detect_first_toxicity_breach(
        time_points: np.ndarray,
        concentrations: np.ndarray,
        threshold: Optional[float],
        solver_events: Optional[List[np.ndarray]] = None,
    ) -> Optional[float]:
        """
        Detect the earliest timestamp where concentration crosses upward into toxicity.
        Prefers exact root-finding events from the ODE solver.
        """
        if threshold is None:
            return None

        # 1. Prefer exact solver root-finding events
        if solver_events:
            for ev_array in solver_events:
                if len(ev_array) > 0 and ev_array[0] is not None:
                    return round(float(ev_array[0]), 4)

        # 2. Discrete grid search fallback
        thresh_val = float(threshold)
        for t, c in zip(time_points, concentrations):
            if c >= thresh_val:
                return round(float(t), 4)

        return None
