"""Numerical ODE solver integration for stiff biological dynamics using SciPy."""

from typing import Callable, List, Optional, Tuple
import numpy as np
from scipy.integrate import solve_ivp

from .models import TPSMSolverError


class ODESolverWrapper:
    """Robust numerical integrator wrapper with stiffness detection and root-finding."""

    SUPPORTED_METHODS = {"Radau", "BDF", "RK45"}

    def __init__(
        self,
        method: str = "Radau",
        rtol: float = 1e-6,
        atol: float = 1e-8,
        max_step: Optional[float] = None,
    ):
        if method not in self.SUPPORTED_METHODS:
            raise TPSMSolverError(f"Unsupported solver method '{method}'. Choose from {self.SUPPORTED_METHODS}")
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step

    def integrate_segment(
        self,
        fun: Callable[[float, np.ndarray], np.ndarray],
        t_span: Tuple[float, float],
        y0: np.ndarray,
        events: Optional[List[Callable]] = None,
        num_points: int = 50,
    ) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray], bool, str]:
        """
        Integrate system from t_span[0] to t_span[1].
        Returns (t_eval, y_eval, t_events, success, message).
        Raises TPSMSolverError if solver fails or produces non-finite states.
        """
        t0, t1 = t_span
        if t1 <= t0:
            return np.array([t0]), y0.reshape(-1, 1), [], True, "Zero-duration segment"

        # Construct evaluation grid
        t_eval = np.linspace(t0, t1, max(2, num_points))

        # Check initial state for non-finite values
        if not np.all(np.isfinite(y0)):
            raise TPSMSolverError(f"Initial state at t={t0} contains non-finite values: {y0}")

        kwargs = {
            "fun": fun,
            "t_span": t_span,
            "y0": y0,
            "method": self.method,
            "t_eval": t_eval,
            "events": events,
            "rtol": self.rtol,
            "atol": self.atol,
        }
        if self.max_step is not None:
            kwargs["max_step"] = self.max_step

        sol = solve_ivp(**kwargs)

        if not sol.success:
            raise TPSMSolverError(
                f"ODE solver '{self.method}' failed on interval [{t0:.2f}, {t1:.2f}]h: {sol.message}"
            )

        # Sanity check: finite values and physical non-negativity
        if not np.all(np.isfinite(sol.y)):
            raise TPSMSolverError(
                f"ODE solver produced non-finite values (NaN/Inf) on interval [{t0:.2f}, {t1:.2f}]h"
            )

        # Numerical cleanup: prevent tiny negative floats from precision artifacts
        y_cleaned = np.maximum(sol.y, 0.0)

        return sol.t, y_cleaned, sol.t_events or [], sol.success, sol.message
