"""Structured serialization for BioMistral and clinical explanation interfaces."""

from typing import Any, Dict
from .models import SimulationResult


class TPSMSerializer:
    """Serializes numerical TPSM results into an auditable prompt fact matrix for BioMistral."""

    @staticmethod
    def to_biomistral_context(result: SimulationResult) -> Dict[str, Any]:
        """
        Produce a schema-enforced dictionary for LLM clinical explanation conditioning.
        The LLM is constrained to explain these findings without recalculating numbers.
        """
        payload: Dict[str, Any] = {
            "meta": {
                "engine": "BioMy Temporal Pharmacokinetic State Machine (TPSM)",
                "status": "computed_biological_state",
                "solver_method": result.solver_method,
                "solver_success": result.solver_success,
                "solver_message": result.solver_message,
            },
            "system_instruction_for_llm": (
                "You are an auditable clinical intelligence assistant. The numerical findings below "
                "are deterministically computed by the TPSM ODE engine and are the authoritative source of truth. "
                "DO NOT recalculate concentrations, DO NOT invent alternative values, and DO NOT propose "
                "unsupervised dosage adjustments. Explain the pharmacological mechanisms, timing of risks, "
                "and advise professional clinician/pharmacist review."
            ),
            "exposure_metrics": {
                "compounds": list(result.concentration_curves.keys()),
                "auc_mg_h_l": result.auc,
                "c_max_mg_l": result.c_max,
                "c_min_mg_l": result.c_min,
                "first_toxicity_breach_timestamp_h": result.first_toxicity_breach_time,
            },
            "enzyme_dynamics": {
                "enzymes_monitored": list(result.enzyme_activity_curves.keys()),
                "final_active_fraction": {
                    enz: result.enzyme_activity_curves[enz][-1]
                    for enz in result.enzyme_activity_curves
                },
            },
            "therapeutic_safety": {
                "latest_tmr": {
                    comp: result.therapeutic_margin_ratio_curves[comp][-1]
                    for comp in result.therapeutic_margin_ratio_curves
                    if comp in result.concentration_curves
                },
                "breach_detected": any(
                    t is not None for t in result.first_toxicity_breach_time.values()
                ),
            },
            "parameter_provenance": result.parameter_provenance,
            "clinical_warnings": result.warnings,
            "modeling_assumptions": result.assumptions,
            "uncertainty_metadata": result.uncertainty_metadata,
        }
        return payload

    @staticmethod
    def to_formatted_prompt(result: SimulationResult, patient_context: str = "") -> str:
        """Generate a complete text prompt for BioMistral consultation generation."""
        context = TPSMSerializer.to_biomistral_context(result)
        lines = [
            "=== BIOMY DETERMINISTIC TPSM FACT MATRIX ===",
            f"Solver: {context['meta']['solver_method']} (Success: {context['meta']['solver_success']})",
            "",
            "EXPOSURE METRICS:",
        ]
        for comp in context["exposure_metrics"]["compounds"]:
            auc = context["exposure_metrics"]["auc_mg_h_l"].get(comp, "N/A")
            cmax = context["exposure_metrics"]["c_max_mg_l"].get(comp, "N/A")
            breach = context["exposure_metrics"]["first_toxicity_breach_timestamp_h"].get(comp)
            breach_str = f"BREACH at t={breach}h" if breach is not None else "Within therapeutic limits"
            lines.append(f"  - {comp}: Cmax={cmax} mg/L | AUC={auc} mg*h/L | Status: {breach_str}")

        lines.append("\nENZYME ACTIVITY FRACTION (END OF WINDOW):")
        for enz, val in context["enzyme_dynamics"]["final_active_fraction"].items():
            lines.append(f"  - {enz}: {val * 100.0:.1f}% active")

        if context["clinical_warnings"]:
            lines.append("\nCLINICAL CAUTIONS:")
            for warn in context["clinical_warnings"]:
                lines.append(f"  - [WARNING] {warn}")

        lines.append("\nMODEL ASSUMPTIONS:")
        for asm in context["modeling_assumptions"]:
            lines.append(f"  - {asm}")

        if patient_context:
            lines.append(f"\nPATIENT CLINICAL CONTEXT:\n{patient_context}")

        lines.append("\n=== TASK FOR BIOMISTRAL ===")
        lines.append(context["system_instruction_for_llm"])
        return "\n".join(lines)

    @staticmethod
    def to_llm_prompt_context(result: SimulationResult, patient_context: str = "") -> str:
        """Alias for to_formatted_prompt for clean LLM prompt context generation."""
        return TPSMSerializer.to_formatted_prompt(result, patient_context=patient_context)
