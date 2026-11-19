import os
from typing import Protocol


class TextGenerator(Protocol):
    def generate(self, prompt: str, max_new_tokens: int = 256) -> str:
        ...


class BioMistralCore:
    """Thin adapter for a local or Kaggle BioMistral model."""

    def __init__(self, model_id: str = "BioMistral/BioMistral-7B") -> None:
        self.model_id = model_id

    def explain(self, structured_result: dict, generator: TextGenerator | None = None) -> str:
        prompt = (
            "Use only the supplied facts. Do not calculate new doses or make diagnosis claims.\n"
            f"Structured evidence: {structured_result}\n"
            "Give an educational summary, evidence limits, and a request for professional review."
        )
        if generator is not None:
            return generator.generate(prompt)
        return (
            "Local model not loaded. The deterministic findings above remain the source of truth.\n"
            "Consult a qualified clinician or pharmacist before clinical use."
        )
