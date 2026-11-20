"""Catalog helper for loading temporal pharmacokinetic parameters."""

from pathlib import Path
from typing import Any, Dict, Tuple

from .tpsm.models import DrugParameters, HerbParameters
from .tpsm.parameters import ParameterCatalog, DEFAULT_CATALOG_PATH


def load_temporal_catalog(
    path: Path = DEFAULT_CATALOG_PATH,
) -> Tuple[Dict[str, DrugParameters], Dict[str, HerbParameters], Dict[str, Any]]:
    """Load and return (drugs, herbs, metadata) dictionary tuple."""
    cat = ParameterCatalog.from_json_file(path)
    return cat.drugs, cat.herbs, cat.metadata
