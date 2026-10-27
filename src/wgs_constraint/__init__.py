"""Shared code for the HMM constraint model and the unified meta-regression.

The pipeline stages in pipelines/ import from here, so each step (GERP lookup,
AlphaMissense collapse, Epi25 filtering, the regression) is implemented once.
"""

from .alphamissense import collapse, load_collapsed
from .config import Config, get_config
from .gerp import GERP_COLUMN, annotate_with_gerp, chromosome_sort_key
from .metareg import (
    GROUP_KEY, MODERATORS, fit_per_gene, haldane_effect_sizes,
    prepare_regression_input,
)

__all__ = [
    "Config", "get_config",
    "annotate_with_gerp", "chromosome_sort_key", "GERP_COLUMN",
    "load_collapsed", "collapse",
    "haldane_effect_sizes", "prepare_regression_input", "fit_per_gene",
    "MODERATORS", "GROUP_KEY",
]

__version__ = "0.3.0"
