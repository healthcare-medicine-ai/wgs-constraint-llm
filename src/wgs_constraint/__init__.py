"""Shared components for the WGS constraint / unified meta-regression analyses.

Importing from here rather than copy-pasting into notebooks is the point: the
GERP annotation existed in three notebooks with the same coordinate defect in
each, and correcting one did nothing for the others.

Typical use from a notebook::

    import sys; sys.path.insert(0, "src")
    from wgs_constraint import get_config, annotate_with_gerp
    cfg = get_config()
    merged = annotate_with_gerp(pred, cfg.data("gerp_bigwig"))
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
