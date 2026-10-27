"""Paths and analysis parameters, read from config.yaml.

Resolution order for the two root directories:

    1. WGS_DATA_DIR / WGS_RESULTS_DIR environment variables
    2. the config file named by WGS_CONFIG
    3. config.yaml at the repository root

Relative paths are taken relative to the repository root.
"""

from __future__ import annotations

import os
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config.yaml"


def _root(env: str, configured: str) -> Path:
    path = Path(os.environ.get(env) or configured)
    return path if path.is_absolute() else REPO_ROOT / path


class Config:
    """Resolved paths and analysis parameters."""

    def __init__(self, path: Path | str | None = None):
        self.path = Path(path or os.environ.get("WGS_CONFIG") or DEFAULT_CONFIG)
        raw = yaml.safe_load(self.path.read_text(encoding="utf-8"))
        self.data_dir = _root("WGS_DATA_DIR", raw["data_dir"])
        self.results_dir = _root("WGS_RESULTS_DIR", raw["results_dir"])
        self._inputs = raw.get("inputs", {})
        self._derived = raw.get("derived", {})
        self.params = raw.get("params", {})

    def data(self, key: str) -> Path:
        """Path to a raw input, by config key."""
        return self._resolve(self._inputs, key, self.data_dir)

    def derived(self, key: str) -> Path:
        """Path to a derived input, by config key."""
        return self._resolve(self._derived, key, self.results_dir)

    def result(self, filename: str) -> Path:
        """Path to an output file under results_dir."""
        return self.results_dir / filename

    @staticmethod
    def _resolve(mapping: dict, key: str, root: Path) -> Path:
        if key not in mapping:
            raise KeyError(f"{key!r} is not in the config (available: {sorted(mapping)})")
        value = Path(mapping[key])
        return value if value.is_absolute() else root / value

    def param(self, key: str):
        if key not in self.params:
            raise KeyError(f"{key!r} is not under params in the config")
        return self.params[key]

    def describe(self) -> str:
        return (f"config      : {self.path}\n"
                f"data_dir    : {self.data_dir}\n"
                f"results_dir : {self.results_dir}")


_cached: Config | None = None


def get_config(reload: bool = False) -> Config:
    global _cached
    if _cached is None or reload:
        _cached = Config()
    return _cached
