"""Path and parameter resolution.

Every path used by the pipelines resolves through here, so that running the
analysis somewhere other than our OAK directory is a config change rather than
an edit to tracked code. Resolution order:

    1. WGS_DATA_DIR / WGS_RESULTS_DIR environment variables
    2. the config file named by WGS_CONFIG
    3. config.yaml at the repository root

PyYAML is optional. If it is unavailable the loader falls back to a minimal
parser sufficient for this file's structure, so a missing dependency cannot
silently change which inputs an analysis reads.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config.yaml"


def _minimal_yaml(text: str) -> dict:
    """Parse the restricted YAML subset used by config.yaml.

    Handles two levels of mapping, inline lists, and scalars. Deliberately
    narrow -- anything it cannot parse raises rather than guessing.
    """
    root: dict = {}
    section: dict | None = None

    for raw in text.splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue

        indented = line.startswith((" ", "\t"))
        if ":" not in line:
            raise ValueError(f"cannot parse config line: {raw!r}")
        key, _, value = line.strip().partition(":")
        key, value = key.strip(), value.strip()

        if not indented:
            if value == "":
                section = {}
                root[key] = section
            else:
                root[key] = _scalar(value)
                section = None
        else:
            if section is None:
                raise ValueError(f"indented entry outside a section: {raw!r}")
            section[key] = _scalar(value)

    return root


def _scalar(value: str):
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [_scalar(v.strip()) for v in inner.split(",")]
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    low = value.lower()
    if low in ("true", "false"):
        return low == "true"
    if low in ("null", "none", "~"):
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def _load_raw(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore
    except ImportError:
        return _minimal_yaml(text)
    return yaml.safe_load(text)


class Config:
    """Resolved paths and analysis parameters."""

    def __init__(self, path: Path | str | None = None):
        self.path = Path(path or os.environ.get("WGS_CONFIG") or DEFAULT_CONFIG)
        if not self.path.exists():
            raise FileNotFoundError(
                f"config not found at {self.path}. Set WGS_CONFIG, or copy "
                f"config.yaml from the repository root."
            )
        raw = _load_raw(self.path)

        self.data_dir = Path(
            os.environ.get("WGS_DATA_DIR") or raw["data_dir"])
        self.results_dir = Path(
            os.environ.get("WGS_RESULTS_DIR") or raw["results_dir"])
        self._inputs = raw.get("inputs", {})
        self._derived = raw.get("derived", {})
        self.params = raw.get("params", {})

    # -- path accessors ----------------------------------------------------
    def data(self, key: str) -> Path:
        """Path to a raw input, by config key."""
        return self._resolve(self._inputs, key, self.data_dir)

    def derived(self, key: str) -> Path:
        """Path to a derived artefact, by config key."""
        return self._resolve(self._derived, key, self.results_dir)

    def result(self, filename: str) -> Path:
        """Path to an output file under results_dir."""
        return self.results_dir / filename

    @staticmethod
    def _resolve(mapping: dict, key: str, root: Path) -> Path:
        if key not in mapping:
            raise KeyError(f"{key!r} is not defined in the config "
                           f"(available: {sorted(mapping)})")
        value = Path(mapping[key])
        return value if value.is_absolute() else root / value

    def param(self, key: str):
        if key not in self.params:
            raise KeyError(f"{key!r} is not defined under params in the config")
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
