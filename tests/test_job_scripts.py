"""The Slurm job scripts are code too, and nothing was checking them.

Two defects motivated this file, both found the expensive way -- by submitting
a job and reading the log afterwards:

  * A gate printed FAIL and exited 0, so a --dependency=afterok job downstream
    launched against a failing gate.
  * A patch to a gate left a literal newline inside a string, so the comparison
    block died with a SyntaxError after the pipeline stage had already run.

Both are checkable statically, in milliseconds, with no cluster.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

JOBS = Path(__file__).resolve().parents[1] / "jobs"
HEREDOC = re.compile(r"<<'PY'\n(.*?)\n(?:PY)(?:\n|$)", re.S)

SCRIPTS = sorted(JOBS.glob("*.sbatch"))
GATES = [p for p in SCRIPTS if p.name.startswith("gate_")]


def _blocks(path: Path):
    return HEREDOC.findall(path.read_text())


def test_there_are_job_scripts_to_check():
    assert SCRIPTS, "no .sbatch files found -- the glob or layout changed"
    assert GATES, "no gate_*.sbatch files found"


@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: p.name)
def test_embedded_python_compiles(path):
    """Every `python - <<'PY' ... PY` block must parse."""
    for i, block in enumerate(_blocks(path)):
        compile(block, f"{path.name}[block {i}]", "exec")


@pytest.mark.parametrize("path", GATES, ids=lambda p: p.name)
def test_a_gate_can_fail(path):
    """A gate that cannot exit nonzero cannot gate anything.

    Either the script exits nonzero itself, or it delegates to a comparison
    helper that does.
    """
    text = path.read_text()
    exits = "sys.exit(" in text or re.search(r"\bexit [1-9]", text)
    delegates = "compare_one_figure" in text or "compare_figures" in text
    assert exits or delegates, (
        f"{path.name} reports a verdict but always exits 0; a "
        "--dependency=afterok job would launch against a failing gate"
    )


@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: p.name)
def test_declares_partition_and_logs(path):
    text = path.read_text()
    assert "--partition=" in text, f"{path.name} has no #SBATCH --partition"
    assert "--output=" in text, f"{path.name} does not capture stdout"
