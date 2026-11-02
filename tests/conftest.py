"""Shared fixtures.

The pipeline tests replay captured ReacNetGenerator output through a stub
binary rather than calling the real thing.  ReacNetGenerator is not
deterministic in which SMILES it writes for a given fragment (see
``docs/reproducibility.md``), so a test that invoked it could not assert on
species labels at all.
"""

from __future__ import annotations

import io
import shutil
import subprocess
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path

import pytest

from ensemblemtd.cli import main

DATA = Path(__file__).parent / "data"
REACNET_FIXTURES = DATA / "reacnet"

STUB_SOURCE = '''#!/usr/bin/env python3
"""Replays captured ReacNetGenerator output for one trajectory."""
import os
import shutil
import sys

fixtures = os.environ["REACNET_FIXTURES"]
argv = sys.argv[1:]
name = argv[argv.index("-i") + 1]
for extension in ("json", "species"):
    shutil.copy(os.path.join(fixtures, f"{name}.{extension}"), f"{name}.{extension}")
shutil.copy(os.path.join(fixtures, "template.html"), f"{name}.html")
'''


@pytest.fixture(scope="session")
def obabel() -> str:
    """Path to Open Babel, skipping the test if it is not installed."""
    found = shutil.which("obabel")
    if not found:
        pytest.skip("Open Babel is required for the InChI collapse")
    return found


@pytest.fixture(scope="session")
def reacnet_stub(tmp_path_factory) -> Path:
    stub = tmp_path_factory.mktemp("stub") / "reacnet_stub.py"
    stub.write_text(STUB_SOURCE)
    stub.chmod(0o755)
    return stub


@pytest.fixture
def stub_env(monkeypatch) -> None:
    monkeypatch.setenv("REACNET_FIXTURES", str(REACNET_FIXTURES))


@pytest.fixture
def ensemble_inputs() -> list[str]:
    """The three healthy synthetic trajectories."""
    return [str(DATA / f"k0.5_a0.6run{n}.trj") for n in (1, 2, 3)]


@dataclass
class CliResult:
    returncode: int
    stdout: str
    stderr: str


def run_cli(argv: list[str]) -> CliResult:
    """Invoke the aggregator in-process and capture what it printed."""
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = main(list(argv))
    return CliResult(returncode=code, stdout=out.getvalue(), stderr=err.getvalue())


def run_cli_subprocess(argv: list[str]) -> subprocess.CompletedProcess:
    """Invoke the installed console script the way a user would."""
    return subprocess.run(
        ["ensemble-mtd-aggregate", *argv], capture_output=True, text=True
    )
