# Contributing

Bug reports, questions about the method, and patches are all welcome. Open an
issue at https://github.com/adilpsn/EnsembleMTD/issues.

## Getting set up

```bash
git clone https://github.com/adilpsn/EnsembleMTD.git
cd EnsembleMTD
pip install -e ".[dev,plots]"
pytest -m unit          # seconds, needs nothing but Python
pytest                  # adds the end-to-end tests, which need Open Babel
```

CI runs the unit tests on Python 3.9–3.13, the full suite with Open Babel and
Graphviz installed, `shellcheck --severity=warning` on the shell scripts, and a
package build. Running `pytest` and `shellcheck reactor/*.sh` locally will catch
almost everything before you push.

## Things worth knowing before you change something

**ReacNetGenerator is not deterministic.** It can write a different SMILES for
the same fragment between runs, so a test that invoked it could not assert on
species labels. The end-to-end tests replay captured output through a stub
binary passed via `--reacnet-bin`; see `tests/conftest.py`. If you are changing
the aggregation, that same stub is how to check your change against the
previous behaviour — `docs/reproducibility.md` explains the procedure.

**The output filenames and column names are load-bearing.** Published figures
and other people's analysis scripts read them. Adding a column is fine;
renaming or reordering one is a breaking change.

**Prevalence is a binary per-run indicator.** If a change makes Π_i depend on
how often or how long a species was present, that is a different metric, not a
bug fix. `docs/interpretation.md` states what these numbers do and do not mean,
and it should stay true.

**Nothing in the repository may contain site-specific configuration.** No
scheduler directives, no `module load`, no absolute paths into somebody's home
directory. `reactor/submit.example.sbatch` is the single place a user adapts.

## Style

Existing code sticks to PEP 8 with type annotations on function signatures,
frozen dataclasses for anything that travels between modules, and modules under
about 400 lines. Comments explain *why* a thing is done — the kpush `/4`, why Li
is stripped before connectivity is inferred, why HMM filtering is off — rather
than restating the code.

## The reactor scripts

They need `crest` and `xtb`, so CI can only lint them. If you change one, run
`bash reactor/smoke_test.sh` somewhere with those installed and say in the pull
request that you did; it checks the generated `rcontrol`, a two-member ensemble,
and that resuming tops an ensemble up rather than overwriting it.
