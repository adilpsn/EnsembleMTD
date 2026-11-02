# Contributing

Report bugs or suggest changes through
[GitHub issues](https://github.com/adilpsn/EnsembleMTD/issues).

## Setup

```bash
git clone https://github.com/adilpsn/EnsembleMTD.git
cd EnsembleMTD
pip install -e ".[dev,plots]"
pytest -m unit
```

The full suite requires Open Babel. CI also runs shellcheck and a package build.
For changes to the reactor scripts, run `bash reactor/smoke_test.sh` with CREST
and xTB installed and include the result in the pull request.

## Code and outputs

Follow the existing Python style and type annotations. Preserve output filenames
and column names used by downstream analysis. Integration tests replay captured
ReacNetGenerator output through a stub; see `tests/conftest.py` and
[reproducibility](docs/reproducibility.md).

Keep cluster accounts, partitions, module loads, and local paths in your copy
of `reactor/submit.example.sbatch`.
