# Reproducibility

ReacNetGenerator can assign different SMILES to the same fragment on repeated
analyses. Species are grouped by fixed-H InChI. Use `collapsed_inchi` to match
species between aggregations; representative SMILES and `species_id` can change.

In three repeat analyses of a six-run ensemble, all 12 collapsed species had
identical prevalence and run-support counts. Representative SMILES and
`member_raw_species_count` varied.

## Compare two analyses

```bash
ensemble-mtd-aggregate --inputs 'k*.trj' --outdir a --nohmm
ensemble-mtd-aggregate --inputs 'k*.trj' --outdir b --nohmm
```

Compare `pct_runs_observed` in the two species tables using `collapsed_inchi`
as the key. The summary JSON records the Open Babel version and analysis settings.

## Test fixtures

Integration tests use captured ReacNetGenerator output passed through a stub
binary via `--reacnet-bin`. To retain output from a new analysis:

```bash
ensemble-mtd-aggregate --inputs 'k*.trj' --outdir out --nohmm \
    --keep-temp --tmp-root scratch
```

Each run's files are saved under `scratch/rng_liagg_*/run_NNN/`.

Simulation parameters are recorded in line 2 of each `.trj` file and in the run
log. Each trajectory uses independently initialised velocities. Convergence
subsampling uses `--seed` (default `20240706`).

## Software used in the study

| Program | Version |
|---|---|
| xTB | 6.7.1 |
| CREST | 3.0.2 |
| ReacNetGenerator | 1.6.16 |
| Packmol | 20.16.0 |

The simulations use GFN2-xTB.
