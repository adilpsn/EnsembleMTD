# Reproducibility

## ReacNetGenerator is not bit-reproducible

Running the pipeline twice on the *same* trajectories does not necessarily give
the same species labels. ReacNetGenerator can write a different SMILES for the
same fragment between runs — most often choosing between a double bond and two
radical centres, e.g.

```
[H][C]([H])=[C]([H])[O][C]([O])=[O]      and      [H][C]([H])[C]([H])[O][C]([O])=[O]
[O]=[C]=[O]                              and      [O][C][O]
```

Each pair is one molecule with one fixed-H InChI. Which string comes back
varies; so, therefore, does which one gets promoted to represent its group.

This is worth stating plainly because it looks alarming and mostly is not. The
fixed-H InChI collapse exists partly to absorb exactly this.

## What is stable and what is not

Repeating a six-run aggregation three times (twice with the original
implementation, once with this one), keyed on `collapsed_inchi`:

| Quantity | Stable? |
|---|---|
| `pct_runs_observed` (Π_i) | Yes — identical in all three |
| `observed_run_support` | Yes — identical in all three |
| `is_initial` | Yes |
| Set of collapsed species | Yes — 12 species every time |
| `member_raw_species_count` | **No** — e.g. 2 / 3 / 1 for the vinyloxy radical |
| Representative SMILES, and hence `species_id` | **No** — follows the member set |

So: **every prevalence and every run-support count was identical**, while the
label attached to a group moved. That is the collapse doing its job.

Two consequences for how you use the output:

1. **Key on `collapsed_inchi`, not on the SMILES or on `S#`.** Species indices
   are assigned per aggregation and are not comparable between aggregations,
   even of the same data. The manuscript says so explicitly for its own
   figures.
2. **Do not read anything into `member_raw_species_count`.** It counts how many
   spellings ReacNetGenerator happened to emit, not chemistry.

## Testing it on your own data

Run the aggregation twice into different directories and compare on the InChI:

```bash
ensemble-mtd-aggregate --inputs 'k*.trj' --outdir a --nohmm
ensemble-mtd-aggregate --inputs 'k*.trj' --outdir b --nohmm

python - <<'PY'
import csv
def load(d):
    rows = csv.DictReader(open(f"{d}/aggregate_reacnet_species.tsv"), delimiter="\t")
    return {r["collapsed_inchi"]: r["pct_runs_observed"] for r in rows}
a, b = load("a"), load("b")
assert a == b, {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)}
print(f"{len(a)} collapsed species, prevalences identical")
PY
```

## Making the pipeline deterministic

Everything downstream of ReacNetGenerator *is* deterministic: same inputs, same
outputs, byte for byte. To test that half in isolation, replay captured
ReacNet output through a stub binary via `--reacnet-bin`. This is how the test
suite works — see `tests/conftest.py` for the stub, which is about ten lines.
It is also how this implementation was checked against the original single-file
script it replaces: with the stub in place, all eleven output artefacts came out
byte-identical across twelve different option combinations, plus the excluded-run,
failed-run, rebuild and dry-run paths.

To capture your own fixtures:

```bash
mkdir scratch
ensemble-mtd-aggregate --inputs 'k*.trj' --outdir out --nohmm \
    --keep-temp --tmp-root scratch
# scratch/rng_liagg_*/run_NNN/ now holds the .json, .species and .html per run
```

## Everything else that is fixed

- **The MTD ensemble.** Each member's velocities come from its own random seed,
  so individual trajectories are not reproducible by design — that is what the
  ensemble measures. The settings that produced each one are recorded in line 2
  of its `.trj`.
- **The convergence analysis** uses a seeded generator
  (`--seed`, default 20240706), so its subsampling band is reproducible.
- **Ordering.** Every table has an explicit sort with a string tie-break, and
  the representative of a collapsed group is chosen by run support with the
  SMILES as tie-break. Nothing depends on Python set or dict iteration order.
- **Open Babel version** is recorded in the summary, because the InChI grouping
  depends on the build.

## Versions used for the published results

| Component | Version |
|---|---|
| `xtb` | 6.7.1 |
| `crest` | 3.0.2 |
| ReacNetGenerator | 1.6.16 |
| Packmol | 20.16.0 |
| Hamiltonian | GFN2-xTB |
