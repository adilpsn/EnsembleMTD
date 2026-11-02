# Output files

`ensemble-mtd-aggregate --outdir rcng` writes the following files into `rcng/`.

## Main files

### `aggregate_reacnet_species.tsv`

One row per collapsed species, reactants first (in the order they were first
seen), then everything else by prevalence.

| Column | Meaning |
|---|---|
| `species_id` | `S1`, `S2`, … Stable within one aggregation, assigned per aggregation |
| `species` | The representative SMILES for the group |
| `collapsed_inchi` | Fixed-H InChI used to identify the group |
| `representative_smiles` | Same as `species`; kept for older readers |
| `member_raw_species_count` | How many distinct raw ReacNet SMILES collapsed into this group |
| `display_run_support` | `max(observed_run_support, first_frame_run_support)` — for display only |
| `observed_run_support` | Number of successful runs containing the species |
| `observed_run_support_json` | The same count from ReacNet's JSON species list, which `--maxspecies` truncates |
| `support_delta_json` | The gap between the two. Non-zero means the JSON list undercounted |
| `support_source` | Which source `observed_run_support` came from (`timeline` by default) |
| `pct_runs_observed` | Percentage of successful runs containing the species |
| `is_initial` | 1 if present in the first analysed frame |
| `first_frame_run_support` | Runs in which it was in the first frame |
| `initial_frame_count` | Summed first-frame multiplicity |
| `involvement_events` | Total reaction events this species took part in, either side |

### `aggregate_reacnet_species_run_membership.tsv`

The presence matrix: for each species, the **names of the runs** it appeared
in, semicolon-separated. Everything else stores only counts, so this is the
file needed for any per-run or resampling analysis — including
`ensemble-mtd-convergence`, which reads it. The run sets here are the same ones
`pct_runs_observed` is computed from, so the counts reproduce exactly.

### `aggregate_reacnet_reactions.tsv`

Collapsed reaction records, ordered by run support then total events.

| Column | Meaning |
|---|---|
| `reaction` | `lhs => rhs`, with stoichiometry written as `2*X` |
| `run_support` | Distinct runs containing this record |
| `pct_runs` | The same as a percentage |
| `total_events` | Summed event count across the ensemble |

These are the **unnetted, directional** weights. `--clean` affects only the
drawn network, never this table.

### `aggregate_reacnet_network.dot`

The consensus network as Graphviz source. Nodes are `S#` (map them through
`aggregate_reacnet_network_nodes.tsv` or the species table); reactants are
outlined in indigo, a highlighted route in red. Every edge is labelled
`e=<event weight>, r=<run support>` and its tooltip carries the full SMILES.
Render with:

```bash
dot -Tpdf rcng/aggregate_reacnet_network.dot -o network.pdf
```

`aggregate_reacnet_network.svg` is also written when Graphviz is installed.
Use `--rebuild-from` to regenerate the network from an existing analysis.

### `aggregate_reacnet_summary.json`

Provenance and run accounting. Includes:

- `n_input_runs`, `n_success_runs`, `n_excluded_runs`, `n_failed_runs` — run
  counts used to calculate prevalence.
- `excluded_runs` — runs dropped because a frame reported zero energy, i.e. a
  collapsed SCF.
- `run_stats` — per run: frames, atoms before and after Li removal, species and
  reaction counts, and the exclusion reason if any.
- `n_raw_species`, `n_collapsed_species`, `collapse_ratio` — how much the InChI
  grouping merged.
- `obabel_version` — the InChI grouping depends on the Open Babel build, so
  this is part of the provenance of every species label.
- `reacnet_settings` — what ReacNetGenerator was actually told, `nohmm`
  included.
- `graph` — the thresholds applied, the edge and node counts, and the
  highlighted route if one was drawn.
- `warnings`, `failures` — same content as the two TSVs below.

## Supporting files

| File | Contents |
|---|---|
| `aggregate_reacnet_inchi_alias_map.tsv` | Mapping from raw ReacNet SMILES to collapsed groups |
| `aggregate_reacnet_network_nodes.tsv` | `S#` → SMILES for the nodes actually drawn, with their support |
| `aggregate_reacnet_mechanistic_path.tsv` | The highlighted route, one row per step. Header only when `--pathway-mode none` |
| `aggregate_reacnet_pathlength_scan.tsv` | Written only with `--scan-pathlength`: the best route of each exact length |
| `aggregate_reacnet_report.html` | Standalone page: run list, top reactions, top species. No dependencies |
| `aggregate_reacnet_reacnetstyle.html` | ReacNetGenerator's own interactive shell, with the aggregate substituted in and run-support labels added to each row |
| `aggregate_reacnet_rngdata.json` | The payload behind that page |
| `aggregate_reacnet_failures.tsv` | Per-run errors. Written only if something failed |
| `aggregate_reacnet_warnings.tsv` | Excluded runs and JSON-vs-timeline support gaps. Written only if there are any |

## Exit codes

| Code | Meaning |
|---|---|
| 0 | At least one run was aggregated |
| 2 | Bad arguments, no matching inputs, or the InChI collapse failed |
| 3 | Ran, but no run succeeded (or, for `--rebuild-from`, something in the rebuild failed) |

Failed trajectories are recorded in `aggregate_reacnet_failures.tsv`.
The remaining successful runs are aggregated.

## Intermediate files

Per-run ReacNetGenerator output lives in a temporary directory that is removed
on exit. `--keep-temp` (with `--tmp-root <dir>` to control where) keeps it, one
subdirectory per run, which is the way to inspect a per-run `.species` timeline
or a ReacNet log. `--save-noli-dir <dir>` keeps just the Li-removed
trajectories.
