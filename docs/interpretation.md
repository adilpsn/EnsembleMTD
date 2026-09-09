# What the numbers mean, and what they do not

Three quantities come out of `ensemble-mtd-aggregate`. They answer different
questions and it is easy to read more into them than they carry.

## Π_i — species prevalence (a node property)

$$\Pi_i = 100\,\% \times \frac{N_i^\mathrm{runs}}{N_\mathrm{total}^\mathrm{runs}}$$

where $N_i^\mathrm{runs}$ counts the successful runs in which collapsed species
$i$ appeared in **at least one** analysed frame.

Π_i **is** the fraction of independent trajectories that produced a species at
all. It is the right quantity for "is this recurrent, or did one trajectory
happen to go there".

Π_i is **not**:

- a **concentration or abundance**. The indicator is binary per run, so a
  species present in every frame of every run and one present in a single frame
  of every run both give Π_i = 100 %.
- a **residence time or lifetime**. Temporal filtering is off; one frame is
  enough. The membership table records which runs, not for how long.
- a **yield**. Nothing is normalised by the amount of reactant consumed.
- a **branching ratio**. Two competing products can both reach Π_i = 100 % if
  every trajectory makes some of each, and the numbers need not sum to anything.
- a **rate or a barrier**. The trajectories are biased; the bias is active from
  t = 0. Nothing here is a kinetic observable.

Because the reference is the number of *successful* runs, excluded and failed
runs change the denominator. `aggregate_reacnet_summary.json` reports
`n_input_runs`, `n_success_runs`, `n_excluded_runs` and `n_failed_runs` so the
denominator is always visible.

### The reactant caveat

A species present in the first analysed frame is flagged `is_initial`. Under a
strong bias a reactant can be gone before the first dump is written, which
would make it look rare, so the *displayed* support is
`max(observed_run_support, first_frame_run_support)`. Both columns are in the
table; `observed_run_support` and `pct_runs_observed` are the unadjusted
numbers, and those are the ones to quote.

## e_uv — edge event weight (an edge property)

A ReacNet record can carry several species on each side, and the consensus
figure is a graph over single species, so each record is split over the ordered
pairs it implies:

$$e_{uv} = \sum_{q:\,u\in L_q,\,v\in R_q,\,u\neq v} \frac{n_q}{|L_q|\,|R_q|}$$

The denominator counts unique species on both sides **before** self-pairs are
dropped. So a record `A + B -> A + C` has denominator 4, the `A -> A` pair is
omitted, and the three surviving edges carry 3/4 of the record's events between
them. Consequently e_uv is generally **not an integer**, and summed edge weight
is not conserved. It is a relative measure of how much reactive traffic a step
carried, nothing more.

## r_uv — edge run support (an edge property)

The number of distinct successful runs contributing at least one record to
`u -> v`. This is **not** the reaction-count matrix element $r_{ij}$ of the
original ReacNetGenerator formulation, despite the letter.

A node's Π_i and an edge's r_uv are independent: a species can be present in
every run (Π_i = 100 %) while no single formation edge is supported by more
than a few, because different runs made it different ways.

## What the collapse costs

Species are grouped by fixed-H InChI after Li is removed. Two things follow:

1. The labels describe the **Li-stripped organic framework**, not charge-resolved
   speciation. A node is a heavy-atom skeleton with a given hydrogen inventory.
   It is not a charge state, and `[F]` in an anion study may in fact be LiF.
2. Li-mediated effects are **not resolved**. Removing Li is what stops Li⁺
   hopping between coordination sites from registering as a reaction, and the
   price is that Li coordination states and any Li catalysis are invisible. Any
   claim about them needs a different analysis.

## Cleaned networks

With `--clean`, an opposing edge pair is drawn as a single edge with
$e^\mathrm{net} = |e_{uv} - e_{vu}|$ in the direction of the larger weight and
$r = \max(r_{uv}, r_{vu})$. An exact tie is kept as a dashed bidirectional edge
with $e = 0$ — that is a statement about the ensemble (both directions were seen
equally often), not a rendering artefact. The unnetted directional weights are
always in `aggregate_reacnet_reactions.tsv`.

## Highlighted paths

`--pathway-mode` draws one route through the network. It is a reading aid so a
figure has something to follow. Neither rule (`support-first`,
`probabilistic`) produces a kinetic path, a mechanism, or a barrier; the
evidence is the edge and node statistics, and the route only picks a way through
them. The default is `none`.

## How many runs

Prevalence estimated from N trajectories carries the binomial uncertainty of N
trajectories, and 20 is not a large N. For the one system where this was tested
at length (2LiEC, extended to 40 runs), the dominant products were stable by
N ≈ 10–15 while species around 20–35 % were still moving by 7–12 percentage
points at N = 20. Twenty runs identify recurrent chemistry; they do not pin
prevalence values or rank minor species. `ensemble-mtd-convergence` measures
this for your own ensemble — see the caveat there about why the subsampling
band's *mean* is not a convergence test.
