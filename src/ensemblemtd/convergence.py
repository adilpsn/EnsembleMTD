"""How many trajectories does an ensemble need?

Reads an aggregate directory and traces the prevalence of chosen species as
trajectories are added, Pi_i(N), from the per-run membership table.  Two curves
per species:

*native*
    Runs in the order they were produced (run1, run2, ...).  This is the curve
    an author would actually have watched while the ensemble was accumulating.

*subsampling band*
    The 5th-95th percentile of Pi_i(N) over many random run orderings, which
    removes any accident of that particular order.

The band's *mean* is not a convergence test: because the runs are exchangeable,
the mean equals the full-ensemble Pi_i at every N by construction, and the band
necessarily closes to zero at N = N_total because sampling is without
replacement.  What the band does say is how much a prevalence could have moved
had the runs arrived in a different order -- so the useful readouts are its
half-width at the production ensemble size, and the smallest N beyond which the
native curve stays inside a stated tolerance.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

DEFAULT_ORDERINGS = 2000
DEFAULT_SEED = 20240706
DEFAULT_TOLERANCE_PP = 5.0
DEFAULT_TOP = 5

BAND_LOW_PERCENTILE = 5
BAND_HIGH_PERCENTILE = 95

SPECIES_TSV = "aggregate_reacnet_species.tsv"
MEMBERSHIP_TSV = "aggregate_reacnet_species_run_membership.tsv"

CSV_COLUMNS = (
    "species_id",
    "smiles",
    "collapsed_inchi",
    "is_initial",
    "pi_full",
    "N",
    "pi_native",
    "pi_band_mean",
    "pi_band_lo",
    "pi_band_hi",
)

RUN_NUMBER_RE = re.compile(r"run(\d+)")


class ConvergenceError(RuntimeError):
    """The aggregate directory does not contain what is needed."""


@dataclass(frozen=True)
class Species:
    species_id: str
    smiles: str
    collapsed_inchi: str
    is_initial: bool
    run_support: int
    prevalence: float
    runs: frozenset


@dataclass(frozen=True)
class Curve:
    species: Species
    native: Tuple[float, ...]
    band_mean: Tuple[float, ...]
    band_lo: Tuple[float, ...]
    band_hi: Tuple[float, ...]

    def half_width_at(self, n: int) -> Optional[float]:
        if not 1 <= n <= len(self.band_lo):
            return None
        return (self.band_hi[n - 1] - self.band_lo[n - 1]) / 2.0

    def settles_at(self, tolerance_pp: float) -> int:
        """Smallest N beyond which the native curve stays within tolerance."""
        target = self.species.prevalence
        total = len(self.native)
        for start in range(total):
            if all(abs(self.native[k] - target) <= tolerance_pp for k in range(start, total)):
                return start + 1
        return total


def run_number(run_name: str) -> int:
    match = RUN_NUMBER_RE.search(run_name)
    if not match:
        raise ConvergenceError(f"Cannot read a run number from {run_name!r}")
    return int(match.group(1))


def read_ensemble(aggregate_dir: Path) -> Tuple[List[Species], List[int]]:
    """Load the species table and the per-run membership matrix."""
    species_path = aggregate_dir / SPECIES_TSV
    membership_path = aggregate_dir / MEMBERSHIP_TSV
    for path in (species_path, membership_path):
        if not path.exists():
            raise ConvergenceError(f"Missing {path.name} in {aggregate_dir}")

    membership: Dict[str, Set[int]] = {}
    all_runs: Set[int] = set()
    with membership_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            field = row["runs"].strip()
            runs = {run_number(r) for r in field.split(";") if r} if field else set()
            membership[row["species_id"]] = runs
            all_runs |= runs

    species: List[Species] = []
    with species_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            species.append(
                Species(
                    species_id=row["species_id"],
                    smiles=row["species"],
                    collapsed_inchi=row.get("collapsed_inchi", ""),
                    is_initial=bool(int(row["is_initial"])),
                    run_support=int(row["observed_run_support"]),
                    prevalence=float(row["pct_runs_observed"]),
                    runs=frozenset(membership.get(row["species_id"], set())),
                )
            )

    if not all_runs:
        raise ConvergenceError(f"No runs recorded in {membership_path.name}")
    return species, sorted(all_runs)


def select_species(
    species: Sequence[Species],
    smiles: Optional[Sequence[str]] = None,
    top: int = DEFAULT_TOP,
    include_initial: bool = False,
) -> List[Species]:
    """Pick the species to trace: an explicit SMILES list, or the top by support.

    Reactants are left out of the automatic selection by default.  A reactant is
    present in the first frame of every run by construction, so its curve is a
    flat line at 100% and says nothing about convergence.
    """
    if smiles:
        by_smiles = {s.smiles: s for s in species}
        missing = [x for x in smiles if x not in by_smiles]
        if missing:
            raise ConvergenceError(
                "Not in the aggregate: " + ", ".join(repr(m) for m in missing)
            )
        return [by_smiles[x] for x in smiles]

    pool = [s for s in species if include_initial or not s.is_initial]
    ranked = sorted(pool, key=lambda s: (-s.run_support, s.species_id))
    return ranked[:top]


def native_curve(present: Sequence[float]) -> Tuple[float, ...]:
    """Running prevalence over the first N runs, in the given order."""
    running = 0.0
    out = []
    for index, value in enumerate(present, start=1):
        running += value
        out.append(100.0 * running / index)
    return tuple(out)


def build_curve(
    species: Species,
    runs: Sequence[int],
    orderings,
    numpy_module,
) -> Curve:
    np = numpy_module
    present = np.array([1.0 if r in species.runs else 0.0 for r in runs])
    cumulative = np.cumsum(present[orderings], axis=1) / np.arange(1, len(runs) + 1)
    band = 100.0 * cumulative
    return Curve(
        species=species,
        native=native_curve(present.tolist()),
        band_mean=tuple(band.mean(axis=0)),
        band_lo=tuple(np.percentile(band, BAND_LOW_PERCENTILE, axis=0)),
        band_hi=tuple(np.percentile(band, BAND_HIGH_PERCENTILE, axis=0)),
    )


def build_curves(
    species: Sequence[Species],
    runs: Sequence[int],
    n_orderings: int = DEFAULT_ORDERINGS,
    seed: int = DEFAULT_SEED,
) -> List[Curve]:
    """Trace every selected species over a shared set of random orderings."""
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - numpy is an extra
        raise ConvergenceError(
            "numpy is required for the convergence analysis: "
            "pip install 'ensemblemtd[plots]'"
        ) from exc

    rng = np.random.default_rng(seed)
    # One shared set of orderings, so the curves are directly comparable.
    orderings = np.array([rng.permutation(len(runs)) for _ in range(n_orderings)])
    return [build_curve(s, runs, orderings, np) for s in species]


def write_csv(path: Path, curves: Sequence[Curve]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)
        for curve in curves:
            s = curve.species
            for index in range(len(curve.native)):
                writer.writerow(
                    [
                        s.species_id,
                        s.smiles,
                        s.collapsed_inchi,
                        int(s.is_initial),
                        f"{s.prevalence:.3f}",
                        index + 1,
                        f"{curve.native[index]:.3f}",
                        f"{curve.band_mean[index]:.3f}",
                        f"{curve.band_lo[index]:.3f}",
                        f"{curve.band_hi[index]:.3f}",
                    ]
                )


def write_plot(path: Path, curves: Sequence[Curve], production_size: Optional[int]) -> None:
    """Draw the curves.  Colour-blind-safe palette (Okabe-Ito)."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - matplotlib is an extra
        raise ConvergenceError(
            "matplotlib is required for --plot: pip install 'ensemblemtd[plots]'"
        ) from exc

    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    fig, ax = plt.subplots(figsize=(6.8, 3.4), constrained_layout=True)

    for index, curve in enumerate(curves):
        color = palette[index % len(palette)]
        n_axis = range(1, len(curve.native) + 1)
        ax.fill_between(n_axis, curve.band_lo, curve.band_hi, color=color, alpha=0.14, linewidth=0)
        ax.plot(n_axis, curve.native, color=color, label=curve.species.species_id)
        ax.axhline(curve.species.prevalence, color=color, linestyle=":", linewidth=0.8, alpha=0.85)

    if production_size:
        ax.axvline(production_size, color="#4B5563", linestyle=(0, (2, 2)), linewidth=0.8)
        ax.text(
            production_size + 0.4,
            4.0,
            f"N = {production_size}",
            ha="left",
            va="bottom",
            color="#4B5563",
            fontsize=8,
        )

    total = len(curves[0].native) if curves else 1
    ax.set_xlim(1, total)
    ax.set_ylim(0, 102)
    ax.set_xlabel("Number of trajectories included, N")
    ax.set_ylabel(r"Species prevalence, $\Pi_i$ (%)")
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)
    ax.legend(title="Species", frameon=False, fontsize=8, title_fontsize=8, ncol=2)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ensemble-mtd-convergence",
        description=(
            "Trace species prevalence against ensemble size using the per-run "
            "membership table written by ensemble-mtd-aggregate."
        ),
    )
    p.add_argument("aggregate_dir", help="An aggregate output directory")
    p.add_argument("--out-csv", default="convergence.csv", help="Where to write the curves")
    p.add_argument("--plot", default=None, help="Also write a PDF figure here")
    p.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="Representative SMILES to trace, in plot order. Defaults to the "
        "most-supported products.",
    )
    p.add_argument(
        "--top", type=int, default=DEFAULT_TOP, help="How many species to trace"
    )
    p.add_argument(
        "--include-initial",
        action="store_true",
        help="Also consider reactants, whose curves are flat by construction",
    )
    p.add_argument("--orderings", type=int, default=DEFAULT_ORDERINGS)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE_PP,
        help="Tolerance in percentage points for the settle-N diagnostic",
    )
    p.add_argument(
        "--production-size",
        type=int,
        default=None,
        help="Ensemble size to mark and report the band half-width at",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        species, runs = read_ensemble(Path(args.aggregate_dir))
        selected = select_species(
            species,
            smiles=args.species,
            top=args.top,
            include_initial=args.include_initial,
        )
        if not selected:
            raise ConvergenceError("No species selected")
        curves = build_curves(selected, runs, args.orderings, args.seed)
    except ConvergenceError as exc:
        print(str(exc))
        return 2

    production = args.production_size or len(runs)
    print(f"{len(runs)} runs, {len(curves)} species traced")
    for curve in curves:
        half_width = curve.half_width_at(production)
        settle = curve.settles_at(args.tolerance)
        band = "n/a" if half_width is None else f"{half_width:.1f} pp"
        print(
            f"  {curve.species.species_id:>4s}  Pi = {curve.species.prevalence:5.1f}%"
            f"  band half-width at N={production}: {band}"
            f"  within +-{args.tolerance:.0f} pp from N = {settle}"
            f"  {curve.species.smiles}"
        )

    write_csv(Path(args.out_csv), curves)
    print(f"Wrote {args.out_csv}")
    if args.plot:
        write_plot(Path(args.plot), curves, args.production_size)
        print(f"Wrote {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
