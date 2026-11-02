# Simulation parameters

```text
setup_reactor.sh <xyz> <charge> <k/N> <alpha> <density> [options]
```

The example uses `k/N = 0.5`, `alpha = 0.6`, and a confinement density parameter
of `5 g cm⁻³`. See the [preprint](https://doi.org/10.26434/chemrxiv.15011301/v1)
for the parameter study.

## Bias and confinement

CREST generates the initial `rcontrol`. The script rescales its Gaussian height:

```text
kpush = kpush_crest * (k/N) * bias_ratio / 4
```

`bias_ratio` is the fraction of biased atoms with `--nohbias`, and 1 otherwise.
The original and rescaled values are printed in the log.

A log-Fermi wall confines the cluster. The density argument is passed to
CREST's `--genpot`; check the generated cavity radius with `--dryrun` when
setting up a new system. `--walltemp` sets the wall temperature parameter.

## Example settings

| Setting | Value | Option |
|---|---|---|
| Trajectory length | 30 ps | `--time 30` |
| Time step | 0.5 fs | `--timestep 0.5` |
| Coordinate output | Every 100 steps (50 fs) | `--mddump 100` |
| Thermostat | Berendsen, 298.15 K | Template default |
| Hydrogen mass | 4 amu | CREST reactor default |
| SHAKE | Enabled for relaxation, disabled for MTD | Set by script |
| Reference images | One per ps, up to 20 retained | `save=20` in `rcontrol` |
| Wall temperature | 1000 K | `--walltemp 1000.0` |
| Ensemble size | 20 | `--nrun 20` |

Each run begins with a 0.5 ps unbiased relaxation from the packed geometry,
using independently initialised velocities.

## Options

| Option | Effect |
|---|---|
| `--nohbias` | Bias non-hydrogen atoms and rescale `kpush` |
| `--atoms "4-9"` | Bias a selected atom range |
| `--alpb <solvent>` | Use xTB ALPB implicit solvent |
| `--nm <prefix>` | Set the output prefix |
| `--noclutter` | Keep the MTD output and final trajectory |
| `--dryrun` | Generate `rcontrol` without running trajectories |

## Restarting and run records

Existing trajectory files are counted before the ensemble starts. Increasing
`--nrun` adds members. Parameters for each run are recorded in line 2 of its
`.trj` and in `mtd_run<N>.out`. CREST rewrites `rcontrol` for each member, so use
the trajectory and log to recover settings for an earlier run.
