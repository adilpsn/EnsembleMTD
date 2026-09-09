# Choosing the bias and the confinement

`setup_reactor.sh` takes five positional arguments:

```
setup_reactor.sh <xyz> <charge> <k/N> <alpha> <density> [options]
```

The production settings used throughout the paper are **k/N = 0.5, α = 0.6,
ρ = 5 g cm⁻³**, calibrated on EC ring-opening. They are a starting point, not a
universal recipe.

## The bias

The RMSD bias is a sum of Gaussians on the reference images collected so far,

$$E_\mathrm{bias}(t) = \sum_{i=1}^{n} k_i \exp\!\left[-\alpha_i \Delta_i^2(t)\right]$$

with $\Delta_i$ the RMSD from the *i*th reference structure. `crest --reactor`
picks a `kpush` from the system size; the script then rescales it:

```
kpush = kpush_crest * (k/N) * bias_ratio / 4
```

The `/4` maps crest's conformer-search default onto the reactive window used
here. `bias_ratio` is 1 unless `--nohbias` is given, in which case it is the
fraction of atoms actually biased. The value crest chose and the value handed
to xtb are both printed, so the arithmetic is visible in every log:

```
kpush: 0.440000 -> .055000 (factor=0.5, ratio=1.0)
```

**α (Gaussian width).** 0.6–0.7 gives reliable sampling, consistent with
Grimme's recommendation for RMSD metadynamics. EC ring-opening is also seen at
α = 0.9, but the product distribution narrows: a tighter Gaussian pushes the
system along reaction-specific distortions, which is the opposite of what a
discovery run is for. α = 0.6 is the production value.

**k/N (Gaussian height per atom).** 0.5 at α = 0.6 reproduced ring-opening
across the scanned grid. Too small and nothing happens within the trajectory
length; too large and the trajectory is driven rather than explored.

## The confinement, and why it is not optional

This is the part that catches people. RMSD is a high-dimensional collective
variable, so the cheapest way for a cluster to become dissimilar from every
reference structure is simply to **fly apart**. The bias gets absorbed by
intermolecular separation and no chemistry happens. A log-Fermi wall keeps the
fragments in contact so the bias has to go into intramolecular distortion
instead.

The diagnostic used for this is the RMSD ratio $R$: the intramolecular share of
the accumulated RMSD. $R \approx 1$ means the bias is deforming molecules;
small $R$ means it is just rearranging them. For LiEC₂:

| ρ (g cm⁻³) | Behaviour of $R$ |
|---|---|
| 2 | Flat at ≈ 0.15 for the whole trajectory — the bias goes almost entirely into intermolecular reorganisation |
| 4, 6 | Starts higher and rises past ≈ 15 ps to ≈ 0.18 — under tighter confinement more of the accumulated RMSD is forced into intramolecular distortion |

ρ = 5 was adopted as a compromise in that regime. `--walltemp` sets the
log-Fermi temperature parameter (1000 K in production); the remaining
confinement parameters are left at their xtb defaults.

Two practical notes:

- Clusters larger than about three carbonate molecules needed
  **system-specific** confinement tuning. That is an observation about these
  systems, not a size limit of the method — but expect to re-check $R$ rather
  than reuse ρ = 5 blindly.
- The density passed here is crest's `--genpot` parameter. crest reports both
  the unscaled and the requested density and the resulting cavity radius, so
  check the radius it derived is physically sensible for your cluster before
  launching twenty runs.

## Trajectory length and sampling

Production settings, for reference:

| Setting | Value | Flag |
|---|---|---|
| Trajectory length | 30 ps | `--time 30` |
| Time step | 0.5 fs | `--timestep 0.5` |
| Coordinate dump | every 50 fs (600 frames) | `--mddump 100` |
| Thermostat | Berendsen, 298.15 K | xtb default in the template |
| Hydrogen mass | 4 amu | crest reactor default |
| SHAKE | disabled for MTD, enabled for the relaxation | handled by the script |
| Reference images | one per ps, at most 20 retained | `save=20` in `rcontrol` |
| Wall temperature | 1000 K | `--walltemp 1000.0` |
| Ensemble size | 20 trajectories | `--nrun 20` |

`--mddump` is in steps, so with a 0.5 fs step `--mddump 100` is a frame every
50 fs.

**Where the ensemble spread comes from.** Every member starts from the same
packed geometry, then gets its own 0.5 ps unbiased room-temperature MD
relaxation with a fresh random seed. The differing velocities and relaxed
geometries are the entire source of the variation the analysis then measures.
There is no other randomisation.

## Controls worth running

- **Unbiased reference.** Same settings and confinement, `kpush = 0`. If a
  product appears in the unbiased ensemble too, the bias did not create it.
- **Neutral control.** The same cluster without the extra electron. For
  reduction-induced cleavage this is the control that matters: if the neutral
  ensemble shows the same chemistry at identical bias, the finding is a bias
  artefact rather than a reduction pathway.

## Restarting and provenance

Existing `<prefix>k<k/N>_a<alpha>run*.trj` files in the working directory are
counted before anything runs, so re-running with a larger `--nrun` tops the
ensemble up. Asking for a size you already have prints "Nothing to do" and
exits cleanly.

**The `rcontrol` left in the directory is a stale artefact of the last member** —
crest rewrites it for every run. The authoritative record of how a trajectory
was produced is line 2 of the `.trj` itself:

```
MTD params: charge=-1 kpush_factor=0.5, alp_value=0.6, Density=5, mddump=100, walltemp=1000.0, solvent=, bias_ratio=1.0
```

which is also echoed into `mtd_run<N>.out`. Audit an ensemble from those lines,
never from `rcontrol`.

## Other options

| Flag | Effect |
|---|---|
| `--nohbias` | Bias only non-hydrogen atoms, and scale `kpush` by the biased fraction |
| `--atoms "4-9"` | Bias an explicit atom range instead of everything |
| `--alpb <solvent>` | Add xtb's ALPB implicit solvent |
| `--nm <prefix>` | Prefix the output names, so several systems can share a directory |
| `--noclutter` | Keep only the MTD output and the final `.trj` |
| `--dryrun` | Generate `rcontrol` and stop — useful for checking the wall radius and `kpush` before committing node hours |

## Cost

For LiEC (11 atoms) on 36 OpenMP threads of a PALMA node: about 35 s per 30 ps
trajectory, ~710 s for a 20-run ensemble, and ~96 s for the ReacNetGenerator
and consensus analysis of those 20 runs. Machine-specific, but the ratio is the
useful part — the analysis is cheap next to the sampling.
