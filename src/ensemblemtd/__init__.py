"""EnsembleMTD -- ensemble RMSD-driven metadynamics for reaction discovery.

The reactor stage (``reactor/genmol.sh``, ``reactor/setup_reactor.sh``) produces
N independent biased GFN2-xTB trajectories from one starting geometry.  This
package is the analysis stage: it strips Li, runs ReacNetGenerator per
trajectory, collapses the resulting species onto fixed-H InChI groups, and
reports which species and which reaction steps recur across the ensemble.

The headline number is the species prevalence

    Pi_i = 100% * (runs in which collapsed species i appears) / (successful runs)

which measures reproducibility across independent trajectories, not abundance
or yield.  See :mod:`ensemblemtd.prevalence`.
"""

__version__ = "1.0.0"

__all__ = ["__version__"]
