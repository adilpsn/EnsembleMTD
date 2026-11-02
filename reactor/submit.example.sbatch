#!/bin/bash
# Example Slurm wrapper for a production ensemble.  Copy it, edit the two
# sections marked SITE, and submit the copy:
#
#   cp reactor/submit.example.sbatch submit.sbatch
#   $EDITOR submit.sbatch
#   sbatch submit.sbatch
#
# Nothing else in this repository contains scheduler directives or module
# commands, so this file is the only thing you should need to adapt.

#SBATCH --job-name=xtb_mtd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --output=mtd_%j.out
#SBATCH --error=mtd_%j.err

# --- SITE 1: partition, account, and anything else your scheduler requires ---
# Uncomment and edit. These vary between clusters and cannot be guessed.
#
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT

# --- SITE 2: make crest, xtb, packmol and Open Babel available ---------------
# Whatever your site needs. Some examples:
#
# module load GCC OpenMPI                    # Packmol often needs libgfortran
# module load xtb crest
# source /path/to/your/env.sh
# conda activate ensemblemtd
#
# The scripts check for their tools and exit with a clear message if any is
# missing, so a wrong guess here fails fast rather than half way through.

set -euo pipefail

# xtb is OpenMP-threaded and takes its thread count from the environment.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=4G
ulimit -s unlimited 2>/dev/null || true

repo=$(cd "$(dirname "$0")/.." && pwd)
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# --- the actual work ---------------------------------------------------------
# Build the cluster, then run the ensemble. Adjust the SMILES, the charge and
# the ensemble size; k/N = 0.5, alpha = 0.6 and density = 5 are the calibrated
# production settings (see docs/parameters.md).

"$repo/reactor/genmol.sh" "[Li+].C1COC(=O)O1.C1COC(=O)O1" liec2 --shape sphere

"$repo/reactor/setup_reactor.sh" liec2.xyz 0 0.5 0.6 5 \
    --time 30 \
    --mddump 100 \
    --nrun 20 \
    --walltemp 1000.0 \
    --nm liec2

# If the job is killed by the wall clock, resubmit unchanged: existing
# trajectories are counted first, so it resumes rather than starting over.
