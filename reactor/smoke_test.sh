#!/bin/bash
# Exercise the reactor scripts end to end on a Li + EC cluster small enough to
# finish in under a minute.  This is the shell counterpart to the pytest suite,
# which covers the analysis stage only.
#
# Run it on a compute node, not a login node:
#
#   sbatch reactor/smoke_test.sh
#
# or directly, if crest and xtb are on PATH and you are somewhere you are
# allowed to burn a couple of CPU minutes:
#
#   bash reactor/smoke_test.sh
#
#SBATCH -J emtd_smoke
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH -o smoke_test.out
#SBATCH -e smoke_test.err

# Cluster-specific: Packmol needs libgfortran, which a bare batch environment
# does not always carry.  Adjust or drop for your site.
if command -v module >/dev/null 2>&1; then
    module load palma/2023b GCC/13.2.0 OpenMPI/4.1.6 2>/dev/null
fi

# Locating the sibling scripts is not as simple as dirname "$0": slurmd copies
# the batch script into its spool directory before running it, so under sbatch
# "$0" points somewhere that contains nothing else.  Prefer the submit
# directory, and let REACTOR_DIR override either way.
for candidate in \
        "${REACTOR_DIR:-}" \
        "${SLURM_SUBMIT_DIR:-}/reactor" \
        "${SLURM_SUBMIT_DIR:-}" \
        "$(cd "$(dirname "$0")" 2>/dev/null && pwd)"; do
    if [ -n "$candidate" ] && [ -f "$candidate/setup_reactor.sh" ]; then
        here=$candidate
        break
    fi
done
if [ -z "${here:-}" ]; then
    echo "Cannot find setup_reactor.sh next to this script." >&2
    echo "Set REACTOR_DIR to the directory holding it and re-run." >&2
    exit 1
fi

for tool in obabel packmol crest xtb; do
    command -v "$tool" >/dev/null || {
        echo "Required tool not on PATH: $tool" >&2
        echo "Load your modules or source your environment before running." >&2
        exit 1
    }
done

work=${SLURM_SUBMIT_DIR:-$PWD}/smoke_test_work
rm -rf "$work" && mkdir -p "$work" && cd "$work" || exit 1

failed=0
say() { echo; echo "########## $* ##########"; }
check() { if [ "$1" -eq 0 ]; then echo "  ok: $2"; else echo "  FAILED: $2"; failed=1; fi }

say "tools"
echo "  scripts   $here"
for tool in obabel packmol crest xtb; do
    printf '  %-9s %s\n' "$tool" "$(command -v $tool)"
done

say "1. genmol.sh builds a cluster"
bash "$here/genmol.sh" "[Li+].C1COC(=O)O1" liec --shape sphere
check $? "genmol.sh wrote liec.xyz"
[ "$(head -1 liec.xyz | tr -d '[:space:]')" = "11" ]
check $? "11 atoms as expected"

say "2. --dryrun generates rcontrol without running any MD"
bash "$here/setup_reactor.sh" liec.xyz -1 0.5 0.6 5 --time 20 --dryrun
check $? "dry run completed"
grep -q '^\$metadyn' rcontrol && grep -q '^\$wall' rcontrol
check $? "rcontrol has the metadyn and wall blocks"
grep -q 'temp=1000' rcontrol
check $? "log-Fermi wall temperature applied"
# crest picks kpush from system size; the driver rescales it by factor/4
grep -qE 'kpush=0?\.0*55000' rcontrol
check $? "kpush rescaled to 0.44 * 0.5 / 4 = 0.055"
[ ! -f xtb.trj ]
check $? "no trajectory produced by a dry run"

say "3. a two-member ensemble runs"
bash "$here/setup_reactor.sh" liec.xyz -1 0.5 0.6 5 \
    --time 1 --mddump 100 --nrun 2 --walltemp 1000.0 --nm test
check $? "ensemble completed"
[ "$(ls testk0.5_a0.6run*.trj 2>/dev/null | wc -l)" -eq 2 ]
check $? "two trajectories written"
sed -n 2p testk0.5_a0.6run1.trj | grep -q 'MTD params: charge=-1'
check $? "parameters recorded in line 2 of the trajectory"

say "4. a larger --nrun tops the ensemble up"
# Log to a file rather than piping into grep -q -- grep exits on its first
# match, and the SIGPIPE that follows would kill the run part way through.
bash "$here/setup_reactor.sh" liec.xyz -1 0.5 0.6 5 \
    --time 1 --mddump 100 --nrun 3 --walltemp 1000.0 --nm test > resume.log 2>&1
check $? "resumed run completed"
grep -qi 'indices from 3' resume.log
check $? "started at run 3 rather than run 1"
[ "$(ls testk0.5_a0.6run*.trj | wc -l)" -eq 3 ]
check $? "three trajectories now present"

say "5. a full ensemble is left alone"
bash "$here/setup_reactor.sh" liec.xyz -1 0.5 0.6 5 --time 1 --nrun 3 --nm test \
    > full.log 2>&1
grep -qi 'nothing to do' full.log
check $? "declined to re-run a complete ensemble"
[ "$(ls testk0.5_a0.6run*.trj | wc -l)" -eq 3 ]
check $? "still three trajectories, nothing overwritten"

say "6. bad input is refused"
bash "$here/setup_reactor.sh" missing.xyz -1 0.5 0.6 5 >/dev/null 2>&1
[ $? -ne 0 ]
check $? "missing geometry exits non-zero"

say "result"
if [ $failed -eq 0 ]; then
    echo "  all checks passed; working files are in $work"
else
    echo "  something failed; see above and the files in $work"
fi
exit $failed
