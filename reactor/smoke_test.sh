#!/bin/bash
# Exercise the reactor scripts end to end on a Li + EC cluster small enough to
# finish in about a minute.  This is the shell counterpart to the pytest suite,
# which covers the analysis stage only.
#
# It needs crest and xtb, so run it where you are allowed a couple of CPU
# minutes -- a compute node or an interactive allocation, not a shared login
# node:
#
#   bash reactor/smoke_test.sh
#   srun --ntasks=8 --time=00:20:00 bash reactor/smoke_test.sh
#
# To submit it instead, use submit.example.sbatch as the wrapper.  There are no
# scheduler directives here, so nothing site-specific to edit in this file.

# Locating the sibling scripts is not as simple as dirname "$0": a scheduler
# may copy the script into a spool directory before running it, so under sbatch
# "$0" resolves somewhere that holds nothing else.  Prefer the submit
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
pass() { echo "  ok: $1"; }
fail() { echo "  FAILED: $1"; failed=1; }

# Run a command, showing its output, and report on its exit status.
step() {
    local desc=$1
    shift
    if "$@"; then pass "$desc"; else fail "$desc"; fi
}

# Assert a condition quietly.  Passing the command as arguments, rather than
# testing and then reading $?, keeps it unambiguous which status is being read.
expect() {
    local desc=$1
    shift
    if "$@" >/dev/null 2>&1; then pass "$desc"; else fail "$desc"; fi
}

# Count trajectories matching a prefix without parsing ls output.
count_trj() {
    local n=0 f
    for f in "$1"*.trj; do
        [ -e "$f" ] && n=$((n + 1))
    done
    echo "$n"
}

atom_count() { head -n 1 "$1" | tr -d '[:space:]'; }

say "tools"
echo "  scripts   $here"
for tool in obabel packmol crest xtb; do
    printf '  %-9s %s\n' "$tool" "$(command -v "$tool")"
done

say "1. genmol.sh builds a cluster"
step "genmol.sh wrote liec.xyz" \
    bash "$here/genmol.sh" "[Li+].C1COC(=O)O1" liec --shape sphere
expect "11 atoms as expected" test "$(atom_count liec.xyz)" = "11"

say "2. --dryrun generates rcontrol without running any MD"
step "dry run completed" \
    bash "$here/setup_reactor.sh" liec.xyz -1 0.5 0.6 5 --time 20 --dryrun
expect "rcontrol has a metadyn block" grep -q '^[$]metadyn' rcontrol
expect "rcontrol has a wall block" grep -q '^[$]wall' rcontrol
expect "log-Fermi wall temperature applied" grep -q 'temp=1000' rcontrol
# crest picks kpush from the system size; the driver rescales it by factor/4
expect "kpush rescaled to 0.44 * 0.5 / 4 = 0.055" \
    grep -qE 'kpush=0?\.0*55000' rcontrol
expect "no trajectory produced by a dry run" test ! -f xtb.trj

say "3. a two-member ensemble runs"
step "ensemble completed" \
    bash "$here/setup_reactor.sh" liec.xyz -1 0.5 0.6 5 \
        --time 1 --mddump 100 --nrun 2 --walltemp 1000.0 --nm test
expect "two trajectories written" test "$(count_trj testk0.5_a0.6run)" = "2"
sed -n 2p testk0.5_a0.6run1.trj > line2.txt
expect "parameters recorded in line 2 of the trajectory" \
    grep -q 'MTD params: charge=-1' line2.txt

say "4. a larger --nrun tops the ensemble up"
# Redirect to a file rather than piping into grep -q -- grep exits on its first
# match, and the SIGPIPE that follows would kill the run part way through.
step "resumed run completed" \
    bash -c "bash '$here/setup_reactor.sh' liec.xyz -1 0.5 0.6 5 \
        --time 1 --mddump 100 --nrun 3 --walltemp 1000.0 --nm test > resume.log 2>&1"
expect "started at run 3 rather than run 1" grep -qi 'indices from 3' resume.log
expect "three trajectories now present" test "$(count_trj testk0.5_a0.6run)" = "3"

say "5. a full ensemble is left alone"
bash "$here/setup_reactor.sh" liec.xyz -1 0.5 0.6 5 --time 1 --nrun 3 --nm test \
    > full.log 2>&1
expect "declined to re-run a complete ensemble" grep -qi 'nothing to do' full.log
expect "still three trajectories, nothing overwritten" \
    test "$(count_trj testk0.5_a0.6run)" = "3"

say "6. bad input is refused"
if bash "$here/setup_reactor.sh" missing.xyz -1 0.5 0.6 5 >/dev/null 2>&1; then
    fail "missing geometry should exit non-zero"
else
    pass "missing geometry exits non-zero"
fi

say "result"
if [ "$failed" -eq 0 ]; then
    echo "  all checks passed; working files are in $work"
else
    echo "  something failed; see above and the files in $work"
fi
exit "$failed"
