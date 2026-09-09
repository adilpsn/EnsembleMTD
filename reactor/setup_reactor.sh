#!/bin/bash
# Run an ensemble of RMSD-driven metadynamics trajectories on one starting
# geometry with GFN2-xTB.
#
# Each ensemble member is an independent MD relaxation followed by a biased MTD
# run, so the members differ only through the thermostat's random velocities and
# the relaxed geometry that comes out of them.  That is the whole source of the
# spread the ensemble analysis then measures.
#
# Per member:
#   1. crest --reactor --genpot rho --genmtd 0.5  writes an rcontrol template
#   2. xtb --md                                   0.5 ps unbiased relaxation
#   3. crest --reactor --genpot rho --genmtd t     rcontrol for the biased run
#   4. xtb --md with $metadyn                      the production MTD trajectory
#
# crest picks kpush from the system size; step 3 then rescales it by the
# requested factor (and by the fraction of biased atoms under --nohbias) before
# xtb sees it.  The scaling is kpush * factor * bias_ratio / 4 -- the /4 is what
# maps crest's conformer-search default onto the reactive window used here.
#
# Existing k<kpush>_a<alp>run*.trj files in the working directory are counted
# first, so re-running with a larger --nrun tops the ensemble up instead of
# starting over.  The MTD parameters are written into line 2 of every .trj:
# that comment, not the rcontrol left behind, is the authoritative record of
# how a trajectory was produced (crest rewrites rcontrol on every member).
#
# This is a plain program with no scheduler directives in it, so it works the
# same whether you run it directly or submit it.  A production ensemble is long
# (20 trajectories x 30 ps is hours to days) and xtb is threaded, so it wants a
# whole node and a generous wall clock.  Copy submit.example.sbatch and put your
# site's partition, account and module loads there.
#
# xtb picks up its thread count from OMP_NUM_THREADS.

usage() {
    echo "Usage: $0 <xyz_file> <charge> <kpush_factor> <alp_value> <density> [--time <ps>] [--mddump <val>] [--nrun <N>] [--atoms <range>] [--walltemp <K>] [--alpb <solvent>] [--timestep <fs>] [--nohbias] [--nm <name>] [--noclutter] [--dryrun]"
    echo "Example: $0 liec2.xyz 0 0.5 0.6 5 --time 30 --mddump 100 --nrun 20 --walltemp 1000.0 --nm liec2"
    echo "         (k/N=0.5, alpha=0.6, density=5 are the production settings)"
    exit 1
}

# ---------- Defaults ----------
simulation_time=30
mddump=1000
num_simulations=1
save=20
step_value=0.5
atoms_list=""
walltemp=1000.0
solvent=""
nohbias=false
dryrun=false
noclutter=false
run_name=""

# ---------- Positional ----------
if [ "$#" -lt 5 ]; then usage; fi
input_file="$1"; shift
charge="$1"; shift
kpush_factor="$1"; shift
alp_value="$1"; shift
genpot_value="$1"; shift

# ---------- Optional flags ----------
while [[ $# -gt 0 ]]; do
    case $1 in
        --time) simulation_time="$2"; shift 2 ;;
        --mddump) mddump="$2"; shift 2 ;;
        --nrun) num_simulations="$2"; shift 2 ;;
        --atoms) atoms_list="$2"; shift 2 ;;
        --walltemp) walltemp="$2"; shift 2 ;;
        --alpb) solvent="$2"; shift 2 ;;
        --timestep) step_value="$2"; shift 2 ;;
        --nohbias) nohbias=true; shift ;;
        --dryrun) dryrun=true; shift ;;
        --noclutter) noclutter=true; shift ;;
        --nm) run_name="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# ---------- Validation ----------
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

echo "----------------------------------------------"
echo "  Input file:        $input_file"
echo "  Charge:            $charge"
echo "  kpush factor:      $kpush_factor"
echo "  alpha (alp):       $alp_value"
echo "  Density:           $genpot_value g/cm³"
echo "  Simulation time:   ${simulation_time} ps"
echo "  mddump:            ${mddump}"
echo "  Num. simulations:  ${num_simulations}"
echo "  Wall temperature:  ${walltemp} K"
[ -n "$solvent" ] && echo "  Solvent:           $solvent"
[ "$nohbias" = true ] && echo "  Excluding H from bias"
[ -n "$atoms_list" ] && echo "  Metadyn atoms:     $atoms_list"
[ -n "$run_name" ] && echo "  Run name prefix:   $run_name"
[ "$noclutter" = true ] && echo "  Noclutter mode:    ON"
echo "----------------------------------------------"

# ---------- Compute atoms if --nohbias ----------
if [ "$nohbias" = true ]; then
    total_atoms=$(awk 'NR==1 {print $1}' "$input_file")
    nonh_indices=$(awk 'NR>2 && $1!="H" {printf NR-2","}' "$input_file" | sed 's/,$//')
    atoms_list="$nonh_indices"
    if [ -n "$atoms_list" ]; then
        num_bias_atoms=$(echo "$atoms_list" | awk -F',' '{print NF}')
        bias_ratio=$(echo "scale=6; $num_bias_atoms / $total_atoms" | bc)
        echo " Biasing $num_bias_atoms / $total_atoms atoms (ratio=$bias_ratio)"
    else
        bias_ratio=1.0
    fi
else
    bias_ratio=1.0
fi

# ---------- Helpers ----------
# delete a $block safely (until next $-header, not including it)
delete_block () {
    local block="$1"
    local start end
    start=$(grep -n "^\$$block" rcontrol | cut -d: -f1)
    if [ -n "$start" ]; then
        end=$(awk -v s="$start" 'NR>s && /^\$[a-zA-Z]/ {print NR; exit}' rcontrol)
        if [ -z "$end" ]; then end=$(wc -l < rcontrol); fi
        sed -i "${start},$((end-1))d" rcontrol
    fi
}

# ensure $wall exists and has temp=walltemp; if missing, insert before $end
ensure_wall () {
    if grep -q "^\$wall" rcontrol; then
        sed -i -E "s/(^|[[:space:]])temp=[0-9.]+/\1temp=${walltemp}/" rcontrol
    else
        endline=$(grep -n '^\$end' rcontrol | cut -d: -f1 | head -n1)
        wall_block="\$wall
  potential=logfermi
  sphere: auto, all
  temp=${walltemp}"
        if [ -n "$endline" ]; then
            awk -v e="$endline" -v w="$wall_block" 'NR==e{print w}1' rcontrol > rcontrol.tmp && mv rcontrol.tmp rcontrol
        else
            printf "%s\n\$end\n" "$wall_block" >> rcontrol
        fi
    fi
}

# ---------- DRY RUN: just generate rcontrol template and exit ----------
if [ "$dryrun" = true ]; then
    echo " Dry run — generating rcontrol for MTD only."
    crest "$input_file" --reactor --genpot "$genpot_value" --genmtd "$simulation_time"
    status=$?
    if [ $status -ne 0 ]; then
        echo " CREST dry-run failed (status $status)."
        exit $status
    fi

    [ ! -f rcontrol ] && { echo "Error: rcontrol not found after CREST!"; exit 1; }

    sed -i "1s;^;\$cma\n;" rcontrol

    current_kpush=$(sed -n 's/.*kpush=\([0-9.eE+-]*\).*/\1/p' rcontrol | head -n 1)
    [ -z "$current_kpush" ] && { echo "Error: could not find kpush in rcontrol"; exit 1; }
    new_kpush=$(echo "scale=6; $current_kpush * $kpush_factor * $bias_ratio / 4" | bc)
    echo "kpush: $current_kpush -> $new_kpush (factor=$kpush_factor, ratio=$bias_ratio)"

    sed -i "s/kpush=$current_kpush/kpush=$new_kpush/" rcontrol
    sed -i -E "s/step=[0-9.]+/step=$step_value/" rcontrol
    sed -i -E "s/alp=[0-9.]+/alp=$alp_value/" rcontrol
    sed -i -E "s/save=[0-9.]+/save=$save/" rcontrol

    if grep -q "^\$set" rcontrol; then
        if grep -q "mddump" rcontrol; then
            sed -i -E "s/(mddump[[:space:]]+)[0-9]+/\1${mddump}/" rcontrol
        else
            sed -i "/^\$set/a\  mddump  ${mddump}" rcontrol
        fi
    fi

    ensure_wall

    if [ -n "$atoms_list" ]; then
        delete_block "metadyn"
        insert_line=$(grep -n '^\$wall' rcontrol | cut -d: -f1 | head -n1)
        if [ -z "$insert_line" ]; then
            insert_line=$(grep -n '^\$end' rcontrol | cut -d: -f1 | head -n1)
            [ -z "$insert_line" ] && insert_line=$(($(wc -l < rcontrol) + 1))
        fi
        metadyn_block="\$metadyn
  atoms: $atoms_list
  save=$save
  kpush=$new_kpush
  alp=$alp_value
  bias_ramp_time=0.01"
        awk -v il="$insert_line" -v blk="$metadyn_block" 'NR==il{print blk}1' rcontrol > rcontrol.tmp && mv rcontrol.tmp rcontrol
    fi

    echo " Dry run complete. rcontrol ready."
    exit 0
fi

# ---------- Determine starting run index based on existing files ----------
# Pattern: <run_name>k<kpush>_a<alp>runX.trj, e.g. k0.5_a0.6run3.trj or
# lieck0.5_a0.6run3.trj.  Only the * is left unquoted, so a prefix containing
# spaces still matches; the -e test stands in for nullglob.
max_existing=0
run_prefix="${run_name}k${kpush_factor}_a${alp_value}run"
existing=()
for candidate in "$run_prefix"*.trj; do
    [ -e "$candidate" ] && existing+=("$candidate")
done

if (( ${#existing[@]} > 0 )); then
    for f in "${existing[@]}"; do
        idx=${f##*run}
        idx=${idx%.trj}
        [[ "$idx" =~ ^[0-9]+$ ]] || continue
        (( idx > max_existing )) && max_existing=$idx
    done
    if (( max_existing > 0 )); then
        echo " Found ${max_existing} existing run(s) for this k, a, and prefix '${run_name}'."
    fi
fi

start_index=$(( max_existing + 1 ))

if (( start_index > num_simulations )); then
    echo " Requested nrun=${num_simulations}, but ${max_existing} runs already exist. Nothing to do."
    exit 0
fi

echo " Will run indices from ${start_index} to ${num_simulations}."

# ---------- MAIN LOOP: MD relaxation + MTD per run ----------
for ((i=start_index; i<=num_simulations; i++)); do
    echo "=============================================="
    echo "Run $i / $num_simulations"
    echo "=============================================="

    # ----- 1) Short MD relaxation for this run -----
    echo " MD pre-relaxation for run $i ..."
    crest "$input_file" --reactor --genpot "$genpot_value" --genmtd 0.5
    status=$?
    if [ $status -ne 0 ]; then
        echo " CREST MD pre-relaxation failed or was killed (status $status). Stopping."
        exit $status
    fi

    [ ! -f rcontrol ] && { echo "Error: rcontrol not found after CREST (MD relax)!"; exit 1; }

    # ensure $cma
    sed -i "1s;^;\$cma\n;" rcontrol
    sed -i -E "s/step=[0-9.]+/step=$step_value/" rcontrol

    # MD: set shake=t
    if grep -q "^shake=" rcontrol; then
        sed -i "s/^shake=.*/shake=t/" rcontrol
    else
        sed -i "/^\$md/a shake=t" rcontrol
    fi

    # keep $wall; just ensure temp
    ensure_wall

    # for relaxation we do *not* want metadynamics
    delete_block "metadyn"

    ns_relax="${run_name}moldy_relax_${i}"
    xtb_cmd=(xtb "$input_file" --md --input rcontrol --namespace "$ns_relax" -c "$charge")
    [ -n "$solvent" ] && xtb_cmd+=(--alpb "$solvent")

    echo " MD relax cmd: ${xtb_cmd[*]}"
    "${xtb_cmd[@]}" > "md_relax_run${i}.out"
    status=$?
    if [ $status -ne 0 ]; then
        echo " xTB MD pre-relaxation failed or was killed (status $status). Stopping."
        exit $status
    fi

    trj_relax="${ns_relax}.xtb.trj"
    md_traj="${run_name}md_run${i}.xyz"

    if [ -f "$trj_relax" ]; then
        # rename MD trajectory to a nicer name and keep it
        mv "$trj_relax" "$md_traj"
        na=$(head -n 1 "$md_traj")
        lpf=$((na + 2))
        tail -n "$lpf" "$md_traj" > "relaxed_run${i}.xyz"
        start_xyz="relaxed_run${i}.xyz"
    else
        echo "Error: '$trj_relax' not found!"
        exit 1
    fi

    # clean MD control + moldy files
    rm -f rcontrol moldy*

    # ----- 2) Generate MTD control for this relaxed structure -----
    echo "Generating MTD control for run $i ..."
    crest "$start_xyz" --reactor --genpot "$genpot_value" --genmtd "$simulation_time"
    status=$?
    if [ $status -ne 0 ]; then
        echo " CREST MTD control generation failed or was killed (status $status). Stopping."
        exit $status
    fi

    [ ! -f rcontrol ] && { echo "Error: rcontrol not found after CREST (MTD)!"; exit 1; }

    # ensure $cma
    sed -i "1s;^;\$cma\n;" rcontrol

    # scale kpush
    current_kpush=$(sed -n 's/.*kpush=\([0-9.eE+-]*\).*/\1/p' rcontrol | head -n 1)
    [ -z "$current_kpush" ] && { echo "Error: could not find kpush in rcontrol"; exit 1; }
    new_kpush=$(echo "scale=6; $current_kpush * $kpush_factor * $bias_ratio / 4" | bc)
    echo "kpush: $current_kpush -> $new_kpush (factor=$kpush_factor, ratio=$bias_ratio)"

    sed -i "s/kpush=$current_kpush/kpush=$new_kpush/" rcontrol
    sed -i -E "s/step=[0-9.]+/step=$step_value/" rcontrol
    sed -i -E "s/alp=[0-9.]+/alp=$alp_value/" rcontrol
    sed -i -E "s/save=[0-9.]+/save=$save/" rcontrol

    # mddump in $set
    if grep -q "^\$set" rcontrol; then
        if grep -q "mddump" rcontrol; then
            sed -i -E "s/(mddump[[:space:]]+)[0-9]+/\1${mddump}/" rcontrol
        else
            sed -i "/^\$set/a\  mddump  ${mddump}" rcontrol
        fi
    fi

    # ensure $wall exists & has temp
    ensure_wall

    # replace $metadyn content (without adding $end!)
    if [ -n "$atoms_list" ]; then
        delete_block "metadyn"
        insert_line=$(grep -n '^\$wall' rcontrol | cut -d: -f1 | head -n1)
        if [ -z "$insert_line" ]; then
            insert_line=$(grep -n '^\$end' rcontrol | cut -d: -f1 | head -n1)
            [ -z "$insert_line" ] && insert_line=$(($(wc -l < rcontrol) + 1))
        fi

        metadyn_block="\$metadyn
  atoms: $atoms_list
  save=$save
  kpush=$new_kpush
  alp=$alp_value
  bias_ramp_time=0.01"

        awk -v il="$insert_line" -v blk="$metadyn_block" 'NR==il{print blk}1' rcontrol > rcontrol.tmp && mv rcontrol.tmp rcontrol
    fi

    # ----- 3) Run MTD -----
    echo " Running MTD simulation $i ..."
    xtb_cmd=(xtb "$start_xyz" --md --input rcontrol -c "$charge")
    [ -n "$solvent" ] && xtb_cmd+=(--alpb "$solvent")

    "${xtb_cmd[@]}" > "mtd_run${i}.out"
    status=$?
    if [ $status -ne 0 ]; then
        echo " xTB MTD run failed or was killed (status $status). Stopping."
        exit $status
    fi

    # ----- 4) Log run parameters into outputs -----
    info="\nkpush factor: ${kpush_factor}\nalp_value: ${alp_value}\nDensity: ${genpot_value}\nmddump: ${mddump}\nwalltemp: ${walltemp}\nsolvent: ${solvent}\nbias_ratio: ${bias_ratio}"
    echo -e "$info" >> "mtd_run${i}.out"
    echo -e "$info" >> "md_relax_run${i}.out"

    # Also write parameters into comment line (2nd line) of xtb.trj
    if [ -f xtb.trj ]; then
        comment="MTD params: charge=${charge} kpush_factor=${kpush_factor}, alp_value=${alp_value}, Density=${genpot_value}, mddump=${mddump}, walltemp=${walltemp}, solvent=${solvent}, bias_ratio=${bias_ratio}"
        awk -v c="$comment" 'NR==1{print; next} NR==2{print c; next} {print}' xtb.trj > xtb_withcomment.trj
        mv xtb_withcomment.trj xtb.trj
    fi

    out="${run_name}k${kpush_factor}_a${alp_value}run${i}.trj"
    cp xtb.trj "$out"
    echo " MTD simulation $i done → $out"

    # clean per-run temp files
    rm -f xtbtop* wbo coord* scoord* xtbmdok mdrestart xtbrestart

    if [ "$noclutter" = true ]; then
        # keep only MTD out + final .trj
        rm -f "md_relax_run${i}.out" "relaxed_run${i}.xyz" "$md_traj" rcontrol moldy* xtb.trj
    else
        # default: keep md trajectory + relaxed structure; remove md.out
        rm -f "md_relax_run${i}.out"
        rm -f moldy* 
        # keep: $md_traj and relaxed_run${i}.xyz, plus mtd_run${i}.out and final .trj
    fi
done

# ---------- Products from LAST run ----------
# (If you want products from each run, you can move this into the loop.)
out="${run_name}k${kpush_factor}_a${alp_value}run${num_simulations}.trj"
obabel -i xyz "$out" -o smi -O traj.smi

echo " All $num_simulations MTD simulations completed successfully."
