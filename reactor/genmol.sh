#!/bin/bash
# Build a starting geometry for the metadynamics reactor from a dot-separated
# SMILES string.  Each fragment is embedded with Open Babel --gen3d and the set
# is then packed by Packmol into a cube or a sphere whose size is set by the
# atom count and a packing parameter rho, so the initial fragment separation is
# comparable across systems of different size.
#
#   ./genmol.sh "[Li+].C1COC(=O)O1" liec
#   ./genmol.sh "C1COC(=O)O1.C1COC(=O)O1" 2ec --shape box --density 0.13
#
# The sphere is the shape used for the production ensembles; the cube is kept
# for the earlier screening runs.

set -u

usage() {
    cat <<'USAGE'
Usage: genmol.sh "SMILES[.SMILES...]" OUTPUT_NAME [--shape sphere|box] [--density RHO] [--tolerance D]

  SMILES        dot-separated SMILES, one fragment per molecule
  OUTPUT_NAME   basename of the packed .xyz written next to the current directory

  --shape       sphere (default) or box
  --density     packing parameter rho (default 0.6 for sphere, 0.13 for box).
                Cube edge is L = (N/rho)^(1/3), so rho is a number density in
                atoms/Ang^3 there.  The sphere uses r = (3N/(pi*rho))^(1/3),
                which is the expression the published ensembles were built with;
                note the missing factor of 4 relative to a true number density,
                so the sphere is looser than rho alone suggests.  Kept as-is for
                reproducibility -- rho is a knob, not a measured density.
  --tolerance   Packmol minimum interatomic distance in Ang
                (default 2.2 for sphere, 2.0 for box)
USAGE
    exit 1
}

[ "$#" -ge 2 ] || usage
smiles_string="$1"; shift
output_name="${1}.xyz"; shift

shape="sphere"
density=""
tolerance=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        --shape)     shape="$2"; shift 2 ;;
        --density)   density="$2"; shift 2 ;;
        --tolerance) tolerance="$2"; shift 2 ;;
        -h|--help)   usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

case "$shape" in
    sphere) : "${density:=0.6}";  : "${tolerance:=2.2}" ;;
    box)    : "${density:=0.13}"; : "${tolerance:=2.0}" ;;
    *) echo "Unknown shape '$shape' (expected sphere or box)" >&2; exit 1 ;;
esac

for tool in obabel packmol bc; do
    command -v "$tool" >/dev/null || { echo "Required tool not on PATH: $tool" >&2; exit 1; }
done

target_dir="$PWD"
workdir=$(mktemp -d "${TMPDIR:-/tmp}/genmol.XXXXXX") || exit 1
trap 'rm -rf "$workdir"' EXIT
cd "$workdir" || exit 1

# Embed each fragment separately.  Packmol needs the atom counts to size the
# container, and reading them back from the .xyz avoids re-parsing the SMILES.
IFS='.' read -ra fragments <<< "$smiles_string"
mol_files=()
total_atoms=0
for i in "${!fragments[@]}"; do
    mol_file="mol_$((i + 1)).xyz"
    echo "Embedding $mol_file from ${fragments[i]}"
    if ! obabel -:"${fragments[i]}" -oxyz -O "$mol_file" --gen3d 2>/dev/null; then
        echo "Open Babel could not embed fragment: ${fragments[i]}" >&2
        exit 1
    fi
    n=$(head -n 1 "$mol_file" | tr -d '[:space:]')
    case "$n" in
        ''|*[!0-9]*) echo "No atom count in $mol_file for ${fragments[i]}" >&2; exit 1 ;;
    esac
    total_atoms=$((total_atoms + n))
    mol_files+=("$mol_file")
done

if [ "$total_atoms" -le 0 ]; then
    echo "Total atom count is zero -- check the SMILES input" >&2
    exit 1
fi

# Container size from N and rho.  bc has no pi, so a(1) = atan(1) = pi/4 and
# 4*a(1) = pi.  See --density in the usage text for why the sphere expression
# carries pi rather than 4*pi.
if [ "$shape" = "sphere" ]; then
    radius=$(echo "scale=2; e(l((3 * $total_atoms) / (4 * a(1) * $density)) / 3)" | bc -l)
    region="inside sphere 0. 0. 0. $radius"
    echo "$total_atoms atoms, sphere radius $radius Ang (rho = $density)"
else
    length=$(echo "scale=2; e(l($total_atoms / $density) / 3)" | bc -l)
    region="inside box 0. 0. 0. $length $length $length"
    echo "$total_atoms atoms, cube edge $length Ang (rho = $density atoms/Ang^3)"
fi

{
    echo "tolerance $tolerance"
    echo "filetype xyz"
    echo "output packed_structure.xyz"
    for mol_file in "${mol_files[@]}"; do
        echo "structure $mol_file"
        echo "  number 1"
        echo "  $region"
        echo "end structure"
    done
} > packmol.inp

if ! packmol < packmol.inp > packmol.log 2>&1; then
    echo "Packmol failed; see the log below" >&2
    tail -n 20 packmol.log >&2
    exit 1
fi

[ -f packed_structure.xyz ] || { echo "Packmol wrote no packed_structure.xyz" >&2; exit 1; }

mv packed_structure.xyz "$target_dir/$output_name"
echo "Wrote $output_name"
