#!/bin/bash
# Run ecrad_ifs_blocked_* on N GPUs of one node: one process per GPU, each
# computing the blocks of its part (ECRAD_NPARTS/ECRAD_PART), then merge the
# per-part outputs into OUTPUT with merge_parts.py.
#
#   usage: ecrad_multigpu.sh [-n N] [--chunk C] [--cores-per-gpu M] EXE NAMELIST INPUT OUTPUT
#
# Each process is bound to M cores and the memory of the NUMA node attached to
# its GPU (from rocm-smi --showtoponuma, else GPUs spread evenly over nodes).
# ECRAD_FULL_BATCH_SIZE and other environment settings pass through unchanged.
# PYTHON selects an interpreter with netCDF4 for the merge (default python3).

set -u

nparts=8
chunk=8
cores_per_gpu=16

usage() {
    sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n)                nparts="$2"; shift 2 ;;
        -n=*|--nparts=*)   nparts="${1#*=}"; shift ;;
        --chunk=*)         chunk="${1#*=}"; shift ;;
        --cores-per-gpu=*) cores_per_gpu="${1#*=}"; shift ;;
        -h|--help)         usage; exit 0 ;;
        --)                shift; break ;;
        -*)                echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
        *)                 break ;;
    esac
done

if [[ $# -ne 4 ]]; then
    usage >&2
    exit 1
fi
exe=$1 namelist=$2 input=$3 output=$4
python=${PYTHON:-python3}
merge=$(dirname "$(readlink -f "$0")")/merge_parts.py

nproma=$(sed -nE 's/^[[:space:]]*nblocksize[[:space:]]*=[[:space:]]*([0-9]+).*/\1/p' "$namelist" | head -1)
if [[ -z "$nproma" ]]; then
    echo "ERROR: no nblocksize in $namelist" >&2
    exit 1
fi

ngpus=$(rocm-smi --showtoponuma 2>/dev/null | grep -c 'Numa Node:')
(( ngpus > 0 )) || ngpus=$nparts
if (( nparts < 1 || nparts > ngpus )); then
    echo "ERROR: -n $nparts but this node has $ngpus GPUs" >&2
    exit 1
fi
nnodes=$(numactl -H | sed -nE 's/^available: ([0-9]+) nodes.*/\1/p')

gpu_node() {
    local node
    node=$(rocm-smi --showtoponuma 2>/dev/null | sed -nE "s/^GPU\[$1\].*Numa Node: ([0-9]+).*/\1/p")
    echo "${node:-$(( $1 * nnodes / ngpus ))}"
}

declare -A used_on_node
pids=()
for ((r = 0; r < nparts; r++)); do
    node=$(gpu_node $r)
    k=${used_on_node[$node]:-0}
    used_on_node[$node]=$((k + 1))
    cpus=$(numactl -H | sed -nE "s/^node $node cpus: (.*)/\1/p" | tr ' ' '\n' \
           | sed -n "$((k * cores_per_gpu + 1)),$(((k + 1) * cores_per_gpu))p" | paste -sd,)
    if [[ -z "$cpus" ]]; then
        echo "ERROR: not enough cores on NUMA node $node for GPU $r" >&2
        exit 1
    fi
    echo "part $r: GPU $r, NUMA node $node, cores $cpus"
    ROCR_VISIBLE_DEVICES=$r ECRAD_NPARTS=$nparts ECRAD_PART=$r ECRAD_PART_CHUNK=$chunk \
        numactl -C "$cpus" -m "$node" "$exe" "$namelist" "$input" "$output.part$r.nc" \
        > "$output.part$r.log" 2>&1 &
    pids+=($!)
done

rc=0
for ((r = 0; r < nparts; r++)); do
    if ! wait "${pids[$r]}"; then
        echo "ERROR: part $r failed, see $output.part$r.log" >&2
        rc=1
    fi
done
(( rc == 0 )) || exit $rc

# Measured iterations only (after warm-up); the slowest part sets the pace.
nwarmup=$(sed -nE 's/^[[:space:]]*nwarmup[[:space:]]*=[[:space:]]*([0-9]+).*/\1/p' "$namelist" | head -1)
for ((r = 0; r < nparts; r++)); do
    med=$(grep 'time elapsed in radiative transfer:' "$output.part$r.log" | tail -n +$(( ${nwarmup:-1} + 1 )) \
          | awk '{print $(NF-1)}' | sort -g | awk '{a[NR]=$1} END {if (NR) print a[int((NR+1)/2)]}')
    echo "part $r: median ${med:-n/a} s per iteration"
done | tee "$output.timing" | sort -k4 -g | tail -1 | awk '{print "slowest part: " $4 " s per iteration"}'

parts=()
for ((r = 0; r < nparts; r++)); do parts+=("$output.part$r.nc"); done
"$python" "$merge" --nparts "$nparts" --nproma "$nproma" --chunk "$chunk" "$output" "${parts[@]}" || exit 1
echo "merged: $output"
