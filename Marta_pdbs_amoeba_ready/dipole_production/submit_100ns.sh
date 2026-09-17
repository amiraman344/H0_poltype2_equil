#!/bin/bash
# submit_100ns.sh — chain counters 1..10 to reach 100 ns total
#
# Unlike the H0_poltype2_equil reference this is based on, no counter has
# run yet here, so by default this chains ALL of 1-10 (each depending on
# the previous one via SLURM afterok), starting from traj/equil.rst.
#
# Usage:
#   bash submit_100ns.sh              # chain all missing runs (1-10)
#   bash submit_100ns.sh <start> <end>  # e.g. bash submit_100ns.sh 5 10 (resume)

START=${1:-1}
END=${2:-10}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SUBMIT="$SCRIPT_DIR/submit_dipole.sh"

echo "============================================================"
echo "  Marta H0 dipole production: chaining counters $START – $END"
echo "  Each run: 10 ns  |  Target: 100 ns total"
echo "  Submit script: $SUBMIT"
echo "============================================================"

if [ ! -f "$SCRIPT_DIR/traj/equil.rst" ]; then
    echo "  ERROR: $SCRIPT_DIR/traj/equil.rst not found."
    echo "  Run/submit equil_npt.py (submit_equil.sh) first."
    exit 1
fi

PREV_JOB=""
for CNT in $(seq $START $END); do

    OUT_DCD="$SCRIPT_DIR/traj/npt_${CNT}.dcd"
    if [ -f "$OUT_DCD" ]; then
        echo "  [skip] Counter $CNT: npt_${CNT}.dcd already exists."
        # Treat the last finished counter as the dependency anchor
        PREV_JOB=""
        continue
    fi

    if [ -z "$PREV_JOB" ]; then
        # First job in the chain (no dependency)
        JOB_ID=$(sbatch --parsable "$SUBMIT" $CNT)
    else
        # Subsequent job depends on previous one succeeding
        JOB_ID=$(sbatch --parsable --dependency=afterok:$PREV_JOB "$SUBMIT" $CNT)
    fi

    if [ $? -ne 0 ]; then
        echo "  ERROR: sbatch failed for counter $CNT. Stopping."
        exit 1
    fi

    echo "  Submitted counter $CNT  →  SLURM job $JOB_ID" \
         $([ -n "$PREV_JOB" ] && echo "(depends on $PREV_JOB)" || echo "(no dependency)")
    PREV_JOB=$JOB_ID
done

echo "============================================================"
echo "  All jobs submitted. Monitor with:"
echo "    squeue -u $USER"
echo ""
echo "  Estimated output (when complete):"
for CNT in $(seq 1 10); do
    NS=$((CNT * 10))
    echo "    After counter $CNT: $NS ns total"
done
echo "============================================================"
