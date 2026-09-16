#!/bin/bash
# Run all analysis scripts in order, each logging to its own results subfolder.
# Usage:
#   bash run_all.sh          — run all scripts (skips 03 if hole2 missing)
#   bash run_all.sh 2        — run only script 02
#
# Results layout:
#   results/01_rmsd_rmsf/
#   results/02_channel_geometry/
#   results/03_pore_radius/      (requires hole2)
#   results/04_rdf_water_his/
#   results/05_water_ocf/
#   results/06_ligand_hbonds/
#   results/07_hydrophobic_contacts/
#   results/08_ligand_diffusion/
#   results/09_induced_dipoles/

source /data/home/aaamir2/miniconda3/etc/profile.d/conda.sh
conda activate openmm_env

cd "$(dirname "$0")"

SCRIPTS=(
    01_rmsd_rmsf.py
    02_channel_geometry.py
    03_pore_radius.py
    04_rdf_water_his.py
    05_water_ocf.py
    06_ligand_hbonds.py
    07_hydrophobic_contacts.py
    08_ligand_diffusion.py
    09_induced_dipoles.py
)

run_script() {
    local SCRIPT="$1"
    local NUM="${SCRIPT:0:2}"
    local BASE="${SCRIPT%.py}"
    local LOGDIR="results/${BASE}"
    mkdir -p "$LOGDIR"
    local LOG="$LOGDIR/run.log"

    echo ""
    echo "======================================"
    echo " Running: $SCRIPT"
    echo " Log:     $LOG"
    echo " Started: $(date)"
    echo "======================================"

    python "$SCRIPT" 2>&1 | tee "$LOG"
    local STATUS=${PIPESTATUS[0]}

    if [ $STATUS -eq 0 ]; then
        echo " DONE: $SCRIPT  ($(date))"
    else
        echo " FAILED: $SCRIPT  exit=$STATUS  — check $LOG"
    fi
}

if [ -n "$1" ]; then
    IDX=$(( $1 - 1 ))
    SCRIPT="${SCRIPTS[$IDX]}"
    [ -z "$SCRIPT" ] && { echo "ERROR: script $1 not found (1-9)"; exit 1; }
    run_script "$SCRIPT"
else
    for SCRIPT in "${SCRIPTS[@]}"; do
        # Skip pore radius if hole2 not available
        if [ "$SCRIPT" = "03_pore_radius.py" ]; then
            if ! python -c "from MDAnalysis.analysis import hole2; hole2.HoleAnalysis" 2>/dev/null; then
                echo "Skipping 03_pore_radius.py (hole2 not available)"
                continue
            fi
        fi
        run_script "$SCRIPT"
    done

    echo ""
    echo "======================================"
    echo "All scripts finished."
    echo "Results:"
    find results/ -name "*.png" -o -name "*.csv" | sort
    echo "======================================"
fi
