#!/bin/bash
# Run Poltype2 AMOEBA parameterization for ligand
# Execute: bash run_poltype.sh

cd "$(dirname "$0")"

source /data/home/aaamir2/miniconda3/etc/profile.d/conda.sh
conda activate poltype2-env-py310

export TINKERDIR=/data/home/aaamir2/tinker/bin/
export PSI_SCRATCH=/tab1/aaamir2/scratch
export GAUSS_SCRDIR=/tab1/aaamir2/scratch

python /data/home/aaamir2/poltype2/master/PoltypeModules/poltype.py > poltype.log 2>&1 &
echo "Poltype2 started with PID $!"
echo "Monitor progress: tail -f poltype.log"
