#!/bin/bash
#SBATCH --job-name=npt_h0p2
#SBATCH --partition=active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem=28G
#SBATCH --time=72:00:00
#SBATCH --nodelist=tab-gpu-10-[16-19,23]
#SBATCH --output=/tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/out/npt_run_%j.log
#SBATCH --error=/tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/out/npt_run_%j.err

# Usage: sbatch submit_npt.sh <counter>
# counter=1 : reads traj/equil.rst → 10 ns production
# counter=N : reads traj/npt_{N-1}.rst → 10 ns production
CNT=${1:-1}

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "Run counter: $CNT"

source /data/home/aaamir2/miniconda3/etc/profile.d/conda.sh
conda activate openmm_env

cd /tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil
python amoeba_npt.py $CNT

echo "Job finished: $(date)"
