#!/bin/bash
#SBATCH --job-name=equil_h0_dipole_prod
#SBATCH --partition=hipri
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=28G
#SBATCH --time=72:00:00
#SBATCH --qos=hipri
#SBATCH --output=/tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/dipole_production/out/slurm_%j.log

# Usage: sbatch submit_dipole.sh <counter>
#   counter=1  : start from traj/equil.rst (end of restrained equilibration)
#   counter=N  : restart from this folder's traj/npt_{N-1}.rst
CNT=${1:-1}

echo "========================================"
echo "Job     : $SLURM_JOB_ID"
echo "Node    : $(hostname)"
echo "GPU     : $CUDA_VISIBLE_DEVICES"
echo "Counter : $CNT"
echo "Started : $(date)"
echo "========================================"

source /data/home/aaamir2/miniconda3/etc/profile.d/conda.sh
conda activate openmm_env

cd /tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/dipole_production

python amoeba_npt_dipole.py $CNT

echo "========================================"
echo "Finished: $(date)"
echo "========================================"
