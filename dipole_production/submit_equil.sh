#!/bin/bash
#SBATCH --job-name=equil_h0_dipole
#SBATCH --partition=active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem=28G
#SBATCH --time=72:00:00
#SBATCH --nodelist=tab-gpu-10-[16-19,23]
#SBATCH --output=/tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/dipole_production/out/equil_%j.log
#SBATCH --error=/tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/dipole_production/out/equil_%j.err

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "GPU: $CUDA_VISIBLE_DEVICES"

source /data/home/aaamir2/miniconda3/etc/profile.d/conda.sh
conda activate openmm_env

cd /tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/dipole_production
python equil_npt.py

echo "Job finished: $(date)"
