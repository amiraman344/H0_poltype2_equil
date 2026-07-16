#!/bin/bash
#SBATCH --job-name=poltype2_popc
#SBATCH --ntasks=1
#SBATCH --time=168:00:00
#SBATCH --output=output%j.out
#SBATCH --error=error%j.out
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=28GB
#SBATCH --partition=hipri
#SBATCH --qos=hipri

cd "$SLURM_SUBMIT_DIR"

source /data/home/aaamir2/miniconda3/etc/profile.d/conda.sh
conda activate poltype2-env-py310

export TINKERDIR=/data/home/aaamir2/tinker/bin/
export PSI_SCRATCH=/tab3/aman/scratch
export GAUSS_SCRDIR=/tab3/aman/scratch

mkdir -p /tab3/aman/scratch

echo "Job started:  $(date)"
echo "Host:         $(hostname)"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Molecule:     POPC (134 atoms, 44 rotatable bonds)"
echo "NOTE: This job may take 24-72 hours due to torsion scanning"

python /data/home/aaamir2/poltype2/master/PoltypeModules/poltype.py > poltype.log 2>&1

echo "Job finished: $(date)"
