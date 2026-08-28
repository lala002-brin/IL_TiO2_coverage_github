#!/bin/bash
#SBATCH --job-name=dtpa_tio2_scf
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=12:00:00
#SBATCH --output=output/scf_slurm-%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# module purge
# module load quantum-espresso
mkdir -p tmp output
echo "Host: $(hostname)"; echo "Start: $(date)"; echo "pw.x: $(which pw.x)"
srun pw.x -in output/01_DTPA_TiO2_scf.in > output/01_DTPA_TiO2_scf.out
echo "End: $(date)"
