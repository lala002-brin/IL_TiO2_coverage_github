#!/bin/bash
#SBATCH --job-name=dtpa_tio2_relax
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --output=output/relax_slurm-%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# module purge
# module load quantum-espresso
mkdir -p tmp_relax output
echo "Host: $(hostname)"; echo "Start: $(date)"; echo "pw.x: $(which pw.x)"
srun pw.x -in output/02_DTPA_TiO2_relax.in > output/02_DTPA_TiO2_relax.out
echo "End: $(date)"
