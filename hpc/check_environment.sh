#!/bin/bash
set -u
echo '=== Python ==='; python3 --version || true
echo '=== pip ==='; python3 -m pip --version || true
echo '=== ASE ==='; python3 -c "import ase; print(ase.__version__)" || true
echo '=== pw.x ==='; which pw.x || true; pw.x --version 2>/dev/null || true
echo '=== Scheduler ==='; which sbatch || true; which qsub || true
