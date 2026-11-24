#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# Activate your DFT environment and load your licensed VASP/MPI installation.
# Prepare the work directory and PAW inputs first; choose a single phase here.
python -m dana_kappa.vasp --phase fc2 --workdir runs/example_dft \
  --vasp-command "srun vasp_std"
