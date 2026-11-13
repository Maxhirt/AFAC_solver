#!/bin/bash
#ASC --vanilla
#SBATCH --job-name=fac_gravity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2           # 9 coarray images
#SBATCH --cpus-per-task=42
#SBATCH --time=00:10:00
#SBATCH --qos=dev_zen4_0768
#SBATCH --output=job_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL    # first have to state the type of event to occur
#SBATCH --mail-user=max.hirtenlehner@student.uibk.ac.at
# Load compiler and MPI environment

echo "job started"
echo "running on host"
hostname
set -euo pipefail

export OMP_NUM_THREADS=42
export FOR_COARRAY_CONFIG_FILE=/home/mh44676/AFAC_solver/config.caf
export I_MPI_HYDRA_BOOTSTRAP=slurm
export I_MPI_DEBUG=5

/home/mh44676/AFAC_solver/build/main/AFAC_solver
echo "finished"
