#!/bin/bash
#ASC --vanilla
#SBATCH --job-name=fac_gravity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2           # 9 coarray images
#SBATCH --cpus-per-task=16
#SBATCH --time=00:10:00
#SBATCH --qos=dev_zen4_0768
#SBATCH --output=job_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL    # first have to state the type of event to occur
#SBATCH --mail-user=max.hirtenlehner@student.uibk.ac.at
module load intel/2024a

echo "job started"
echo "running on host"
hostname

set -euo pipefail

# Explicitly set LD_LIBRARY_PATH for runtime using the CVMFS paths from the module
export LD_LIBRARY_PATH="/cvmfs/software.asc.ac.at/versions/2025.06/software/linux/x86_64/amd/zen4/software/intel-compilers/2024.2.0/compiler/2024.2/lib:\
/cvmfs/software.asc.ac.at/versions/2025.06/software/linux/x86_64/amd/zen4/software/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib:\
/cvmfs/software.asc.ac.at/versions/2025.06/software/linux/x86_64/amd/zen4/software/imkl/2024.2.0/mkl/2024.2/lib:\
${LD_LIBRARY_PATH:-}"

echo "LD_LIBRARY_PATH is now: $LD_LIBRARY_PATH"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_NUM_THREADS=16
export FOR_COARRAY_CONFIG_FILE=/data/fs201272/mh44676/AFAC_solver/config.caf
ulimit -s unlimited
export KMP_STACKSIZE=256m
# Run your coarray solver executable
/data/fs201272/mh44676/AFAC_solver/build/main/AFAC_solver
echo "finished"