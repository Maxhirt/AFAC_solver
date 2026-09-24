#!/bin/bash
#SBATCH --job-name=fac_gravity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3            # 9 coarray images
#SBATCH --cpus-per-task=40
#SBATCH --time=00:10:00
#SBATCH --output=job_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL    # first have to state the type of event to occur
# Load compiler and MPI environment

echo "job started"
echo "running on host"
hostname
ulimit -l unlimited
ulimit -s unlimited

# Optional: pin OpenMP threads if hybrid parallelism is used
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export MPIR_CVAR_CH4_OFI_ENABLE_RMA=0

# Run the distributed-memory coarray executable
export FOR_COARRAY_CONFIG_FILE=config.caf
./build/main/AFAC_solver

echo "finished"
