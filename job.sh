#!/bin/bash
#SBATCH --job-name=fac_gravity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3            # 9 coarray images
#SBATCH --cpus-per-task=40
#SBATCH --time=00:10:00
#SBATCH -A p70652
#SBATCH --output=job_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL    # first have to state the type of event to occur
#SBATCH --mail-user=max.hirtenlehner@student.uibk.ac.at   # and then your email address
# Load compiler and MPI environment
module purge
module load cmake/3.22.0-intel-2021.5.0-25fymvk
module load intel-oneapi-compilers/2022.0.2-gcc-11.2.0-yzi4tsu
module load intel-oneapi-mpi/2021.4.0-intel-2021.5.0-jjcwtuf

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
