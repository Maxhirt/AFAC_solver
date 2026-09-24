#!/bin/bash

# Exit immediately if no argument is provided
if [ -z "$1" ]; then
  echo "Error: No option specified."
  echo "Usage: ./results.sh <option>"
  exit 1
fi

OPTION="$1"

# Create required directories
mkdir -p "build_${OPTION}"
mkdir -p "${OPTION}"

# Create and populate the SLURM job script
cat << EOF > "job_${OPTION}.sh"
#!/bin/bash
#SBATCH --job-name=${OPTION}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=40
#SBATCH --time=00:10:00
#SBATCH -A p70652
#SBATCH --output=job_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=max.hirtenlehner@student.uibk.ac.at
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
export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export MPIR_CVAR_CH4_OFI_ENABLE_RMA=0

# Run the distributed-memory coarray executable
export FOR_COARRAY_CONFIG_FILE=config.caf
./build_${OPTION}/main/AFAC_solver
echo "finished"
EOF

# Make the created job script executable
chmod +x "job_${OPTION}.sh"

# Overwrite config.caf with the updated path
echo "-n 2 ./${OPTION}/main/AFAC_solver" > config.caf

# Enter build folder, run CMake & Make, then return to parent directory
cd "build_${OPTION}" || exit 1
FC=ifort cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
cd ..

# Submit the generated SLURM job script
sbatch "job_${OPTION}.sh"

echo "Setup complete, built in build_${OPTION}, and submitted job_${OPTION}.sh."
