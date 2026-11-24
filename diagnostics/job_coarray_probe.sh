#!/bin/bash
#ASC --vanilla
#SBATCH --job-name=coarray_probe
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2           # 9 coarray images
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --qos=dev_zen4_0768
#SBATCH --output=coarray_probe%j.out
#SBATCH --mail-type=BEGIN,END,FAIL    # first have to state the type of event to occur
#SBATCH --mail-user=max.hirtenlehner@student.uibk.ac.at
# Load compiler and MPI environment

set -euxo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:?Submit this job from the repository directory}"
BUILD_DIR="$REPO_DIR/build/diagnostics"
PROBE="$BUILD_DIR/coarray_probe"
CONFIG="$BUILD_DIR/coarray_probe.caf"

cd "$REPO_DIR"
module load intel/2024a

echo "=== Loaded modules ==="
module list 2>&1
echo "=== MPI launcher version ==="
mpiexec -version
launch_vars="FOR_COARRAY_CONFIG_FILE LD_LIBRARY_PATH UCX_TLS FI_PROVIDER FI_PROVIDER_PATH I_MPI_ROOT I_MPI_FABRICS I_MPI_OFI_PROVIDER I_MPI_HYDRA_BOOTSTRAP I_MPI_HYDRA_BOOTSTRAP_EXEC I_MPI_PMI_LIBRARY I_MPI_PMI MPIR_CVAR_CH4_OFI_ENABLE_RMA LD_PRELOAD LD_AUDIT PMI_RANK PMI_SIZE PMI_FD PMI_PORT PMIX_RANK PMIX_NAMESPACE OMP_NUM_THREADS OMP_PROC_BIND OMP_PLACES OMP_DYNAMIC KMP_STACKSIZE SLURM_JOB_ID SLURM_NTASKS SLURM_CPUS_PER_TASK"
echo "=== Relevant variables before cleanup ==="
for name in $launch_vars; do
    if [[ -v $name ]]; then
        printf '%s=%s\n' "$name" "${!name}"
    else
        printf '%s=<unset>\n' "$name"
    fi
done
env | LC_ALL=C sort | grep -E '^(I_MPI_HYDRA_|HYDRA_)' || true
unset UCX_TLS FI_PROVIDER MPIR_CVAR_CH4_OFI_ENABLE_RMA LD_PRELOAD LD_AUDIT
export LD_LIBRARY_PATH="/cvmfs/software.asc.ac.at/versions/2025.06/software/linux/x86_64/amd/zen4/software/intel-compilers/2024.2.0/compiler/2024.2/lib:\
/cvmfs/software.asc.ac.at/versions/2025.06/software/linux/x86_64/amd/zen4/software/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib:\
/cvmfs/software.asc.ac.at/versions/2025.06/software/linux/x86_64/amd/zen4/software/imkl/2024.2.0/mkl/2024.2/lib"

if [[ ! -x "$PROBE" ]]; then
    echo "Missing probe executable: $PROBE" >&2
    echo "Build with: mpiifx -O0 -g -traceback -warn all -coarray=distributed diagnostics/coarray_probe.f90 -o build/diagnostics/coarray_probe" >&2
    exit 2
fi

printf '%s\n' "-n 2 $PROBE" > "$CONFIG"
export FOR_COARRAY_CONFIG_FILE="$CONFIG"
export I_MPI_HYDRA_BOOTSTRAP="${I_MPI_HYDRA_BOOTSTRAP:-slurm}"
export I_MPI_DEBUG=5
export OMP_NUM_THREADS=1

echo "=== Coarray probe starting ==="
echo "=== Effective launch variables ==="
for name in FOR_COARRAY_CONFIG_FILE LD_LIBRARY_PATH UCX_TLS FI_PROVIDER FI_PROVIDER_PATH I_MPI_ROOT I_MPI_FABRICS I_MPI_OFI_PROVIDER I_MPI_HYDRA_BOOTSTRAP I_MPI_HYDRA_BOOTSTRAP_EXEC I_MPI_PMI_LIBRARY I_MPI_PMI I_MPI_DEBUG OMP_NUM_THREADS LD_PRELOAD LD_AUDIT MPIR_CVAR_CH4_OFI_ENABLE_RMA SLURM_JOB_ID SLURM_NTASKS SLURM_CPUS_PER_TASK; do
    if [[ -v $name ]]; then
        printf '%s=%s\n' "$name" "${!name}"
    else
        printf '%s=<unset>\n' "$name"
    fi
done
env | LC_ALL=C sort | grep -E '^(FOR_COARRAY_|I_MPI_|MPI_|OMPI_|MPICH_|MPIR_|FI_|UCX_|PMI_|PMIX_|OMP_|KMP_|LD_|SLURM_(JOB_ID|NTASKS|CPUS_PER_TASK)=)' || true
if "$PROBE"; then
    status=0
else
    status=$?
fi
echo "COARRAY_PROBE_EXIT_STATUS=$status"
exit "$status"
