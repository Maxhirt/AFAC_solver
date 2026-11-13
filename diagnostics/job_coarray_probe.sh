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
module load intel/2023b

echo "=== Loaded modules ==="
module list 2>&1
echo "=== MPI launcher version ==="
mpiexec -version

if [[ ! -x "$PROBE" ]]; then
    echo "Missing probe executable: $PROBE" >&2
    echo "Build with: mpiifx -O0 -g -traceback -warn all -coarray=distributed diagnostics/coarray_probe.f90 -o build/diagnostics/coarray_probe" >&2
    exit 2
fi

printf '%s\n' "-n 2 $PROBE" > "$CONFIG"
export FOR_COARRAY_CONFIG_FILE="$CONFIG"
export I_MPI_HYDRA_BOOTSTRAP=slurm
export I_MPI_DEBUG=5
export OMP_NUM_THREADS=1
unset I_MPI_FABRICS FI_PROVIDER MPIR_CVAR_CH4_OFI_ENABLE_RMA

echo "=== Coarray probe starting ==="
if "$PROBE"; then
    status=0
else
    status=$?
fi
echo "COARRAY_PROBE_EXIT_STATUS=$status"
exit "$status"
