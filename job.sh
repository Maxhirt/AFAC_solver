#!/bin/bash
#ASC --vanilla
#SBSTCH -A p201447
#SBATCH --job-name=fac_gravity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9           # 9 coarray images
#SBATCH --cpus-per-task=16
#SBATCH --time=00:10:00
#SBATCH --qos=dev_zen4_0768
#SBATCH --output=job_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL    # first have to state the type of event t>
#SBATCH --mail-user=max.hirtenlehner@student.uibk.ac.at
# Load compiler and MPI environment

echo "job started"
echo "running on host"
hostname
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
ulimit -s unlimited
export KMP_STACKSIZE=256m
export MPIR_CVAR_CH4_OFI_ENABLE_RMA=0
/data/fs201272/mh44676/AFAC_solver/build/main/AFAC_solver
echo "finished"