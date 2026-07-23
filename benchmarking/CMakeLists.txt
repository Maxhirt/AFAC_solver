set(RESOLUTIONS 16 32 64 128 256)

# Declare custom target to run/manage benchmarks
add_custom_target(run_benchmarks)

# Optional: Master submission script to submit all slurm jobs with one command
set(SUBMIT_ALL_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/submit_all_jobs.sh")
file(WRITE ${SUBMIT_ALL_SCRIPT} "#!/bin/bash\n# Master Slurm Submission Script\n\n")

foreach(RES IN LISTS RESOLUTIONS)
    set(TARGET_NAME benchmark_${RES})

    # Define executable using absolute path to top-level src/ directory
    add_executable(${TARGET_NAME}
        ${CMAKE_SOURCE_DIR}/src/setup.f90
        ${CMAKE_SOURCE_DIR}/src/communication.f90
        ${CMAKE_SOURCE_DIR}/src/afac.f90
        ${CMAKE_SOURCE_DIR}/src/boundary.f90
        ${CMAKE_SOURCE_DIR}/src/setup_initial_condition.f90
        ${CMAKE_SOURCE_DIR}/src/multigrid_solver.f90
        benchmark.f90
    )

    # Link forbenchmark library if used inside benchmark.f90
    if(TARGET forbenchmark)
        target_link_libraries(${TARGET_NAME} PRIVATE forbenchmark)
    endif()

    set_target_properties(${TARGET_NAME} PROPERTIES
        Fortran_PREPROCESS ON
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mods_${RES}"
    )

    # Pass N_GRID macro definition to compiler
    target_compile_definitions(${TARGET_NAME} PRIVATE GRID_N=${RES})

    add_dependencies(run_benchmarks ${TARGET_NAME})

    # 1. Generate config_<RES>.caf file for Intel Coarrays
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/config_${RES}.caf"
"-n 9 ./${TARGET_NAME}
"
    )

    # 2. Generate Slurm batch script (Unindented to prevent bash syntax errors)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/slurm_${RES}.sh"
"#!/bin/bash
#SBATCH --job-name=fac_mg_${RES}
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=40
#SBATCH --time=01:30:00
#SBATCH -A p70652
#SBATCH --output=job_${RES}_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=max.hirtenlehner@student.uibk.ac.at

# Load environment
module purge
module load cmake/3.22.0-intel-2021.5.0-25fymvk
module load intel-oneapi-compilers/2022.0.2-gcc-11.2.0-yzi4tsu
module load intel-oneapi-mpi/2021.4.0-intel-2021.5.0-jjcwtuf

ulimit -l unlimited
ulimit -s unlimited

# OpenMP / Parallel settings
export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export MPIR_CVAR_CH4_OFI_ENABLE_RMA=0

# Intel Coarray Runtime Config File Selection
export FOR_COARRAY_CONFIG_FILE=config_${RES}.caf

# Execute benchmark binary
./${TARGET_NAME}
"
    )

    # Append to master submission script
    file(APPEND ${SUBMIT_ALL_SCRIPT} "sbatch slurm_${RES}.sh\n")

endforeach()

# Make master submission script executable
execute_process(COMMAND chmod +x ${SUBMIT_ALL_SCRIPT})
