add_library(AFAC_solver_lib
    afac.f90
    analytical_solution.f90
    boundary.f90
    communication.f90
    multigrid_solver.f90
    output.f90
    setup.f90
)

set_target_properties(AFAC_solver_lib PROPERTIES
    Fortran_PREPROCESS ON
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)
target_include_directories(AFAC_solver_lib PUBLIC ${CMAKE_BINARY_DIR}/modules)
