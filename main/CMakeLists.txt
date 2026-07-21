add_executable(AFAC_solver
    main.f90
)
target_link_libraries(AFAC_solver PRIVATE AFAC_solver_lib)
