add_executable(AFAC_solver
    main.f90
)
target_link_libraries(AFAC_solver PRIVATE AFAC_solver_lib)
target_compile_options(AFAC_solver PRIVATE ${CMAKE_Fortran_FLAGS})
