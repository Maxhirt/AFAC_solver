program main
    use setup, only: max_iterations, epsilon, global_domain_length, domain_length, N, hloc
    use afac, only: setup_afac, error_projection, residual, l2_norm_residual, multigrid_afac
    use output, only: output_results_global
    use setup_initial_condition, only: initial_condition
    use analytical_solution, only: analytical_solution_calculator
    implicit none(type, external)
    integer :: num_runs
    double precision :: l2_norm, old_l2_norm
    double precision :: h

    h = global_domain_length/N
    hloc = h/(2.d0**(THIS_IMAGE() - 1))
    domain_length = global_domain_length/(2.d0**(THIS_IMAGE() - 1))


    call initial_condition()

    call setup_afac()
    sync all
    call residual()
    call l2_norm_residual(l2_norm)
    num_runs = 0
    old_l2_norm = 20*l2_norm
    do while ((abs(l2_norm-old_l2_norm) > epsilon) .and. (num_runs < max_iterations))
        call multigrid_afac
        sync all
        call error_projection

        call residual()
        old_l2_norm =l2_norm

        call l2_norm_residual(l2_norm)
        sync all
        num_runs = num_runs + 1
        print *, "num_runs: ", num_runs

    end do
    sync all


end program main
