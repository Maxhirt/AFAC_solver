program main
    use setup, only: max_iterations, epsilon
    use afac, only: setup_afac, error_projection, residual, l2_norm_residual, multigrid_afac
    use output, only: output_results_global
    use setup_initial_condition, only: initial_condition

    implicit none(type, external)
    integer :: num_runs
    double precision :: l2_norm

    call initial_condition()

    call setup_afac()
    sync all
    if (THIS_IMAGE() == 1) then
        print *, "Beginning first residual calculation "
    end if
    call residual()
    call l2_norm_residual(l2_norm)
    if (THIS_IMAGE() == 1) then
        print *, "Starting calculation with l2_error:  ", l2_norm
    end if
    num_runs = 0
    do while ((l2_norm > epsilon) .and. (num_runs < max_iterations))
        call multigrid_afac
        sync all
        call error_projection

        call residual()

        call l2_norm_residual(l2_norm)

        if (THIS_IMAGE() == 1) then
            print *, "l2_error:  ", l2_norm
        end if
        sync all
        num_runs = num_runs + 1

    end do
    sync all
    call output_results_global(num_runs)
end program main
