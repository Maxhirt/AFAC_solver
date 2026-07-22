program benchmark
    use setup, only: max_iterations, epsilon, global_domain_length, domain_length, N, hloc
    use afac, only: setup_afac, error_projection, residual, l2_norm_residual, multigrid_afac
    use setup_initial_condition, only: initial_condition
    implicit none(type, external)
    integer :: num_runs, t_start, t_end, t_rate, unit_num
    double precision :: l2_norm, elapsed_seconds
    double precision :: h

    h = global_domain_length/N
    hloc = h/(2.d0**(THIS_IMAGE() - 1))
    domain_length = global_domain_length/(2.d0**(THIS_IMAGE() - 1))

    call initial_condition()
    call setup_afac()
    sync all
    call system_clock(t_start, t_rate)
    call residual()
    call l2_norm_residual(l2_norm)
    num_runs = 0
    do while ((l2_norm > epsilon) .and. (num_runs < max_iterations))
        call multigrid_afac
        sync all
        call error_projection

        call residual()

        call l2_norm_residual(l2_norm)
        sync all
        num_runs = num_runs + 1
    end do
    call system_clock(t_end)
    elapsed_seconds = real(t_end - t_start, kind=8)/t_rate
    print '(A, I0, A, F8.4, A)', " [N = ", N, "] Execution Time: ", elapsed_seconds, " seconds"
    open (newunit=unit_num, file="benchmark_results.csv", status="unknown", position="append")
    write (unit_num, '(I0, ",", F12.6)') N, elapsed_seconds
    close (unit_num)
end program benchmark
