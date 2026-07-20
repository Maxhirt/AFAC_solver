module afac
    use setup, only: N, x, b, domain_length_global, domain_length
    use setup, only: hloc, centerloc_grid, grid
    use boundary, only: init_boundary
    use communication, only: interface_exchange, pack_error, unpack_error
    use multigrid_solver, only: rbgs_smoother, multigrid_setup, copy_afac_to_multigrid, copy_multigrid_to_afac, restriction_operator, prolongation_operator, multigrid_residual
    implicit none(type, external)

    private
    public :: setup_afac

contains
    !> Handles grid setup and boundary condition.
    !!
    !! Calculates the grid cell size for each local grid.
    !! And the center location of the grid as well as initializing
    !! the solution and right-hand side arrays with zeros or the boundary values.
    !!
    subroutine setup_afac()
        implicit none(type, external)
        integer :: i, k, l
        double precision :: h

        h = domain_length_global/N
        hloc = h/(2.d0**(THIS_IMAGE() - 1))
        centerloc_grid = (domain_length/(2.d0**(THIS_IMAGE() - 1)))/2.d0

        call init_boundary()

    end subroutine setup_afac

    !> Calculates the residual of the composite grid.
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine residual()
        implicit none(type, external)
        double precision :: hloc_grid_sq_inv
        integer :: i, k, l

        call interface_exchange()

        if (this_image() > 1) then
            call ghost_cells()
        end if

        if (this_image() < num_images()) then
            call zombie_cells()
        end if

        !$OMP parallel do collapse(3) private(i,k,l,local_laplacian) schedule(static)
        do l = 2, N + 1
            do k = 2, N + 1
                do i = 2, N + 1

                    if (levels%bitmask(i, k, l) == 1) then
                        levels%composite_res(i, k, l) = 0.d0

                    else
                        local_laplacian = (x(i + 1, k, l) + x(i - 1, k, l) + &
                                           x(i, k + 1, l) + x(i, k - 1, l) + &
                                           x(i, k, l + 1) + x(i, k, l - 1) - &
                                           (6.d0*x(i, k, l)))*hloc_grid_sq_inv
                        res(i, k, l) = b(i, k, l) - local_laplacian
                    end if

                end do
            end do
        end do
        !$OMP end parallel do

    end subroutine residual

    !> populates Ghost cells
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine ghost_cells()
        implicit none(type, external)

        double precision :: hloc_coarse, hloc_coarse_inv2, sz, sy, sx, gradient_y, gradient_z, gradient_x
        integer :: i, l, k, l_coarse, k_coarse, i_coarse

        hloc_coarse = (domain_length/N)/(2.d0**(this_image() - 2))
        hloc_coarse_inv2 = 1/(2*hloc_coarse)

        !$OMP PARALLEL PRIVATE(l,k,i,gradient_x,gradient_y,gradient_z,sx,sy,sz,l_coarse,k_coarse,i_coarse)

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        do l = 2, N + 1
            do k = 2, N + 1
                sz = merge(1.d0, -1.d0, mod(l, 2) /= 0)
                sy = merge(1.d0, -1.d0, mod(k, 2) /= 0)

                l_coarse = (l + N/2 + 2)/2 - 1
                k_coarse = (k + N/2 + 2)/2 - 1

                i = 1
                i_coarse = (i + N/2 + 2)/2

                gradient_y = (coarse_comp_x(i_coarse, k_coarse + 1, l_coarse) - coarse_comp_x(i_coarse, k_coarse - 1, l_coarse)) &
                             *hloc_coarse_inv2
                gradient_z = (coarse_comp_x(i_coarse, k_coarse, l_coarse + 1) - coarse_comp_x(i_coarse, k_coarse, l_coarse - 1)) &
                             *hloc_coarse_inv2

                x(i, k, l) = 1.d0/3.d0*x(i + 1, k, l) + &
                             2.d0/3.d0*(coarse_comp_x(i_coarse, k_coarse, l_coarse) + hloc_coarse/2.d0*(sy*gradient_y + &
                                                                                                        sz*gradient_z))

                i = N + 2

                i_coarse = (i + N/2 + 2)/2
                gradient_y = (coarse_comp_x(i_coarse, k_coarse + 1, l_coarse) - coarse_comp_x(i_coarse, k_coarse - 1, l_coarse)) &
                             *hloc_coarse_inv2
                gradient_z = (coarse_comp_x(i_coarse, k_coarse, l_coarse + 1) - coarse_comp_x(i_coarse, k_coarse, l_coarse - 1)) &
                             *hloc_coarse_inv2

                x(i, k, l) = 1.d0/3.d0*x(i - 1, k, l) + &
                             2.d0/3.d0*(coarse_comp_x(i_coarse, k_coarse, l_coarse) + hloc_coarse/2.d0*(sy*gradient_y &
                                                                                                        + sz*gradient_z))

            end do
        end do
        !$OMP END DO

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        do l = 2, N + 1
            do i = 2, N + 1
                sz = merge(1.d0, -1.d0, mod(l, 2) /= 0)
                sx = merge(1.d0, -1.d0, mod(i, 2) /= 0)

                l_coarse = (l + N/2 + 2)/2 - 1
                i_coarse = (i + N/2 + 2)/2

                ! --- k = 1 ---
                k = 1
                k_coarse = (k + N/2 + 2)/2 - 1
                gradient_x = (coarse_comp_x(i_coarse + 1, k_coarse, l_coarse) - coarse_comp_x(i_coarse - 1, k_coarse, l_coarse)) &
                             *hloc_coarse_inv2
                gradient_z = (coarse_comp_x(i_coarse, k_coarse, l_coarse + 1) - coarse_comp_x(i_coarse, k_coarse, l_coarse - 1)) &
                             *hloc_coarse_inv2

                x(i, k, l) = 1.d0/3.d0*x(i, k + 1, l) + &
                             2.d0/3.d0*(coarse_comp_x(i_coarse, k_coarse, l_coarse) + hloc_coarse/2.d0*(sx*gradient_x + &
                                                                                                        sz*gradient_z))

                ! --- k = N+2 ---
                k = N + 2
                k_coarse = (k + N/2 + 2)/2 - 1

                gradient_x = (coarse_comp_x(i_coarse + 1, k_coarse, l_coarse) - coarse_comp_x(i_coarse - 1, k_coarse, l_coarse)) &
                             *hloc_coarse_inv2
                gradient_z = (coarse_comp_x(i_coarse, k_coarse, l_coarse + 1) - coarse_comp_x(i_coarse, k_coarse, l_coarse - 1)) &
                             *hloc_coarse_inv2

                x(i, k, l) = 1.d0/3.d0*x(i, k - 1, l) + &
                             2.d0/3.d0*(coarse_comp_x(i_coarse, k_coarse, l_coarse) + hloc_coarse/2.d0*(sx*gradient_x + &
                                                                                                        sz*gradient_z))
            end do
        end do
        !$OMP END DO

        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
        do k = 2, N + 1
            do i = 2, N + 1
                sy = merge(1.d0, -1.d0, mod(k, 2) /= 0)
                sx = merge(1.d0, -1.d0, mod(i, 2) /= 0)

                k_coarse = (k + N/2 + 2)/2 - 1
                i_coarse = (i + N/2 + 2)/2

                ! --- l = 1 ---
                l = 1
                l_coarse = (l + N/2 + 2)/2 - 1
                gradient_x = (coarse_comp_x(i_coarse + 1, k_coarse, l_coarse) - &
                              coarse_comp_x(i_coarse - 1, k_coarse, l_coarse)) &
                             *hloc_coarse_inv2
                gradient_y = (coarse_comp_x(i_coarse, k_coarse + 1, l_coarse) - &
                              coarse_comp_x(i_coarse, k_coarse - 1, l_coarse)) &
                             *hloc_coarse_inv2

                x(i, k, l) = 1.d0/3.d0*x(i, k, l + 1) + &
                             2.d0/3.d0*(coarse_comp_x(i_coarse, k_coarse, l_coarse) + hloc_coarse/2.d0* &
                                        (sx*gradient_x + sy*gradient_y))

                ! --- l = N+2 ---
                l = N + 2
                l_coarse = (l + N/2 + 2)/2 - 1
                gradient_x = (coarse_comp_x(i_coarse + 1, k_coarse, l_coarse) - coarse_comp_x(i_coarse - 1, k_coarse, l_coarse)) &
                             *hloc_coarse_inv2
                gradient_y = (coarse_comp_x(i_coarse, k_coarse + 1, l_coarse) - coarse_comp_x(i_coarse, k_coarse - 1, l_coarse)) &
                             *hloc_coarse_inv2

                x(i, k, l) = 1.d0/3.d0*x(i, k, l - 1) + &
                             2.d0/3.d0*(coarse_comp_x(i_coarse, k_coarse, l_coarse) + hloc_coarse/2.d0 &
                                        *(sx*gradient_x + sy*gradient_y))
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine ghost_cells

    !> Populates the zombie cells.
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine zombie_cells()
        implicit none(type, external)

        integer :: i, k, l

        !$OMP parallel do collapse(2) private(l,k,i) schedule(static)
        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = N/4 + 2
                x(i, k, l) = 4.d0/3.d0*restricted_interface(i, k, l, 1) - &
                             1.d0/3.d0*x(i - 1, k, l)
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) private(l,k,i) schedule(static)
        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = 3*N/4 + 1
                x(i, k, l) = 4.d0/3.d0*restricted_interface(i, k, l, 2) - &
                             1.d0/3.d0*x(i + 1, k, l)
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) private(l,k,i) schedule(static)
        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = N/4 + 2
                x(i, k, l) = 4.d0/3.d0*restricted_interface(i, k, l, 3) - &
                             1.d0/3.d0*x(i, k - 1, l)
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) private(l,k,i) schedule(static)
        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = 3*N/4 + 1
                x(i, k, l) = 4.d0/3.d0*restricted_interface(i, k, l, 4) - &
                             1.d0/3.d0*x(i, k + 1, l)
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) private(l,k,i,index_zero) schedule(static)
        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = N/4 + 2
                x(i, k, l) = 4.d0/3.d0*restricted_interface(i, k, l, 5) - &
                             1.d0/3.d0*x(i, k, l - 1)
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) private(l,k,i,index_zero) schedule(static)
        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = 3*N/4 + 1
                x(i, k, l) = 4.d0/3.d0*restricted_interface(i, k, l, 6) - &
                             1.d0/3.d0*x(i, k, l + 1)
            end do
        end do
        !$OMP end parallel do

    end subroutine zombie_cells

    !> Calculates the L2 Norm of the Residual
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine l2_norm_residual(L2_error)
        implicit none(type, external)
        double precision, intent(OUT) :: L2_error
        integer :: i, k, l, num_cells
        double precision :: local_l2_error, c_sum

        c_sum = 0.d0
        num_cells = 0
        L2_error = 0.d0

        !$OMP parallel do collapse(3) reduction(+:num_cells, c_sum) private(i,k,l,local_l2_error) schedule(static)
        do l = 2, N + 1
            do k = 2, N + 1
                do i = 2, N + 1
                    if (bitmask(i, k, l) == 1) then
                        local_l2_error = 0.d0
                    else
                        local_l2_error = res(i, k, l)
                        c_sum = c_sum + local_l2_error**2
                        num_cells = num_cells + 1
                    end if
                end do
            end do
        end do
        !$OMP end parallel do
        call co_sum(c_sum)
        call co_sum(num_cells)
        L2_error = sqrt(c_sum/num_cells)

    end subroutine l2_norm_residual

    !> Runs Multigrid cycles for the afac cycle
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine multigrid_afac()
        implicit none(type, external)
        integer :: runs, i
        runs = 0
        call setup_multigrid(grid)

        do while (runs < multigrid_max_iterations)
            call rbgs_smoother(grid(1))
            call multigrid_residual(grid(1))
            do i = 2, multigrid_levels
                call restriction_operator(i - 1)
                call rbgs_smoother(grid(i))
                call multigrid_residual(grid(i))
            end do

            do i = (multigrid_levels - 1), 1, -1
                call prolongation_operator(i + 1)
                call rbgs_smoother(grid(i))
                call multigrid_residual(grid(i))
            end do
            runs = runs + 1
        end do

    end subroutine multigrid_afac

    !> Projects the errors from the different levels to the current one.
    !!
    !! Projects the errors image by image starting from the coarsest level
    !! and moving towards the finest level. Then correcting the approximation with that error.
    subroutine error_projection()
        implicit none(type, external)
        integer :: i
        do i = 2, num_images()
        if (THIS_IMAGE() == i - 1) then
            call pack_error()
        end if
        if (THIS_IMAGE() == i) then
            error_buffer_recv = error_buffer[THIS_IMAGE() - 1] (:)
            call unpack_error()
            call prolongate_error()
        end if
        sync all
        end do

        call error_reconciliation()

    end subroutine error_projection

    !> Prolongates the coarse error to the fine grid.
    !!
    !! Prolongates the coarse error the the fine grid using the error_copy holding only the communicated values.
    subroutine prolongate_error()
        implicit none(type, external)
        integer :: i, k, l, coarse_i, coarse_k, coarse_l, icoarse_neigh, kcoarse_neigh, lcoarse_neigh
        double precision :: prolongated, wx, wy, wz

        wx = 0.75d0
        wy = 0.75d0
        wz = 0.75d0
        !$OMP parallel do collapse(3)  &
        !$OMP private(i,k,l,coarse_i,coarse_k,coarse_l,icoarse_neigh,kcoarse_neigh,lcoarse_neigh,prolongated) &
        !$OMP schedule(static)
        do l = 2, N + 1
        do k = 2, N + 1
        do i = 2, N + 1
            coarse_i = N/4 + 1 + i/2
            coarse_k = N/4 + 1 + k/2
            coarse_l = N/4 + 1 + l/2
            icoarse_neigh = coarse_i + merge(-1, 1, mod(i, 2) == 0)
            kcoarse_neigh = coarse_k + merge(-1, 1, mod(k, 2) == 0)
            lcoarse_neigh = coarse_l + merge(-1, 1, mod(l, 2) == 0)
            prolongated = wx*wy*wz*error_copy(coarse_i, coarse_k, coarse_l) + &
                          (1.d0 - wx)*wy*wz*error_copy(icoarse_neigh, coarse_k, coarse_l) + &
                          wx*(1.d0 - wy)*wz*error_copy(coarse_i, kcoarse_neigh, coarse_l) + &
                          wx*wy*(1.d0 - wz)*error_copy(coarse_i, coarse_k, lcoarse_neigh) + &
                          (1.d0 - wx)*(1.d0 - wy)*wz*error_copy(icoarse_neigh, kcoarse_neigh, coarse_l) + &
                          (1.d0 - wx)*wy*(1.d0 - wz)*error_copy(icoarse_neigh, coarse_k, lcoarse_neigh) + &
                          wx*(1.d0 - wy)*(1.d0 - wz)*error_copy(coarse_i, kcoarse_neigh, lcoarse_neigh) + &
                          (1.d0 - wx)*(1.d0 - wy)*(1.d0 - wz)*error_copy(icoarse_neigh, kcoarse_neigh, lcoarse_neigh)
            err(i, k, l) = err(i, k, l) + prolongated
        end do
        end do
        end do
        !$OMP end parallel do

    end subroutine prolongate_error

    !> Reconciles the errors to the appropriate discretization.
    !!
    !! Corrects the solution by adding the errors to the solution.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine error_reconciliation()
        implicit none(type, external)
        integer :: i, k, l
        !$OMP parallel do collapse(3) private(i,k,l) schedule(static)
        do l = 2, N + 1
        do k = 2, N + 1
        do i = 2, N + 1
            x(i, k, l) = x(i, k, l) + err(i, k, l)
        end do
        end do
        end do
        !$OMP end parallel do

    end subroutine error_reconciliation

end module afac
