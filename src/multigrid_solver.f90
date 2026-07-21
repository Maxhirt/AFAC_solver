module multigrid_solver
    use setup, only: grid, x, b, res, err, grid_level, multigrid_levels, N, domain_length
    implicit none(type, external)
    private
    public :: rbgs_smoother, multigrid_setup, copy_afac_to_multigrid, &
              copy_multigrid_to_afac, restriction_operator, prolongation_operator, multigrid_residual, &
              apply_multigrd_correction

contains

    !> Red-Black Gauss-Seidel solver for each Multigrid level
    !!
    !! Runs the Red-Black Gauss-Seidel solver for each Multigrid level using the level index as an input. The solver
    !! iteratively solves for the unknowns on a grid, alternating between red
    !! and black points, it runs a fixed amount of times.
    !!
    !! @param[in] j The Multigrid level index to solve on.
    subroutine rbgs_smoother(j)
        implicit none(type, external)
        integer, intent(in) :: j
        double precision, parameter :: one_over_six = 1.0d0/6.0d0
        integer :: color, i, k, l, N_local, start_i
        double precision :: hloc_sqr
        hloc_sqr = grid(j)%hloc_grid*grid(j)%hloc_grid

        N_local = grid(j)%N_grid
        !$OMP parallel private(i,k,l,start_i,color)
        do color = 0, 1
            !$OMP do collapse(2) schedule(static)
            do l = 2, N_local + 1
            do k = 2, N_local + 1
                start_i = 2 + iand(k + l + color, 1)
                !$OMP SIMD
                do i = start_i, N_local + 1, 2
                    grid(j)%x(i, k, l) = ((grid(j)%x(i + 1, k, l) + grid(j)%x(i - 1, k, l) + &
                                           grid(j)%x(i, k + 1, l) + grid(j)%x(i, k - 1, l) + grid(j)%x(i, k, l + 1) + &
                                           grid(j)%x(i, k, l - 1)) - grid(j)%b(i, k, l)*hloc_sqr)*one_over_six

                end do
                !$OMP END SIMD
            end do
            end do
            !$OMP end do
        end do
        !$OMP end parallel

    end subroutine rbgs_smoother

    !> Multigrid Setup
    !!
    !! Allocates the multigrid grids with the correct size and initializes them to zero.
    !!
    !! @param[inout] grid Allocatable Multigrid levels.
    subroutine multigrid_setup(grid)
        implicit none(type, external)
        type(grid_level), intent(inout), allocatable :: grid(:)
        integer :: i, N_local

        if (.not. allocated(grid)) then
            print *, "Allocating multigrid grid levels"
            allocate (grid(multigrid_levels))
        end if

        do i = 1, multigrid_levels
            print *, "Allocating grid level ", i
            if (.not. allocated(grid(i)%res)) then
                N_local = N/(2**(i - 1))
                grid(i)%N_grid = N_local
                allocate (grid(i)%res(N_local + 2, N_local + 2, N_local + 2))
                allocate (grid(i)%x(N_local + 2, N_local + 2, N_local + 2))
                allocate (grid(i)%err(N_local + 2, N_local + 2, N_local + 2))
                allocate (grid(i)%b(N_local + 2, N_local + 2, N_local + 2))
                grid(i)%hloc_grid = domain_length/N_local
            end if

            grid(i)%err = 0.d0
            grid(i)%res = 0.d0
            grid(i)%x = 0.d0
            grid(i)%b = 0.d0
        end do

    end subroutine multigrid_setup

    !> Copies the values from the afac level to the multigrid
    !!
    !! Copies the residual from the AFAC to the grid%b and sets the grid%x to zero.
    subroutine copy_afac_to_multigrid()
        implicit none(type, external)
        integer :: i, k, l

        !$OMP parallel do collapse(3) private(i,k,l) schedule(static)
        do l = 1, N + 2
            do k = 1, N + 2
                do i = 1, N + 2
                    grid(1)%b(i, k, l) = res(i, k, l)
                    grid(1)%x(i, k, l) = 0.d0
                end do
            end do
        end do
        !$OMP end parallel do

    end subroutine copy_afac_to_multigrid

    !> Copies the multigrd results to the afac level
    !!
    !! Copies the multigrid error to the AFAC level.
    subroutine copy_multigrid_to_afac()
        implicit none(type, external)
        integer :: i, k, l

        !$OMP parallel do collapse(3) private(i,k,l) schedule(static)
        do l = 2, N + 1
        do k = 2, N + 1
        do i = 2, N + 1
            err(i, k, l) = grid(1)%x(i, k, l)
        end do
        end do
        end do
        !$OMP end parallel do

    end subroutine copy_multigrid_to_afac

    !> Restriction operator multigrid
    !!
    !! Restricts the residual from a finer grid (j) to a coarser grid (j+1).
    !!
    !! @param[in] j index of the finer level
    subroutine restriction_operator(j)
        implicit none(type, external)
        integer, intent(in) :: j
        integer :: i, k, l, N_coarse, if0, if1, kf0, kf1, lf0, lf1

        N_coarse = grid(j + 1)%N_grid

        !$OMP parallel do collapse(3) private(i,k,l,if0,if1,kf0,kf1,lf0,lf1) schedule(static)
        do l = 2, N_coarse + 1
        do k = 2, N_coarse + 1
        do i = 2, N_coarse + 1
            if0 = 2*i - 2
            if1 = 2*i - 1
            kf0 = 2*k - 2
            kf1 = 2*k - 1
            lf0 = 2*l - 2
            lf1 = 2*l - 1
            grid(j + 1)%b(i, k, l) = (grid(j)%res(if0, kf0, lf0) + grid(j)%res(if0, kf0, lf1) + &
                                      grid(j)%res(if0, kf1, lf0) + grid(j)%res(if0, kf1, lf1) + &
                                      grid(j)%res(if1, kf0, lf0) + grid(j)%res(if1, kf0, lf1) + &
                                      grid(j)%res(if1, kf1, lf0) + grid(j)%res(if1, kf1, lf1))/8.d0
        end do
        end do
        end do
        !$OMP end parallel do
    end subroutine restriction_operator

    !> Prolongation operator multigrid
    !!
    !! Prolongates via trilinear interpolation the error of a coarser level (j) to a finer level (j-1)
    !!
    !! @param[in] j index of the coarser level
    subroutine prolongation_operator(j)
        implicit none(type, external)
        integer, intent(in) :: j
        integer i, k, l, N_fine, i_coarse, k_coarse, l_coarse, k_coarse_neigh, i_coarse_neigh, l_coarse_neigh
        double precision :: wx, wy, wz

        wx = 0.75d0
        wy = 0.75d0
        wz = 0.75d0

        N_fine = grid(j - 1)%N_grid

        do l = 2, N_fine + 1
            do k = 2, N_fine + 1
                do i = 2, N_fine + 1
                    i_coarse = 2*i - 1
                    k_coarse = 2*k - 1
                    l_coarse = 2*l - 1

                    i_coarse_neigh = i_coarse + merge(-1, 1, mod(i, 2) == 0)
                    k_coarse_neigh = k_coarse + merge(-1, 1, mod(k, 2) == 0)
                    l_coarse_neigh = l_coarse + merge(-1, 1, mod(l, 2) == 0)

                    grid(j - 1)%err(i, k, l) = wx*wy*wz*grid(j)%x(i_coarse, k_coarse, l_coarse) + &
                                               (1.d0 - wx)*wy*wz*grid(j)%x(i_coarse_neigh, k_coarse, l_coarse) + &
                                               wx*(1.d0 - wy)*wz*grid(j)%x(i_coarse, k_coarse_neigh, l_coarse) + &
                                               wx*wy*(1.d0 - wz)*grid(j)%x(i_coarse, k_coarse, l_coarse_neigh) + &
                                               (1.d0 - wx)*(1.d0 - wy)*wz*grid(j)%x(i_coarse_neigh, k_coarse_neigh, l_coarse) + &
                                               (1.d0 - wx)*wy*(1.d0 - wz)*grid(j)%x(i_coarse_neigh, k_coarse, l_coarse_neigh) + &
                                               wx*(1.d0 - wy)*(1.d0 - wz)*grid(j)%x(i_coarse, k_coarse_neigh, l_coarse_neigh) + &
                                       (1.d0 - wx)*(1.d0 - wy)*(1.d0 - wz)*grid(j)%x(i_coarse_neigh, k_coarse_neigh, l_coarse_neigh)

                end do
            end do
        end do
    end subroutine prolongation_operator

    !> Calculates the residual on a multigrid level
    !!
    !! Calculates the residual for the gravitational potential.
    !!
    !! @param[in] j multigrid level index
    subroutine multigrid_residual(j)
        implicit none(type, external)
        integer, intent(in) :: j
        integer :: i, k, l, N_local
        double precision :: hloc_sqr_inv

        N_local = grid(j)%N_grid
        hloc_sqr_inv = 1.d0/(grid(j)%hloc_grid*grid(j)%hloc_grid)
        !$OMP parallel do collapse(3) private(i,k,l) schedule(static)
        do l = 2, N_local
        do k = 2, N_local
        do i = 2, N_local
            grid(j)%res(i, k, l) = grid(j)%b(i, k, l) - &
                              (grid(j)%x(i - 1, k, l) + grid(j)%x(i + 1, k, l) + grid(j)%x(i, k + 1, l) + grid(j)%x(i, k - 1, l) + &
                                    grid(j)%x(i, k, l + 1) + grid(j)%x(i, k, l - 1) - (6.d0*grid(j)%x(i, k, l)))*hloc_sqr_inv
        end do
        end do
        end do
        !$OMP end parallel do

    end subroutine multigrid_residual

    !> Applies the correction
    !!
    !! Applies the correction of the coarser grid to the finer grid.
    !!
    !! @param[in] j multigrid level index
    subroutine apply_multigrd_correction(j)
        implicit none(type, external)
        integer, intent(in) :: j
        integer :: i, k, l, N_local

        N_local = grid(j)%N_grid

        !$OMP parallel do collapse(3) private(i,k,l) schedule(static)
        do l = 2, N_local + 1
            do k = 2, N_local + 1
                do i = 2, N_local + 1
                    grid(j)%x(i, k, l) = grid(j)%x(i, k, l) + grid(j)%err(i, k, l)
                end do
            end do
        end do
        !$OMP end parallel do

    end subroutine apply_multigrd_correction
end module multigrid_solver
