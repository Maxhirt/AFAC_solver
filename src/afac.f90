module afac
    use setup, only: N, x, b, domain_length_global, domain_length
    use setup, only: hloc, centerloc_grid
    use boundary, only: set_boundary
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

        call set_boundary()

    end subroutine setup_afac

    !> Calculates the residual of the composite grid.
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine residual()
        implicit none(type, external)

    end subroutine residual

    !> Exchanges and populates Ghost and Zombie cells
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine ghost_zombie_cells()
        implicit none(type, external)

    end subroutine ghost_zombie_cells

    !> Calculates the L2 Norm of the Residual
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine l2_norm_residual()
        implicit none(type, external)

    end subroutine l2_norm_residual

    !> Runs Multigrid cycles for the afac cycle
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine multigrid_afac()
        implicit none(type, external)

    end subroutine multigrid_afac

    !> Projects the errors from the different levels to the current one.
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine error_projection()
        implicit none(type, external)

    end subroutine error_projection

    !> Reconciles the errors to the appropriate discretization.
    !!
    !! Detailed description of what it does.
    !!
    !! @param[in] param1 Description of param1.
    !! @param[out] param2 Description of param2.
    subroutine error_reconciliation()
        implicit none(type, external)

    end subroutine error_reconciliation

end module afac
