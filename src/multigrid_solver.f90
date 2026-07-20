module multigrid_solver
    use setup, only: grid_level
    implicit none(type, external)
    private

contains
    !> Red-Black Gauss-Seidel solver for each Multigrid level
    !!
    !! Runs the Red-Black Gauss-Seidel solver for each Multigrid level using the level as an input. The solver
    !! iteratively solves for the unknowns on a grid, alternating between red
    !! and black points, it runs a fixed amount of times.
    !!
    !! @param[in] grid The Multigrid level to solve on.
    subroutine rbgs_smoother(grid)
        implicit none(type, external)
        type(grid_level), intent(in) :: grid(:)
        integer :: i, j, k

    end subroutine rbgs_smoother

    !> Multigrid Setup
    !!
    !! Copies values from the level to the grid and allocates the grid if it did not alrady exist.
    !!
    !! @param[inout] grid Allocatable Multigrid levels.
    subroutine multgrid_setup(grid)
        implicit none(type, external)
        type(grid_level), intent(inout) :: grid(:)
        integer :: i, k, l

    end subroutine multigrid_setup
end module multigrid_solver
