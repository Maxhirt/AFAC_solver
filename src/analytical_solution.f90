module analytical_solution

    use setup, only: N, fourpi, G, rho, boundary_type, ana_solution, relative_error, x
    use boundary, only: boundary_value

    private
    public :: analytical_solution_calculator
contains

!> Calculates the analytical solution
!!
!! calculates the analytical solution for the given boundary conditions.
    subroutine analytical_solution_calculator()
        implicit none(type, external)
        integer :: i, k, l

        do l = 1, N + 2
        do k = 1, N + 2
        do i = 1, N + 2
            ana_solution(i, k, l) = boundary_value(i, k, l)
            relative_error(i, k, l) = abs(x(i, k, l) - ana_solution(i, k, l))/abs(ana_solution(i, k, l))
        end do
        end do
        end do

    end subroutine analytical_solution_calculator

end module analytical_solution
