module setup
    implicit none(type, external)
    private

    public :: pi, fourpi, NSLAE, global_domain_length, multigrid_levels, rsp_sphere, rho_sphere, rsp1, &
              rho1, rsp2, rho2, offset_x1, offset_x2, x, b, res, err, &
              rho, bitmask, analytical_solution, domain_length, hloc, grid, G, boundary_type, semi_x, &
              semi_z, restricted_interface, restricted_interface_buffer, coarse_cell_buffer, &
              restricted_interface_buffer_recv, error_copy, error_buffer, error_buffer_recv

    double precision, PARAMETER :: pi = 3.14159265358973238462d0
    double precision, PARAMETER :: fourpi = 4*pi
    integer, PARAMETER :: boundary_type = 0
    double precision, PARAMETER :: G = 1.d0

#ifndef GRID_N
#define GRID_N 64
#endif
    ! Number of grid points and domain definition
    integer, PARAMETER :: N = Grid_N
    integer, PARAMETER :: NSLAE = (N + 2)*(N + 2)*(N + 2)
    integer, PARAMETER :: global_domain_length = 1.d0
    integer, PARAMETER :: multigrid_levels = 3

    ! Constants
    double precision, PARAMETER :: rsp_sphere = 0.25d0
    double precision, PARAMETER :: rho_sphere = 1.d0
    double precision, PARAMETER :: rsp1 = 1.d0
    double precision, PARAMETER :: rho1 = 1.d0
    double precision, PARAMETER :: rsp2 = 1.d0
    double precision, PARAMETER :: rho2 = 1.d0
    double precision, PARAMETER :: offset_x1 = 1.d0
    double precision, PARAMETER :: offset_x2 = 1.d0
    double precision, PARAMETER :: semi_x = 1.d0
    double precision, PARAMETER :: semi_z = 0.5d0

    ! Arrays
    double precision :: x(N + 2, N + 2, N + 2)
    double precision :: b(N + 2, N + 2, N + 2)
    double precision :: res(N + 2, N + 2, N + 2)
    double precision :: err(N + 2, N + 2, N + 2)
    double precision :: rho(N + 2, N + 2, N + 2)
    integer :: bitmask(N + 2, N + 2, N + 2)
    double precision :: analytical_solution(N + 2, N + 2, N + 2)

    ! descriptive variables
    double precision :: domain_length
    double precision :: hloc

    ! Communication arrays

    ! Multigrid derived typ
    type :: grid_level
        double precision, Allocatable :: res(:, :, :)
        double precision, Allocatable :: err(:, :, :)
        double precision, Allocatable :: x(:, :, :)
        double precision, Allocatable :: b(:, :, :)
        integer :: N_grid
        double precision :: holc_grid
        double precision :: centerloc_grid
    end type grid_level

    type(grid_level), Allocatable :: grid(:)

    ! communication variables
    double precision :: restricted_interface((N + 2), (N + 2), (N + 2), 6)
    DOUBLE PRECISION :: restricted_interface_buffer(6*(N/2)*(N/2)) [*]
    DOUBLE PRECISION :: coarse_cell_buffer(6*(N/2 + 2)*(N/2 + 2)) [*]
    double precision :: restricted_interface_buffer_recv(6*(N/2)*(N/2))
    double precision :: coarse_comp_x(N + 2, N + 2, N + 2)
    DOUBLE PRECISION :: coarse_cell_buffer_recv(6*(N/2 + 2)*(N/2 + 2))
    double precision :: err_copy(N/4 + 1:3*N/4 + 2, N/4 + 1:3*N/4 + 2, N/4 + 1:3*N/4 + 2)
    double precision :: error_buffer((N/2 + 2)*(N/2 + 2)*(N/2 + 2)) [*]
    double precision :: error_buffer_recv((N/2 + 2)*(N/2 + 2)*(N/2 + 2))

end module setup
