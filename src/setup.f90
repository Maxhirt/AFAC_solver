module setup
    implicit none(type, external)
    private

    public :: pi, fourpi, NSLAE, global_domain_length, multigrid_levels, rsp_sphere, rho_sphere, rsp1, &
              rho1, rsp2, rho2, offset_x1, offset_x2, x, b, res, err, &
              rho, bitmask, analytical_solution, domain_length, hloc, grid

    double precision, PARAMETER :: pi = 3.14159265358973238462d0
    double precision, PARAMETER :: fourpi = 4*pi

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

end module setup
