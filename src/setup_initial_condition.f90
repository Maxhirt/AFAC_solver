module setup_initial_condition
    use setup, only: rho, rsp1, rsp2, domain_length, boundary_type, hloc, N, rsp_sphere, rho_sphere, rsp1, &
                     rho1, rsp2, rho2, offset_x1, offset_x2, semi_x, semi_z, r_c
    implicit none(type, external)
    private
    public :: initial_condition

contains

    !> Puts in the initial_condition for the density
    !!
    !! Sets the density according to the problem setup.
    subroutine initial_condition()
        implicit none(type, external)
        double precision :: xpos, ypos, zpos, radius, radius1, radius2
        integer :: i, k, l

        select case (boundary_type)!sphere
        case (0)
            do l = 1, N + 2
                do k = 1, N + 2
                    do i = 1, N + 2
                        xpos = (i - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        ypos = (k - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        zpos = (l - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        radius = sqrt((xpos)*(xpos) + ypos*ypos + zpos*zpos)
                        rho(i, k, l) = 1.d-7
                        if (radius <= rsp_sphere) then
                            rho(i, k, l) = rho_sphere
                        end if
                    end do
                end do
            end do

        case (1)!binary

            do l = 1, N + 2
                do k = 1, N + 2
                    do i = 1, N + 2
                        xpos = (i - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        ypos = (k - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        zpos = (l - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        radius1 = sqrt((xpos - offset_x1)*(xpos - offset_x1) + ypos*ypos + zpos*zpos)
                        radius2 = sqrt((xpos - offset_x2)*(xpos - offset_x2) + ypos*ypos + zpos*zpos)
                        rho(i, k, l) = 1.d-7
                        if (radius1 <= rsp1) then
                            rho(i, k, l) = rho1
                        end if
                        if (radius2 <= rsp2) then
                            rho(i, k, l) = rho2
                        end if

                    end do
                end do
            end do

        case (2)!condensed sphere
            do l = 1, N + 2
                do k = 1, N + 2
                    do i = 1, N + 2
                        xpos = (i - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        ypos = (k - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        zpos = (l - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        radius = sqrt((xpos)*(xpos) + ypos*ypos + zpos*zpos)
                        rho(i, k, l) = 1.d-7
                        if (radius <= rsp_sphere) then
                            rho(i, k, l) = rho_sphere/(1 + (radius/r_c)**2)
                        end if
                    end do
                end do
            end do

        case (3)
            do l = 1, N + 2
                do k = 1, N + 2
                    do i = 1, N + 2
                        xpos = (i - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        ypos = (k - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        zpos = (l - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                        radius = sqrt((xpos)*(xpos) + ypos*ypos + zpos*zpos)
                        rho(i, k, l) = 1.d-7
                        if (radius**2/semi_x**2 + zpos**2/semi_z**2 <= 1.d0) then
                            rho(i, k, l) = rho_sphere
                        end if
                    end do
                end do
            end do

        case default
            rho = 1.d-7

        end select
    end subroutine initial_condition
end module setup_initial_condition
