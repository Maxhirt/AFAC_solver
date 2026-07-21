module boundary
    use setup, only: N, x, b, boundary_type, &
                     G, rho, rsp_sphere, domain_length, &
                     hloc, pi, fourpi, rho_sphere, &
                     offset_x1, offset_x2, rsp1, &
                     rho1, rsp2, rho2, semi_x, semi_z, r_c

    implicit none(type, external)

    private
    public :: init_boundary, boundary_value

contains

    DOUBLE PRECISION function boundary_value(i, k, l)
        implicit none(type, external)

        integer, intent(in) :: i, k, l
        double precision :: xpos, ypos, zpos, radius, radius1, radius2, grav1, grav2
        double precision :: prom1, prom2, lambda, I1, AA1, AA2, AA3, E, sqroot

        xpos = (i - 1)*hloc - domain_length/2.d0 - hloc/2.d0
        ypos = (k - 1)*hloc - domain_length/2.d0 - hloc/2.d0
        zpos = (l - 1)*hloc - domain_length/2.d0 - hloc/2.d0
        radius = sqrt(xpos**2 + ypos**2 + zpos**2)
        radius1 = sqrt((xpos - offset_x1)*(xpos - offset_x1) + ypos*ypos + zpos*zpos)
        radius2 = sqrt((xpos - offset_x2)*(xpos - offset_x2) + ypos*ypos + zpos*zpos)

        select case (boundary_type)

        case (0) !homogenous sphere
            if (radius <= rsp_sphere) then
                boundary_value = -2*pi*G*rho_sphere*(rsp_sphere**2 - radius**2/3.d0)

            else
                boundary_value = -4.d0*pi*G*rho_sphere*rsp_sphere**3/(3.d0*radius)
            end if

        case (1) !binary
            if (radius1 <= rsp1) then
                grav1 = -2*pi*G*rho1*(rsp1**2 - radius1**2/3.d0)
            else
                grav1 = -4*pi*G*rho1*rsp1**3/(3.d0*radius1)
            end if
            if (radius2 <= rsp2) then
                grav2 = -2*pi*G*rho2*(rsp2**2 - radius2**2/3.d0)
            else
                grav2 = -4*pi*G*rho2*rsp2**3/(3.d0*radius2)
            end if
            boundary_value = grav1 + grav2

        case (2) ! condensed sphere

            if (radius <= 1.d-12) then
                boundary_value = 0.d0
            else if (radius <= rsp_sphere) then
                boundary_value = fourpi*G*rho_sphere*r_c**2* &
                                 (atan(radius/r_c)/(radius/r_c) + log((1.d0 + &
                                                                       (radius/r_c)**2)/(1.d0 + (rsp_sphere/r_c)**2))/2.d0 - 1.d0)

            else
                boundary_value = -fourpi*G*rho_sphere*r_c**3*(rsp_sphere/r_c - atan(rsp_sphere/r_c))/radius

            end if

        case (3) ! homogenous ellipsoid

            E = sqrt(1.d0 - (semi_z/semi_x)**2)
            AA1 = (sqrt(1.d0 - E**2)/E**3)*asin(E) - (1.d0 - E**2)/E**2
            AA3 = (2.d0/E**2) - 2.d0*(sqrt(1 - E**2)/E**3)*asin(E)

            if (radius**2/semi_x**2 + zpos**2/semi_z**2 <= 1.d0) then
                boundary_value = -pi*G*rho_sphere*(semi_x**2*AA1 + semi_x**2*AA1 + &
                                                   semi_z**2*AA3 - AA1*radius**2 - AA3*zpos**2)
            else

                prom1 = radius**2 + zpos**2 - semi_x**2 - semi_z**2

                prom2 = 4.d0*(semi_x**2*semi_z**2 - radius**2*semi_z**2 - zpos**2*semi_x**2)

                lambda = (prom1 + sqrt(max(0.d0, prom1**2 - prom2)))/2.d0

                I1 = pi/sqrt(semi_x**2 - semi_z**2) - &
                     2.d0/sqrt(semi_x**2 - semi_z**2)*atan(sqrt((semi_z**2 + lambda)/(semi_x**2 - semi_z**2)))
                sqroot = sqrt(semi_z**2 + lambda)

                boundary_value = -pi*G*rho_sphere*semi_x**2*semi_z*((1.d0 + radius**2/(2.d0*(semi_z**2 - semi_x**2)) - &
                                                                     zpos**2/(semi_z**2 - semi_x**2))*I1 - &
                                                                 radius**2*sqroot/((semi_z**2 - semi_x**2)*(semi_x**2 + lambda)) - &
                                                                    zpos**2*(2.d0/((semi_x**2 + lambda)*sqroot) - &
                                                                        2.d0*sqroot/((semi_z**2 - semi_x**2)*(semi_x**2 + lambda))))

            end if

        case default
            boundary_value = 0.d0

        end select

    end function boundary_value

    !> Sets the boundary values on the outermost level.
    !!
    !! Sets the inital values for x and b as well as setting
    !! the boundary values for the outermost level.
    subroutine init_boundary()

        integer :: i, k, l

        if (THIS_IMAGE() == 1) then
            do l = 1, N + 2
            do k = 1, N + 2
            do i = 1, N + 2
            if ((i == 1) .or. (i == N + 2) .or. (k == 1) .or. (k == N + 2) .or. (l == 1) .or. (l == N + 2)) then
                x(i, k, l) = boundary_value(i, k, l)
            else
                x(i, k, l) = 0.d0
            end if
            b(i, k, l) = fourpi*G*rho(i, k, l)
            end do
            end do
            end do

        else
            do l = 1, N + 2
            do k = 1, N + 2
            do i = 1, N + 2
                x(i, k, l) = 0.d0
                b(i, k, l) = fourpi*G*rho(i, k, l)
            end do
            end do
            end do
        end if

    end subroutine init_boundary

end module boundary
