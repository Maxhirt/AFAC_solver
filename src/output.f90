module output
    use setup, only: N, relative_error, x, ana_solution, hloc, domain_length

    private
    public :: output_results_global

contains

    !> Outputs the results and the relative_error
    !!
    !! Outputs results as well as the relative error.
    !!
    !! @param[in] run The run number to use in the filename.
    subroutine output_results_global(run)
        implicit none(type, external)
        integer, intent(in) :: run
        integer :: i, k, l
        integer :: mid_start, mid_end, u_comb, img
        double precision :: xpos, ypos, zpos
        character(len=32) :: filename

        write (filename, "(A,I0,A)") "global_it_", run, ".dat"
        if (this_image() == 1) then
            open (newunit=u_comb, file=filename, status="replace", action="write")
            write (u_comb, '(A)') '# Level  Xpos           Ypos           Zpos           Numerical      Analytical     Rel_Error'
            close (u_comb)
        end if
        sync all
        mid_start = (N/4) + 2
        mid_end = (3*N/4) + 1

        do img = 1, num_images() - 1
            if (this_image() == img) then
                open (newunit=u_comb, file=trim(filename), position="append", status="old", action="write")

                do l = 2, N + 1
                do k = 2, N + 1
                do i = 2, N + 1

                    if ((i >= mid_start .and. i <= mid_end) .and. &
                        (l >= mid_start .and. l <= mid_end) .and. &
                        (k >= mid_start .and. k <= mid_end)) cycle

                    xpos = (i - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                    ypos = (k - 1)*hloc - domain_length/2.d0 - hloc/2.d0
                    zpos = (l - 1)*hloc - domain_length/2.d0 - hloc/2.d0
           write (u_comb, '(I5, 2X, 6(E14.7, 2X), I2)') this_image(), xpos, ypos, zpos, x(i, k, l), ana_solution(i, k, l), relative_error(i, k, l)

                end do
                end do
                end do
                close (u_comb)
            end if
            sync all
        end do

        if (this_image() == num_images()) then
            open (newunit=u_comb, file=trim(filename), position="append", status="old", action="write")

            do l = 2, N + 1
            do k = 2, N + 1
            do i = 2, N + 1

                xpos = (i - 1)*hloc - (domain_length/(2.d0**(THIS_IMAGE() - 1)))/2.d0 - hloc/2.d0
                ypos = (k - 1)*hloc - (domain_length/(2.d0**(THIS_IMAGE() - 1)))/2.d0 - hloc/2.d0
                zpos = (l - 1)*hloc - (domain_length/(2.d0**(THIS_IMAGE() - 1)))/2.d0 - hloc/2.d0
           write (u_comb, '(I5, 2X, 6(E14.7, 2X), I2)') this_image(), xpos, ypos, zpos, x(i, k, l), ana_solution(i, k, l), relative_error(i, k, l)

            end do
            end do
            end do
            close (u_comb)
        end if

    end subroutine output_results_global

end module output
