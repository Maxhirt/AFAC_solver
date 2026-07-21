module communication
    use setup, only: x, b, err, res, N, restricted_interface, restricted_interface_buffer, &
                     coarse_cell_buffer, restricted_interface_buffer_recv, coarse_comp_x, coarse_cell_buffer_recv, &
                     error_buffer, error_buffer_recv, error_copy
    implicit none(type, external)
    PRIVATE
    public :: interface_exchange, pack_error, unpack_error

contains

    !> Does the data exchange between images for the interfaces.
    !!
    !! Does fine and coarse interface data exchange modeled after the momentum conserving interface handling.
    !!
    subroutine interface_exchange
        implicit none(type, external)

        if (this_image() > 1) then
            call restrict_layer
            call pack_fine
        end if
        if (THIS_IMAGE() < num_images()) then
            call pack_coarse
        end if
        sync all
        if (this_image() > 1) then
            coarse_cell_buffer_recv(:) = coarse_cell_buffer[this_image() - 1] (:)
            call unpack_coarse()
        end if
        if (this_image() < num_images()) then
            restricted_interface_buffer_recv(:) = restricted_interface_buffer[this_image() + 1] (:)
            call unpack_fine
        end if
        sync all

    end subroutine interface_exchange

    !> Restricts the first active fine layer (plus from ATHENA paper)
    !!
    !! Detailed description of what it does.
    subroutine restrict_layer()
        implicit none(type, external)
        double precision :: restricted
        integer :: i, k, l, i_f0, i_f1, k_f0, k_f1, l_f0, l_f1

        !$OMP parallel do collapse(2) schedule(static) &
        !$OMP private(i,k,l,restricted, i_f0,k_f0,k_f1,l_f0,l_f1)
        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = N/4 + 2
                i_f0 = 2*i - N/2 - 2
                k_f0 = 2*k - N/2 - 2
                k_f1 = 2*k - N/2 - 1
                l_f0 = 2*l - N/2 - 2
                l_f1 = 2*l - N/2 - 1

                restricted = (x(i_f0, k_f0, l_f0) + x(i_f0, k_f0, l_f1) + &
                              x(i_f0, k_f1, l_f0) + x(i_f0, k_f1, l_f1))/4.d0

                restricted_interface(i, k, l, 1) = restricted
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) schedule(static) &
        !$OMP private(i,k,l,restricted, i_f0,k_f0,k_f1,l_f0,l_f1)
        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = 3*N/4 + 1
                i_f0 = 2*i - N/2 - 1
                k_f0 = 2*k - N/2 - 2
                k_f1 = 2*k - N/2 - 1
                l_f0 = 2*l - N/2 - 2
                l_f1 = 2*l - N/2 - 1

                restricted = (x(i_f0, k_f0, l_f0) + x(i_f0, k_f0, l_f1) + &
                              x(i_f0, k_f1, l_f0) + x(i_f0, k_f1, l_f1))/4.d0

                restricted_interface(i, k, l, 2) = restricted
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) schedule(static) &
        !$OMP private(i,k,l,restricted, i_f0,k_f0,i_f1,l_f0,l_f1)
        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = N/4 + 2
                i_f0 = 2*i - N/2 - 2
                i_f1 = 2*i - N/2 - 1
                k_f0 = 2*k - N/2 - 2
                l_f0 = 2*l - N/2 - 2
                l_f1 = 2*l - N/2 - 1

                restricted = (x(i_f0, k_f0, l_f0) + x(i_f0, k_f0, l_f1) + &
                              x(i_f1, k_f0, l_f0) + x(i_f1, k_f0, l_f1))/4.d0
                restricted_interface(i, k, l, 3) = restricted
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) schedule(static) &
        !$OMP private(i,k,l,restricted, i_f0,k_f0,i_f1,l_f0,l_f1)
        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = 3*N/4 + 1
                i_f0 = 2*i - N/2 - 2
                i_f1 = 2*i - N/2 - 1
                k_f0 = 2*k - N/2 - 1
                l_f0 = 2*l - N/2 - 2
                l_f1 = 2*l - N/2 - 1

                restricted = (x(i_f0, k_f0, l_f0) + x(i_f0, k_f0, l_f1) + &
                              x(i_f1, k_f0, l_f0) + x(i_f1, k_f0, l_f1))/4.d0
                restricted_interface(i, k, l, 4) = restricted
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) schedule(static) &
        !$OMP private(i,k,l,restricted, i_f0,k_f0,k_f1,l_f0,i_f1)
        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = N/4 + 2
                i_f0 = 2*i - N/2 - 2
                i_f1 = 2*i - N/2 - 1
                k_f0 = 2*k - N/2 - 2
                l_f0 = 2*l - N/2 - 2
                k_f1 = 2*k - N/2 - 1

                restricted = (x(i_f0, k_f0, l_f0) + x(i_f0, k_f1, l_f0) + &
                              x(i_f1, k_f0, l_f0) + x(i_f1, k_f1, l_f0))/4.d0

                restricted_interface(i, k, l, 5) = restricted
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) schedule(static) &
        !$OMP private(i,k,l,restricted, i_f0,k_f0,k_f1,l_f0,i_f1)
        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = 3*N/4 + 1
                i_f0 = 2*i - N/2 - 2
                i_f1 = 2*i - N/2 - 1
                k_f0 = 2*k - N/2 - 2
                l_f0 = 2*l - N/2 - 1
                k_f1 = 2*k - N/2 - 1

                restricted = (x(i_f0, k_f0, l_f0) + x(i_f0, k_f1, l_f0) + &
                              x(i_f1, k_f0, l_f0) + x(i_f1, k_f1, l_f0))/4.d0
                restricted_interface(i, k, l, 6) = restricted
            end do
        end do
        !$OMP end parallel do

    end subroutine restrict_layer

    !> Packs the restricted interface array into a 1D coarray.
    !!
    !! Packages restricted_interface into the 1D coarray
    !!for communication with the coarser grid.
    subroutine pack_fine
        implicit none(type, external)

        integer :: i, k, l, idx

        idx = 1
        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = N/4 + 2
                restricted_interface_buffer(idx) = restricted_interface(i, k, l, 1)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = 3*N/4 + 1
                restricted_interface_buffer(idx) = restricted_interface(i, k, l, 2)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = N/4 + 2
                restricted_interface_buffer(idx) = restricted_interface(i, k, l, 3)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = 3*N/4 + 1
                restricted_interface_buffer(idx) = restricted_interface(i, k, l, 4)
                idx = idx + 1
            end do
        end do

        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = N/4 + 2
                restricted_interface_buffer(idx) = restricted_interface(i, k, l, 5)
                idx = idx + 1
            end do
        end do

        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = 3*N/4 + 1
                restricted_interface_buffer(idx) = restricted_interface(i, k, l, 6)
                idx = idx + 1
            end do
        end do

    end subroutine pack_fine

    !> Packs the coarse cell into 1D coarray.
    !!
    !! Uses the coarse cell buffer to communcate with the finer grid.
    subroutine pack_coarse
        implicit none(type, external)
        integer :: i, k, l, idx

        idx = 1

        do l = N/4 + 1, 3*N/4 + 2
            do k = N/4 + 1, 3*N/4 + 2
                i = N/4 + 1
                coarse_cell_buffer(idx) = x(i, k, l)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 1, 3*N/4 + 2
            do k = N/4 + 1, 3*N/4 + 2
                i = 3*N/4 + 2
                coarse_cell_buffer(idx) = x(i, k, l)
                idx = idx + 1

            end do
        end do

        do l = N/4 + 1, 3*N/4 + 2
            k = N/4 + 1
            do i = N/4 + 1, 3*N/4 + 2
                coarse_cell_buffer(idx) = x(i, k, l)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 1, 3*N/4 + 2
            k = 3*N/4 + 2
            do i = N/4 + 1, 3*N/4 + 2
                coarse_cell_buffer(idx) = x(i, k, l)
                idx = idx + 1

            end do
        end do

        l = N/4 + 1
        do k = N/4 + 1, 3*N/4 + 2
            do i = N/4 + 1, 3*N/4 + 2
                coarse_cell_buffer(idx) = x(i, k, l)
                idx = idx + 1
            end do
        end do

        l = 3*N/4 + 2
        do k = N/4 + 1, 3*N/4 + 2
            do i = N/4 + 1, 3*N/4 + 2
                coarse_cell_buffer(idx) = x(i, k, l)
                idx = idx + 1
            end do
        end do

    end subroutine pack_coarse

    !> Stores the received interpolated fine face values into the restricted interface array.
    !!
    !! Stores the received interpolated fine face values from the buffer into the restricted interface array.
    subroutine unpack_fine
        implicit none(type, external)
        integer :: i, k, l, idx
        restricted_interface = 0.d0

        idx = 1
        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = N/4 + 2
                restricted_interface(i, k, l, 1) = restricted_interface_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 2, 3*N/4 + 1
            do k = N/4 + 2, 3*N/4 + 1
                i = 3*N/4 + 1
                restricted_interface(i, k, l, 2) = restricted_interface_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = N/4 + 2
                restricted_interface(i, k, l, 3) = restricted_interface_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                k = 3*N/4 + 1
                restricted_interface(i, k, l, 4) = restricted_interface_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = N/4 + 2
                restricted_interface(i, k, l, 5) = restricted_interface_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do k = N/4 + 2, 3*N/4 + 1
            do i = N/4 + 2, 3*N/4 + 1
                l = 3*N/4 + 1
                restricted_interface(i, k, l, 6) = restricted_interface_buffer_recv(idx)
                idx = idx + 1
            end do
        end do
    end subroutine unpack_fine

    !> unpacks the receiving buffer
    !!
    !! Unpacks the receiving buffer into the coarse component x.
    subroutine unpack_coarse
        implicit none(type, external)
        integer :: i, k, l, idx

        idx = 1

        do l = N/4 + 1, 3*N/4 + 2
            do k = N/4 + 1, 3*N/4 + 2
                i = N/4 + 1
                coarse_comp_x(i, k, l) = coarse_cell_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 1, 3*N/4 + 2
            do k = N/4 + 1, 3*N/4 + 2
                i = 3*N/4 + 2
                coarse_comp_x(i, k, l) = coarse_cell_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 1, 3*N/4 + 2
            k = N/4 + 1
            do i = N/4 + 1, 3*N/4 + 2
                coarse_comp_x(i, k, l) = coarse_cell_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        do l = N/4 + 1, 3*N/4 + 2
            k = 3*N/4 + 2
            do i = N/4 + 1, 3*N/4 + 2
                coarse_comp_x(i, k, l) = coarse_cell_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        l = N/4 + 1
        do k = N/4 + 1, 3*N/4 + 2
            do i = N/4 + 1, 3*N/4 + 2
                coarse_comp_x(i, k, l) = coarse_cell_buffer_recv(idx)
                idx = idx + 1
            end do
        end do

        l = 3*N/4 + 2
        do k = N/4 + 1, 3*N/4 + 2
            do i = N/4 + 1, 3*N/4 + 2
                coarse_comp_x(i, k, l) = coarse_cell_buffer_recv(idx)
                idx = idx + 1
            end do
        end do
    end subroutine unpack_coarse

    !> Packs the error into a 1D array
    !!
    !! Packs the error for transfer to a finer grid of the center part covered by a finer grid.
    subroutine pack_error()
        implicit none(type, external)
        integer :: i, k, l, idx, N_start, N_end, L_size

        N_start = N/4 + 1
        N_end = 3*N/4 + 2
        L_size = N_end - N_start + 1
        !$OMP parallel do collapse(3) private(l,k,i,idx)
        do l = N_start, N_end
            do k = N_start, N_end
                do i = N_start, N_end
                    idx = 1 + (i - N_start) + &
                          (k - N_start)*L_size + &
                          (l - N_start)*L_size*L_size
                    error_buffer(idx) = err(i, k, l)
                end do
            end do
        end do
        !$OMP end parallel do
    end subroutine pack_error

    !> Unpacks the error from the coarser grid.
    !!
    !! Unpacks the error from the coarse grid packaged as a 1D array.
    !! This is done via a parallelized indexing.
    subroutine unpack_error()
        implicit none(type, external)
        integer :: i, k, l, idx, N_start, N_end, L_size

        N_start = N/4 + 1
        N_end = 3*N/4 + 2
        L_size = N_end - N_start + 1
        !$OMP parallel do collapse(3) private(l,k,i,idx)
        do l = N_start, N_end
            do k = N_start, N_end
                do i = N_start, N_end
                    idx = 1 + (i - N_start) + &
                          (k - N_start)*L_size + &
                          (l - N_start)*L_size*L_size
                    error_copy(i, k, l) = error_buffer_recv(idx)
                end do
            end do
        end do
        !$OMP end parallel do
    end subroutine unpack_error

end module communication
