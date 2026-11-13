program coarray_probe
    use iso_fortran_env, only: output_unit
    implicit none

    write (output_unit, '(A,I0,A,I0)') 'COARRAY_PROBE_STARTED image=', this_image(), ' images=', num_images()
    flush (output_unit)

    sync all

    if (this_image() == 1) then
        write (output_unit, '(A,I0)') 'COARRAY_PROBE_ALL_IMAGES_REACHED_SYNC images=', num_images()
        flush (output_unit)
    end if

    sync all

    write (output_unit, '(A,I0)') 'COARRAY_PROBE_FINISHED image=', this_image()
    flush (output_unit)
end program coarray_probe
