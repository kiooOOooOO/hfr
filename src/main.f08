program hfr_main
    use hf
    implicit none

    type(situation) :: s
    real(8) :: energy

    integer :: i

    write (*,*) "Hello World!"

    s%num_nucleuses = 2
    s%num_basis = 2
    allocate(s%basis_functions(2))
    allocate(s%nucleuses(2))
    do i=1,2
        s%basis_functions(i) = stong_new(3, &
            (/ &
                pgto_new(3.425250914d0, 0d0,0d0,1.4d0*(i-1), 0,0,0), &
                pgto_new(0.6239137298d0, 0d0,0d0,1.4d0*(i-1), 0,0,0), &
                pgto_new(0.168855404d0, 0d0,0d0,1.4d0*(i-1), 0,0,0) &
            /), &
            (/ &
                0.1543289673d0, &
                0.5353281423d0, &
                0.4446345422d0 &
            /) )

        s%nucleuses(i)%cx = 0d0
        s%nucleuses(i)%cy = 0d0
        s%nucleuses(i)%cz = (i-1)*1.4d0
        s%nucleuses(i)%charge = 1d0
    end do

    call hf_run_situation(s, energy)
    write (*,*) energy


end program
