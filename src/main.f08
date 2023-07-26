program hfr_main
    use hf
    use pgto
    implicit none

    type(situation) :: s
    type(pgto) :: g1, g2, g3
    real(8) :: energy, val

    integer :: i

    write (*,*) "Hello World!"

    s%num_nucleuses = 2
    allocate(s%nucleuses(2))

    ! Lithium
    s%nucleuses(1)%cx = 0d0
    s%nucleuses(1)%cy = 0d0
    s%nucleuses(1)%cz = -1.3228d0
    s%nucleuses(1)%charge = 1d0

    ! Hydrogen
    s%nucleuses(2)%cx = 0d0
    s%nucleuses(2)%cy = 0d0
    s%nucleuses(2)%cz = 1.3228d0
    s%nucleuses(2)%charge = 1d0

    s%num_basis = 2
    allocate(s%basis_functions(2))

    ! H-1s
    s%basis_functions(1) = stong_new(3, &
        (/ &
            pgto_new(3.425250914d0,  0d0,0d0,-1.3228d0, 0,0,0), &
            pgto_new(0.6239137298d0, 0d0,0d0,-1.3228d0, 0,0,0), &
            pgto_new(0.168855404d0,  0d0,0d0,-1.3228d0, 0,0,0) &
        /), &
        (/ &
            0.1543289673d0, &
            0.5353281423d0, &
            0.4446345422d0 &
        /) )

!    write (*,*) "sto-ng norm", s%basis_functions(1)%norm
!    do i=1,3
!    write (*,*) "pgto norm", i, s%basis_functions(1)%pgtos(i)%norm
!    end do
!    stop

    s%basis_functions(2) = stong_new(3, &
        (/ &
            pgto_new(3.425250914d0,  0d0,0d0,1.3228d0, 0,0,0), &
            pgto_new(0.6239137298d0, 0d0,0d0,1.3228d0, 0,0,0), &
            pgto_new(0.168855404d0,  0d0,0d0,1.3228d0, 0,0,0) &
        /), &
        (/ &
            0.1543289673d0, &
            0.5353281423d0, &
            0.4446345422d0 &
        /) )

!    ! Li-1s
!    s%basis_functions(2) = stong_new(3, &
!        (/ &
!            pgto_new(0.1611957475d2, 0d0,0d0,0d0, 0,0,0), &
!            pgto_new(0.2936200663d1, 0d0,0d0,0d0, 0,0,0), &
!            pgto_new(0.7946504870d0, 0d0,0d0,0d0, 0,0,0) &
!        /), &
!        (/ &
!            0.1543289673d0, &
!            0.5353281423d0, &
!            0.4446345422d0 &
!        /) )
!
!    ! Li-2s
!    s%basis_functions(1) = stong_new(3, &
!        (/ &
!            pgto_new(0.1611957475d2, 0d0,0d0,0d0, 0,0,0), &
!            pgto_new(0.2936200663d1, 0d0,0d0,0d0, 0,0,0), &
!            pgto_new(0.7946504870d0, 0d0,0d0,0d0, 0,0,0) &
!        /), &
!        (/ &
!            -0.9996722919d-1, &
!            0.3995128261d0, &
!            0.7001154689d0 &
!        /) )

    allocate(s%eri_table(s%num_basis**4))
    s%eri_table = 0d0

    call hf_run_situation(s, energy)
    write (*,*) "molecular energy", energy


end program
