program hfr_main
    use hf
    use pgto
    implicit none

    type(situation) :: s
    real(8) :: energy

    integer :: i

    write (*,*) "Hello World!"

    s%num_nucleuses = 3
    allocate(s%nucleuses(3))
    s%nucleuses(1)%cx = 0d0
    s%nucleuses(1)%cy = 0d0
    s%nucleuses(1)%cz = 0d0
    s%nucleuses(1)%charge = 8d0

    s%nucleuses(2)%cx = -1.8066d0
    s%nucleuses(2)%cy = 0d0
    s%nucleuses(2)%cz = 0d0
    s%nucleuses(2)%charge = 1d0

    s%nucleuses(3)%cx = 0.4523d0
    s%nucleuses(3)%cy = 1.7491d0
    s%nucleuses(3)%cz = 0d0
    s%nucleuses(3)%charge = 1d0

    s%num_basis = 6
    allocate(s%basis_functions(6))

    ! H-1s (left)
    s%basis_functions(1) = stong_new(3, &
        (/ &
            pgto_new(3.425250914d0,  -1.8066d0,0d0,0d0, 0,0,0), &
            pgto_new(0.6239137298d0, -1.8066d0,0d0,0d0, 0,0,0), &
            pgto_new(0.168855404d0,  -1.8066d0,0d0,0d0, 0,0,0) &
        /), &
        (/ &
            0.1543289673d0, &
            0.5353281423d0, &
            0.4446345422d0 &
        /) )

    ! H-1s (right)
    s%basis_functions(2) = stong_new(3, &
        (/ &
            pgto_new(3.425250914d0,  0.4523d0,1.7491d0,0d0, 0,0,0), &
            pgto_new(0.6239137298d0, 0.4523d0,1.7491d0,0d0, 0,0,0), &
            pgto_new(0.168855404d0,  0.4523d0,1.7491d0,0d0, 0,0,0) &
        /), &
        (/ &
            0.1543289673d0, &
            0.5353281423d0, &
            0.4446345422d0 &
        /) )
    ! O-2s
    s%basis_functions(3) = stong_new(3, &
        (/ &
            pgto_new(5.03315d0, 0d0,0d0,0d0, 0,0,0), &
            pgto_new(1.11695d0, 0d0,0d0,0d0, 0,0,0), &
            pgto_new(0.38038d0, 0d0,0d0,0d0, 0,0,0) &
        /), &
        (/ &
            -0.09996d0, &
            0.399512d0, &
            0.70011d0 &
        /) )
    ! O-2px
    s%basis_functions(4) = stong_new(3, &
        (/ &
            pgto_new(5.03315d0, 0d0,0d0,0d0, 1,0,0), &
            pgto_new(1.11695d0, 0d0,0d0,0d0, 1,0,0), &
            pgto_new(0.38038d0, 0d0,0d0,0d0, 1,0,0) &
        /), &
        (/ &
            0.15591d0, &
            0.60768d0, &
            0.39195d0 &
        /) )
    ! O-2py
    s%basis_functions(5) = stong_new(3, &
        (/ &
            pgto_new(5.03315d0, 0d0,0d0,0d0, 0,1,0), &
            pgto_new(1.11695d0, 0d0,0d0,0d0, 0,1,0), &
            pgto_new(0.38038d0, 0d0,0d0,0d0, 0,1,0) &
        /), &
        (/ &
            0.15591d0, &
            0.60768d0, &
            0.39195d0 &
        /) )
    ! O-2pz
    s%basis_functions(6) = stong_new(3, &
        (/ &
            pgto_new(5.03315d0, 0d0,0d0,0d0, 0,0,1), &
            pgto_new(1.11695d0, 0d0,0d0,0d0, 0,0,1), &
            pgto_new(0.38038d0, 0d0,0d0,0d0, 0,0,1) &
        /), &
        (/ &
            0.15591d0, &
            0.60768d0, &
            0.39195d0 &
        /) )

    allocate(s%eri_table(6**4))
    s%eri_table = 0d0

!    energy = stong_eri(s%basis_functions(6), s%basis_functions(5), s%basis_functions(5), s%basis_functions(5))

    call hf_run_situation(s, energy)
    write (*,*) energy


end program
