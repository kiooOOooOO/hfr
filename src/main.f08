program hfr_main
    use hf
    use pgto
    implicit none

    write (*,*) "Hello World!"
    call main()

    contains

        subroutine main()
            type(situation) :: s
            type(situation_result) :: sr

            call load_situation(s)

            allocate(s%eri_table(s%num_basis**4))
            s%eri_table = 0d0
            call hf_run_situation(s, sr)

            call hf_dump_situation_result(s, sr)
        end subroutine

        subroutine load_situation(situ)
            type(situation), intent(inout) :: situ
            integer :: i, j
            integer :: n_nuc, n_basis, n_pgto, n_elec
            integer :: nx, ny, nz
            real(8) :: expo, coef, charge, x, y, z
            real(8), allocatable, dimension(:) :: coefs
            type(pgto), allocatable, dimension(:) :: pgtos

            open(17, file='./inputs/h2.dat', status='old')

            read (17,*) n_nuc
            situ%num_nucleuses = n_nuc
            allocate(situ%nucleuses(n_nuc))
            do i=1,n_nuc
                read (17,*) x, y, z, charge
                situ%nucleuses(i)%cx = x
                situ%nucleuses(i)%cy = y
                situ%nucleuses(i)%cz = z
                situ%nucleuses(i)%charge = charge
            end do

            read (17,*) n_basis
            situ%num_basis = n_basis
            allocate(situ%basis_functions(n_basis))
            do i=1,n_basis
                read (17,*) n_pgto
                allocate(coefs(n_pgto))
                allocate(pgtos(n_pgto))

                do j=1,n_pgto
                    read (17,*) expo, nx, ny, nz, coef

                    pgtos(j) = pgto_new(expo, &
                        situ%nucleuses(i)%cx, situ%nucleuses(i)%cy, situ%nucleuses(i)%cz, &
                        nx, ny, nz)
                    coefs(j) = coef
                end do

                situ%basis_functions(i) = stong_new(n_pgto, pgtos, coefs)

                deallocate(coefs)
                deallocate(pgtos)
            end do

            read (17,*) n_elec
            situ%num_electrons = n_elec

            close(17)
        end subroutine

end program
