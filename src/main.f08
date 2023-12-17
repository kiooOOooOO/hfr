program hfr_main
    use hf
    use pgto
    implicit none

    integer :: ierr
    integer :: wnum, rank
    type(situation) :: dummy

    call mpi_init(ierr)
    if ( ierr .ne. 0 ) then
        write (*,*) "mpi_init error", ierr
        stop
    end if

    call mpi_comm_size(MPI_COMM_WORLD, wnum, ierr)
    if ( ierr .ne. 0 ) then
        write (*,*) "mpi_comm_size error", ierr
        call mpi_finalize(ierr)
        stop
    end if

    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    if ( ierr .ne. 0 ) then
        write (*,*) "mpi_comm_rank error", ierr
        call mpi_finalize(ierr)
        stop
    end if

    call exec(wnum, rank)

    call mpi_finalize(ierr)

    contains

        subroutine exec(wnum, rank)
            integer :: wnum, rank
            intrinsic :: command_argument_count, get_command_argument
            integer :: stat, length

            !    character, allocatable, dimension(:) :: in_file_name
            character(:), allocatable :: in_file_name, out_file_name, dumpstr
            integer :: argc
            logical :: dump

            if ( rank .eq. 0 ) then
                write (*,*) "Hello World!"
            end if

            argc = command_argument_count()
            if ( argc .lt. 1 ) then
                write (*,*) "need input file name and output file name"
                return
            end if
            if ( argc .lt. 2 ) then
                write (*,*) "need output file name"
                return
            end if

            call get_command_argument(1, status=stat, length=length)
            if ( stat .ne. 0 ) then
                return
            end if
            allocate(character(length) :: in_file_name)
            call get_command_argument(1, value=in_file_name, status=stat)
            if ( stat .ne. 0 ) then
                return
            end if

            call get_command_argument(2, status=stat, length=length)
            if ( stat .ne. 0 ) then
                return
            end if
            allocate(character(length) :: out_file_name)
            call get_command_argument(2, value=out_file_name, status=stat)
            if ( stat .ne. 0 ) then
                return
            end if

            call get_command_argument(3, status=stat, length=length)
            if ( stat .eq. 0 ) then
                allocate(character(length) :: dumpstr)
                call get_command_argument(3, value=dumpstr, status=stat)
                if ( stat .eq. 0 ) then
                    dump = (dumpstr .eq. 'yes')
                else
                    dump = .false.
                end if
            else
                dump = .false.
            end if


            if ( rank .eq. 0 ) then
                write (*,*) "IN_FILE: ", in_file_name
                write (*,*) "OUT_FILE: ", out_file_name
                write (*,*) "DUMP_ITERATIONS: ", dump
            end if
            call main(in_file_name, out_file_name, dump, wnum, rank)

            deallocate(in_file_name)
        end subroutine

        subroutine main(in_file_name, out_file_name, dump, wnum, rank)
            character(:), allocatable :: in_file_name, out_file_name
            logical :: dump
            integer :: wnum, rank, i

            type(situation) :: s
            type(situation_result) :: sr

            call load_situation(s, in_file_name)

            allocate(s%eri_table(s%num_basis**4))

            do i=1,LOOP_COUNT
            s%eri_table = 0d0

            call hf_do_electron_repulsion_integrals_mpi(s, wnum, rank)
            end do

            if ( rank .eq. 0 ) then
                call hf_run_situation(s, sr, dump)

                call hf_dump_situation_result(s, sr, out_file_name)
            end if
        end subroutine

        subroutine load_situation(situ, in_file_name)
            type(situation), intent(inout) :: situ
            character(:), intent(in), allocatable :: in_file_name
            integer :: i, j, nuc_idx
            integer :: n_nuc, n_basis, n_pgto, n_elec
            integer :: nx, ny, nz
            real(8) :: expo, coef, charge, x, y, z
            real(8), allocatable, dimension(:) :: coefs
            type(pgto), allocatable, dimension(:) :: pgtos

            open(17, file=in_file_name, status='old')

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
                read (17,*) nuc_idx
                read (17,*) n_pgto
                allocate(coefs(n_pgto))
                allocate(pgtos(n_pgto))

                do j=1,n_pgto
                    read (17,*) expo, nx, ny, nz, coef

                    pgtos(j) = pgto_new(expo, &
                        situ%nucleuses(nuc_idx)%cx, situ%nucleuses(nuc_idx)%cy, situ%nucleuses(nuc_idx)%cz, &
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
