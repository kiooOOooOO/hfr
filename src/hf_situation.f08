module hf_situation
    use sto_ng
    use timer
    implicit none

    #include "mpif.h"

    type nucleus
        real(8) :: cx, cy, cz, charge
    end type

    type situation
        type(nucleus), allocatable, dimension(:) :: nucleuses
        integer :: num_nucleuses

        integer :: num_electrons

        type(sto_ng), allocatable, dimension(:) :: basis_functions
        integer :: num_basis

        real(8), allocatable, dimension(:) :: eri_table
    end type

    type situation_result
        real(8) :: molecular_energy

        real(8) :: nucleus_potential_energy

        real(8) :: electron_energy

        integer :: num_orbitals

        ! (i,j) -> coefficient to j-th basis function in i-th molecular orbital
        real(8), allocatable, dimension(:,:) :: orbital_coefficients

        real(8), allocatable, dimension(:) :: orbital_energies
    end type

    contains
        subroutine hf_dump_situation_result(s, sr, out_file_name)
            type(situation), intent(in) :: s
            type(situation_result), intent(in) :: sr
            character(:), allocatable :: out_file_name

            integer :: i, j

            open(21, file=out_file_name, status='replace')
            write (21,'(A)') "[Metadata]"
            write (21,'(A)') "Title=SituationResult"
            write (21,'("NumElectrons=",I0)') s%num_electrons
            write (21,'("NumBasis=",I0)') s%num_basis
            write (21,'("NumOrbitals=",I0)') sr%num_orbitals
            write (21,'("NumNucleus=",I0)') s%num_nucleuses
            write (21,'("MolecularEnergy=",E16.8)') sr%molecular_energy
            write (21,'("NucleusPotentialEnergy=",E16.8)') sr%nucleus_potential_energy
            write (21,'("ElectronEnergy=",E16.8)') sr%electron_energy

            write (21,'(A)') ""

            write (21,'(A)') "[Nucleus]"
            do i=1,s%num_nucleuses
                write (21,'(4(E15.4e2))') s%nucleuses(i)%cx, s%nucleuses(i)%cy, s%nucleuses(i)%cz, s%nucleuses(i)%charge
            end do
            write (21,'(A)') ""

            write (21,'(A)') "[Basis]"
            do i=1,s%num_basis
                write (21,'(I3,E15.4e2)') s%basis_functions(i)%n, s%basis_functions(i)%norm
                do j=1,s%basis_functions(i)%n
                    write (21,'(3(E15.4),3(I5),2(E15.4))') &
                        s%basis_functions(i)%pgtos(j)%cx, &
                        s%basis_functions(i)%pgtos(j)%cy, &
                        s%basis_functions(i)%pgtos(j)%cz, &
                        s%basis_functions(i)%pgtos(j)%nx, &
                        s%basis_functions(i)%pgtos(j)%ny, &
                        s%basis_functions(i)%pgtos(j)%nz, &
                        s%basis_functions(i)%pgtos(j)%norm * s%basis_functions(i)%coefs(j), &
                        s%basis_functions(i)%pgtos(j)%expo
                end do
            end do
            write (21,'(A)') ""

            write (21,'(A)') "[Orbitals]"
            do i=1,sr%num_orbitals
                write (21,'("E ",E15.4)') sr%orbital_energies(i)
                do j=1,s%num_basis
                    write (21,'(E15.4)') sr%orbital_coefficients(j, i)
                end do
            end do

            close(21)
        end subroutine

        subroutine hf_pick_calculation(s, calculation, calc_count)
            type(situation), intent(in) :: s
            integer, allocatable, dimension(:) :: calculation
            integer :: calc_count

            integer :: e1, e2, e3, e4, i, tmp, idx
            real :: rtmp
            logical, allocatable, dimension(:) :: calculated

            allocate(calculated(s%num_basis**4))

            calculated = .false.
            calculation = 0
            calc_count = 0

            do e1=1,s%num_basis
            do e2=1,s%num_basis
            do e3=1,s%num_basis
            do e4=1,s%num_basis
                if ( calculated(hf_eri_index(s,e1,e2,e3,e4)) ) then
                    cycle
                end if
                calculated(hf_eri_index(s,e1,e2,e3,e4)) = .true.
                calc_count = calc_count + 1
                calculation(calc_count) = hf_eri_index(s,e1,e2,e3,e4)
                
                calculated(hf_eri_index(s, e2, e1, e3, e4)) = .true.
                calculated(hf_eri_index(s, e2, e1, e4, e3)) = .true.
                calculated(hf_eri_index(s, e3, e4, e1, e2)) = .true.
                calculated(hf_eri_index(s, e3, e4, e2, e1)) = .true.
                calculated(hf_eri_index(s, e4, e3, e1, e2)) = .true.
                calculated(hf_eri_index(s, e4, e3, e2, e1)) = .true.
            end do
            end do
            end do
            end do

            do i=1,calc_count
                call random_number(rtmp)
                idx = 1 + (rtmp*calc_count)

                tmp = calculation(calc_count-i+1)
                calculation(calc_count-i+1) = calculation(idx)
                calculation(idx) = tmp
            end do

        end subroutine

        subroutine hf_do_electron_repulsion_integrals_mpi(s, wnum, rank)
            type(situation), intent(inout) :: s
            integer, allocatable, dimension(:) :: calculation
            integer :: calc_count, i

            integer :: wnum, rank, ierr
            integer, allocatable, dimension(:) :: calc_each
            integer :: calc_each_count
            real(8), allocatable, dimension(:) :: each_eri, eri_all
            real(8) :: duration
            integer :: e1, e2, e3, e4, each_calculated
            type(timer_instance) :: tim
            character(MPI_MAX_PROCESSOR_NAME) :: procName
            integer :: nameLen

            allocate(calculation(s%num_basis**4))
            call hf_pick_calculation(s, calculation, calc_count)
            if ( mod(calc_count, wnum) .eq. 0 ) then
                calc_each_count = calc_count/wnum
            else
                calc_each_count = 1 + (calc_count/wnum)
            end if
            allocate(eri_all(calc_count))
            allocate(calc_each(calc_each_count))
            allocate(each_eri(calc_each_count))

            if ( rank .eq. 0 ) then
                write (*,*) "====== ERI calculation ========================="
                write (*,*) "Total ERIs to calculate", calc_count
                write (*,*) "# of Peers", wnum
                write (*,*) "ERIs per Peer", calc_each_count
            end if

            calc_each = 0
            call mpi_scatter(calculation, calc_each_count, MPI_INTEGER, calc_each, &
                calc_each_count, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            call mpi_get_processor_name(procName, nameLen, ierr)

            call timer_start(tim)
            each_calculated = 0
            do i=1,calc_each_count
                if ( calc_each(i) .eq. 0 ) then
                    exit
                end if

                call hf_eri_index_to_number(s, calc_each(i), e1, e2, e3, e4)
                each_eri(i) = stong_eri( &
                    s%basis_functions(e1), &
                    s%basis_functions(e2), &
                    s%basis_functions(e3), &
                    s%basis_functions(e4))
                each_calculated = each_calculated + 1
            end do
            call timer_stop(tim, duration)
            write (*,'("rank",I5," on ", A, " calculated ", I5, " eris in ", F10.3, "sec")') rank, procName(1:nameLen), each_calculated, duration

            call mpi_gather(each_eri, calc_each_count, MPI_DOUBLE_PRECISION, &
                eri_all, calc_each_count, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

            if ( rank .eq. 0 ) then
                write (*,*) "================================================"
            end if

            do i=1,calc_count
                call hf_eri_index_to_number(s, calculation(i), e1, e2, e3, e4)
                s%eri_table(hf_eri_index(s, e1, e2, e3, e4)) = eri_all(i)
                s%eri_table(hf_eri_index(s, e2, e1, e3, e4)) = eri_all(i)
                s%eri_table(hf_eri_index(s, e2, e1, e4, e3)) = eri_all(i)
                s%eri_table(hf_eri_index(s, e3, e4, e1, e2)) = eri_all(i)
                s%eri_table(hf_eri_index(s, e3, e4, e2, e1)) = eri_all(i)
                s%eri_table(hf_eri_index(s, e4, e3, e1, e2)) = eri_all(i)
                s%eri_table(hf_eri_index(s, e4, e3, e2, e1)) = eri_all(i)
            end do

            deallocate(calculation)
            deallocate(calc_each)
            deallocate(each_eri)
        end subroutine

        pure function hf_eri_index(s, e1, e2, e3, e4)
            integer, intent(out) :: hf_eri_index
            type(situation), intent(in) :: s
            integer, intent(in) :: e1, e2, e3, e4

            integer :: ret, base, i

            base = 1

            ret = (e1-1)*base
            base = base*s%num_basis
            ret = (e2-1)*base + ret
            base = base*s%num_basis
            ret = (e3-1)*base + ret
            base = base*s%num_basis
            ret = (e4-1)*base + ret
            base = base*s%num_basis

            hf_eri_index = ret+1
        end function

        subroutine hf_eri_index_to_number(s, idx, e1, e2, e3, e4)
            type(situation), intent(in) :: s
            integer, intent(in) :: idx
            integer, intent(out) :: e1, e2, e3, e4

            integer :: tmp

            tmp = idx - 1
            e1 = mod(tmp, s%num_basis) + 1
            tmp = tmp / s%num_basis
            e2 = mod(tmp, s%num_basis) + 1
            tmp = tmp / s%num_basis
            e3 = mod(tmp, s%num_basis) + 1
            e4 = 1 + (tmp / s%num_basis)
        end subroutine

        pure function hf_eri_cache(s, e1, e2, e3, e4)
            real(8), intent(out) :: hf_eri_cache
            type(situation), intent(in) :: s
            integer, intent(in) :: e1, e2, e3, e4

            integer :: idx

            idx = hf_eri_index(s, e1, e2, e3, e4)

            hf_eri_cache = s%eri_table(idx)
        end function
end module
