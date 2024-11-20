module mulliken
    use hf
    use hf_situation
    use matrix
    implicit none

    type mulliken_report
        real(8), allocatable, dimension(:) :: partial_charges
    end type

    contains

        subroutine mulliken_analysis(s, sres, ret)
            type(situation), intent(in) :: s
            type(situation_result), intent(in) :: sres
            type(mulliken_report), intent(out) :: ret

            real(8), allocatable, dimension(:,:) :: mat_s, mat_p, mat_d
            real(8), allocatable, dimension(:) :: vec_charge
            integer :: i, j, nidx

            allocate(mat_s(s%num_basis, s%num_basis))
            allocate(mat_p(s%num_basis, s%num_basis))
            allocate(mat_d(s%num_basis, s%num_basis))
            allocate(ret%partial_charges(s%num_nucleuses))
            allocate(vec_charge(s%num_nucleuses))

            call hf_overlap_matrix(s, mat_s)

            call _mulliken_create_density_matrix(mat_d, sres%orbital_coefficients, s%num_electrons, s%num_basis)

            call matrix_mult_normal_normal(mat_d, mat_s, s%num_basis, s%num_basis, s%num_basis, mat_p)

            vec_charge = 0d0
            do i=1,s%num_basis
                nidx = s%basis_nuc_index(i)
                vec_charge(nidx) = vec_charge(nidx) + mat_p(i,i)
            end do

            do i=1,s%num_nucleuses
                ret%partial_charges(i) = s%nucleuses(i)%charge - vec_charge(i)
            end do

            do i=1,s%num_basis
                write (*,*) i, mat_p(i,i)
            end do

            write (*,*) "DENSITY MATRIX"
                write (*,*) mat_d(1,1:2)
                write (*,*) mat_d(2,1:2)

            deallocate(mat_s)
            deallocate(mat_p)
            deallocate(mat_d)
            deallocate(vec_charge)
        end subroutine

        ! create density matrix for CLOSED system,
        ! which mean num_elec/2 orbitals with least energy are fully occupied
        subroutine _mulliken_create_density_matrix(mat_d, mat_c, num_elec, num_basis)
            real(8), dimension(num_basis,num_basis), intent(out) :: mat_d
            real(8), dimension(num_basis,num_basis), intent(in) :: mat_c
            integer, intent(in) :: num_elec, num_basis

            integer :: num_occupied
            integer :: mu, nu, i
            real(8) :: rtmp

            num_occupied = num_elec/2
            mat_d = 0d0
            do mu=1,num_basis
            do nu=1,num_basis
                rtmp = 0d0
                do i=1,num_occupied
                    rtmp = rtmp + mat_c(mu,i) * mat_c(nu,i)
                end do
                mat_d(mu,nu) = 2d0 * rtmp
            end do
            end do

        end subroutine

end module
