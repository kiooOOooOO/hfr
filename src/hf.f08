module hf
    use matrix
    use sto_ng
    implicit none

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

    contains

        pure function _hf_eri_index(s, e1, e2, e3, e4)
            integer, intent(out) :: _hf_eri_index
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

            _hf_eri_index = ret+1
        end function

        subroutine hf_run_situation(s, energy)
            type(situation) :: s
            real(8) :: energy

            real(8), allocatable, dimension(:,:) :: mat_s, mat_sd, mat_u, mat_p, mat_ch, mat_f
            real(8), allocatable, dimension(:,:) :: mat_tmp, mat_x, mat_fd, mat_e, mat_cd, mat_c
            integer :: i, iterations
            real(8) :: energy_diff, last_energy

            write (*,*) "creating eri table...."
            call _hf_create_eri_table(s)
            write (*,*) "done"

            allocate(mat_s(s%num_basis, s%num_basis))
            allocate(mat_sd(s%num_basis, s%num_basis))
            allocate(mat_u(s%num_basis, s%num_basis))
            allocate(mat_p(s%num_basis, s%num_basis))
            allocate(mat_ch(s%num_basis, s%num_basis))
            allocate(mat_f(s%num_basis, s%num_basis))
            allocate(mat_tmp(s%num_basis, s%num_basis))
            allocate(mat_x(s%num_basis, s%num_basis))
            allocate(mat_fd(s%num_basis, s%num_basis))
            allocate(mat_e(s%num_basis, s%num_basis))
            allocate(mat_cd(s%num_basis, s%num_basis))
            allocate(mat_c(s%num_basis, s%num_basis))

            call hf_overlap_matrix(s, mat_s)
            call matrix_sym_diagonalize(mat_s, s%num_basis, mat_sd, mat_u)
            do i=1,s%num_basis
                mat_sd(i,i) = mat_sd(i,i)**-0.5d0
            end do

            call matrix_mult_normal_normal(mat_u, mat_sd, s%num_basis, s%num_basis, s%num_basis, mat_tmp)
            call matrix_mult_normal_transpose(mat_tmp, mat_u, s%num_basis, s%num_basis, s%num_basis, mat_x)

            call hf_core_hamiltonian_matrix(s, mat_ch)

            mat_c = 0d0
            mat_p = 0d0

            iterations = 0
            energy_diff = 1d0
            last_energy = 0d0
            do while ( energy_diff .gt. 1e-6 )
                iterations = iterations + 1
                write (*,*) "# of iterations", iterations
                call hf_fock_matrix(s, mat_ch, mat_c, mat_f)
                write (*,*) "fock matrix created"
                write (*,*) mat_f(1,1:2)
                write (*,*) mat_f(2,1:2)


                call matrix_mult_transpose_normal(mat_x, mat_f, s%num_basis, s%num_basis, s%num_basis, mat_tmp)
                call matrix_mult_normal_normal(mat_tmp, mat_x, s%num_basis, s%num_basis, s%num_basis, mat_fd)

                call hf_assert_matrix_symetric(s%num_basis, mat_fd)

                write (*,*) "transposed fock matrix"
                write (*,*) mat_fd(1,1:2)
                write (*,*) mat_fd(2,1:2)

                call matrix_sym_diagonalize(mat_fd, s%num_basis, mat_e, mat_cd)

                call hf_assert_hfr_answer(s%num_basis, mat_fd, mat_cd, mat_e)
                write (*,*) "transposed coefficient matrix"
                write (*,*) mat_cd(1,1:2)
                write (*,*) mat_cd(2,1:2)

                call matrix_mult_normal_normal(mat_x, mat_cd, s%num_basis, s%num_basis, s%num_basis, mat_c)
                write (*,*) "coefficient matrix"
                write (*,*) mat_c(1,1:2)
                write (*,*) mat_c(2,1:2)
                write (*,*) "energy matrix"
                write (*,*) mat_e(1,1:2)
                write (*,*) mat_e(2,1:2)
                call hf_density_matrix(s%num_basis, mat_c, mat_p)

                call hf_assert_hf_answer(s%num_basis, mat_f, mat_c, mat_s, mat_e)

                energy = hf_electron_energy(s, mat_ch, mat_p)
                if ( last_energy .ne. 0d0 ) then
                    energy_diff = abs(energy - last_energy)
                end if
                last_energy = energy
            end do

            write (*,*) "electron energy", energy
            energy = energy + hf_potential_energy(s)

            write (*,*) "########"
                write (*,*) "transposed fock matrix"
                write (*,*) mat_fd(1,1:2)
                write (*,*) mat_fd(2,1:2)

                write (*,*) "transposed coefficient matrix"
                write (*,*) mat_cd(1,1:2)
                write (*,*) mat_cd(2,1:2)

                write (*,*) "F'C'"
                call matrix_mult_normal_normal(mat_fd, mat_cd, 2, 2, 2, mat_tmp)
                write (*,*) mat_tmp(1,1:2)
                write (*,*) mat_tmp(2,1:2)
            write (*,*) "########"

            deallocate(mat_s)
            deallocate(mat_sd)
            deallocate(mat_u)
            deallocate(mat_p)
            deallocate(mat_ch)
            deallocate(mat_f)
            deallocate(mat_tmp)
            deallocate(mat_x)
            deallocate(mat_fd)
            deallocate(mat_e)
            deallocate(mat_cd)
            deallocate(mat_c)
        end subroutine

        subroutine hf_assert_matrix_symetric(n, mat)!{{{
            integer, intent(in) :: n
            real(8), dimension(:,:) :: mat

            integer :: r, c
            logical :: assert

            assert = .false.

            do r=1,n
            do c=1,n
                if ( abs(mat(r,c) - mat(c,r)) .gt. 1e-4 ) then
                    write (*,*) "symmetric matrix not symmetric"
                    write (*,*) r, c, mat(r,c)
                    write (*,*) c, r, mat(c,r)
                    assert = .true.
                end if
            end do
            end do

            if ( assert ) then
                stop
            end if
        end subroutine!}}}

        subroutine hf_assert_hf_answer(n, mat_f, mat_c, mat_s, mat_e)!{{{
            integer, intent(in) :: n
            real(8), dimension(:,:) :: mat_f, mat_c, mat_s, mat_e

            real(8), allocatable, dimension(:,:) :: mat_left, mat_right, mat_tmp
            integer :: r, c

            allocate(mat_left(n,n))
            allocate(mat_right(n,n))
            allocate(mat_tmp(n,n))

            call matrix_mult_normal_normal(mat_f, mat_c, n, n, n, mat_left)
            call matrix_mult_normal_normal(mat_s, mat_c, n, n, n, mat_tmp)
            call matrix_mult_normal_normal(mat_tmp, mat_e, n, n, n, mat_right)
            do r=1,n
            do c=1,n
                if ( abs(mat_left(r,c) - mat_right(r,c)) .gt. 1e-4 ) then
                    write (*,*) "hf error"
                    write (*,*) "left"
                    write (*,*) mat_left

                    write (*,*) "right"
                    write (*,*) mat_right

                    stop
                end if
            end do
            end do

            write (*,*) "Hartree-Fock assertion passed"

            deallocate(mat_left)
            deallocate(mat_right)
            deallocate(mat_tmp)
        end subroutine!}}}

        subroutine hf_assert_hfr_answer(n, mat_fd, mat_cd, mat_e)!{{{
            integer, intent(in) :: n
            real(8), dimension(:,:) :: mat_fd, mat_cd, mat_e

            real(8), allocatable, dimension(:,:) :: mat_left, mat_right
            integer :: r, c, i

            allocate(mat_left(n,n))
            allocate(mat_right(n,n))

            call matrix_mult_normal_normal(mat_fd, mat_cd, n, n, n, mat_left)
            call matrix_mult_normal_normal(mat_cd, mat_e, n, n, n, mat_right)

            do r=1,n
            do c=1,n
                if ( abs(mat_left(r,c) - mat_right(r,c)) .gt. 1e-4 ) then
                    write (*,*) "hfr error"
                    write (*,*) "left"
                    write (*,*) mat_left

                    write (*,*) "right"
                    write (*,*) mat_right

                    stop
                end if
            end do
            end do

            deallocate(mat_left)
            deallocate(mat_right)
        end subroutine!}}}

        subroutine hf_density_matrix(n, mat_c, mat)
            integer, intent(in) :: n
            real(8), intent(in), dimension(:,:) :: mat_c
            real(8), intent(out), dimension(:,:) :: mat

            integer :: r, c, i, j
            real(8) :: val

            mat = 0d0
            do concurrent(r=1:n, c=1:n) local(val, j)
                val = 0d0
                do j=1,n/2
                    val = val + 2*mat_c(c,j)*mat_c(r,j)
                end do

                mat(r,c) = val
            end do
        end subroutine

        function hf_electron_energy(s, mat_ch, mat_p) result(energy)
            real(8), intent(out) :: energy
            type(situation), intent(inout) :: s
            real(8), intent(in), dimension(:,:) :: mat_ch, mat_p

            integer :: nb
            integer :: k, l, m, n

            real(8) :: val

            nb = s%num_basis
            val = 0d0

            do k=1,nb
            do l=1,nb
                val = val + mat_p(l,k)*mat_ch(k,l)
            end do
            end do

            do k=1,nb
            do l=1,nb
            do m=1,nb
            do n=1,nb
                val = val + 0.5d0*mat_p(l,k)*mat_p(n,m)*(_hf_eri_cache(s, k, m, l, n) - 0.5d0*_hf_eri_cache(s, k, m, n, l))
            end do
            end do
            end do
            end do

            energy = val
        end function

        function hf_potential_energy(s) result(energy)
            real(8), intent(out) :: energy
            type(situation), intent(in) :: s

            integer :: i1, i2
            real(8) :: val, dist

            val = 0d0

            do i1=1,s%num_basis
            do i2=i1+1,s%num_basis
                dist = sqrt((s%nucleuses(i1)%cx-s%nucleuses(i2)%cx)**2 + (s%nucleuses(i1)%cy-s%nucleuses(i2)%cy)**2 &
                    + (s%nucleuses(i1)%cz-s%nucleuses(i2)%cz)**2)
                val = val + s%nucleuses(i1)%charge * s%nucleuses(i2)%charge / dist
            end do
            end do

            energy = val
        end function

        subroutine hf_overlap_matrix(s, mat)
            type(situation), intent(in) :: s
            real(8), intent(out), dimension(:,:) :: mat

            integer :: r, c, i
            real(8) :: val

            do i=0,s%num_basis**2 - 1
                r = (i/s%num_basis) + 1
                c = mod(i,s%num_basis) + 1

                val = stong_overlap(s%basis_functions(r), s%basis_functions(c))
                mat(r,c) = val
            end do
        end subroutine

        subroutine hf_core_hamiltonian_matrix(s, mat)
            type(situation), intent(in) :: s
            real(8), intent(out), dimension(:,:) :: mat

            integer :: r, c, i, ni
            real(8) :: val

            do i=0,s%num_basis**2 - 1
                r = (i/s%num_basis) + 1
                c = mod(i,s%num_basis) + 1

                val = stong_kinetic_energy(s%basis_functions(r), s%basis_functions(c))
                do ni=1,s%num_nucleuses
                    val = val - s%nucleuses(ni)%charge*stong_nuclear_attr( &
                        s%basis_functions(r), s%basis_functions(c), &
                        s%nucleuses(ni)%cx, &
                        s%nucleuses(ni)%cy, &
                        s%nucleuses(ni)%cz &
                        )
                end do

                mat(r,c) = val
            end do
        end subroutine

        subroutine hf_fock_matrix(s, ch, matC, mat)
            type(situation), intent(inout) :: s
            real(8), intent(in), dimension(:,:) :: ch, matC
            real(8), intent(out), dimension(:,:) :: mat

            integer :: k, l, j, m, n
            real(8) :: val

            do k=1,s%num_basis
            do l=1,s%num_basis
                val = ch(k,l)

                do j=1,s%num_electrons/2
                    do m=1,s%num_basis
                    do n=1,s%num_basis
                        val = val + (2d0*_hf_eri_cache(s,k,m,l,n)-_hf_eri_cache(s,k,m,n,l))*matC(m,j)*matC(n,j)
                    end do
                    end do
                end do

                mat(k,l) = val
            end do
            end do

!            integer :: r, c, i, j
!            real(8) :: val
!
!            do r=1,s%num_basis
!            do c=1,s%num_basis
!                val = ch(r,c)
!                do i=1,s%num_basis
!                do j=1,s%num_basis
!                    val = val + matP(i,j)*(_hf_eri_cache(s, r, c, i, j) - 0.5d0*_hf_eri_cache(s,r,j,i,c))
!                end do
!                end do
!
!!                val = ch(r,c) &
!!                    + _hf_eri_cache(s, r, c, r, r)*matP(r,r) &
!!                    + _hf_eri_cache(s, r, c, r, c)*matP(r,c) &
!!                    + _hf_eri_cache(s, r, c, c, r)*matP(c,r) &
!!                    + _hf_eri_cache(s, r, c, c, c)*matP(c,c) &
!!                    - 0.5d0*_hf_eri_cache(s, r, r, r, c)*matP(r,r) &
!!                    - 0.5d0*_hf_eri_cache(s, r, r, c, c)*matP(c,r) &
!!                    - 0.5d0*_hf_eri_cache(s, r, c, r, c)*matP(r,c) &
!!                    - 0.5d0*_hf_eri_cache(s, r, c, c, c)*matP(c,c)
!
!                mat(r,c) = val
!            end do
!            end do
        end subroutine

#define BF(x) s%basis_functions(x)
        subroutine hf_test_eri(s)
            type(situation), intent(inout) :: s

            real(8) :: val

            val = stong_eri(BF(3), BF(2), BF(6), BF(1))
            write (*,*) 3261, val
        end subroutine

        subroutine _hf_create_eri_table(s)
            type(situation), intent(inout) :: s
            integer :: e1, e2, e3, e4, idx, n
            real(8) :: val

            n = s%num_basis
            do concurrent(e1=1:n, e2=1:n, e3=1:n, e4=1:n) local(idx, val)
                idx = _hf_eri_index(s, e1, e2, e3, e4)
                val = stong_eri(BF(e1), BF(e2), BF(e3), BF(e4))

                s%eri_table(idx) = val
            end do

            do e1=1,n
            do e2=1,n
            do e3=1,n
            do e4=1,n
                idx = _hf_eri_index(s, e1, e2, e3, e4)
                write (*,'(I6,I6,I6,I6,I6,E16.8)') idx, e1, e2, e3, e4, s%eri_table(idx)
            end do
            end do
            end do
            end do
        end subroutine

        pure function _hf_eri_cache(s, e1, e2, e3, e4)
            real(8), intent(out) :: _hf_eri_cache
            type(situation), intent(in) :: s
            integer, intent(in) :: e1, e2, e3, e4

            integer :: idx

            idx = _hf_eri_index(s, e1, e2, e3, e4)

            _hf_eri_cache = s%eri_table(idx)
        end function
#undef BF
end module
