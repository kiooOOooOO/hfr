#define SZABO
#define CHECKED
module hf
    use hf_situation
    use matrix
    use sto_ng
    implicit none

    integer, parameter :: HF_MAX_ITERATION = 1000

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

        subroutine hf_run_situation(s, res, dump)
            type(situation) :: s
            type(situation_result) :: res
            logical :: dump
            real(8) :: energy

            real(8), allocatable, dimension(:,:) :: mat_s, mat_sd, mat_u, mat_p, mat_ch, mat_f, mat_pnext
            real(8), allocatable, dimension(:,:) :: mat_tmp, mat_x, mat_fd, mat_e, mat_cd, mat_c
            integer :: i, iterations

!            write (*,*) "creating eri table...."
!            call _hf_create_eri_table(s)
!            write (*,*) "done"

            allocate(mat_s(s%num_basis, s%num_basis))
            allocate(mat_sd(s%num_basis, s%num_basis))
            allocate(mat_u(s%num_basis, s%num_basis))
            allocate(mat_ch(s%num_basis, s%num_basis))
            allocate(mat_f(s%num_basis, s%num_basis))
            allocate(mat_x(s%num_basis, s%num_basis))
            allocate(mat_fd(s%num_basis, s%num_basis))
            allocate(mat_e(s%num_basis, s%num_basis))
            allocate(mat_p(s%num_basis, s%num_basis))
            allocate(mat_pnext(s%num_basis, s%num_basis))
            allocate(mat_cd(s%num_basis, s%num_basis))
            allocate(mat_c(s%num_basis, s%num_basis))

            allocate(mat_tmp(s%num_basis, s%num_basis))

            call hf_overlap_matrix(s, mat_s)
            call matrix_sym_diagonalize(mat_s, s%num_basis, mat_sd, mat_u)
            do i=1,s%num_basis
                mat_sd(i,i) = mat_sd(i,i)**-0.5d0
            end do

            call matrix_mult_normal_normal(mat_u, mat_sd, s%num_basis, s%num_basis, s%num_basis, mat_tmp)
            call matrix_mult_normal_transpose(mat_tmp, mat_u, s%num_basis, s%num_basis, s%num_basis, mat_x)

            call hf_core_hamiltonian_matrix(s, mat_ch)
            write (*,*) "CORE HAMILTONIAL MATRIX"
                write (*,*) mat_ch(1,1:2)
                write (*,*) mat_ch(2,1:2)

            mat_c = 0d0
            mat_p = 0d0

            iterations = 0
            do while ( .true. )
                iterations = iterations + 1

                if ( iterations .gt. HF_MAX_ITERATION ) then
                    write (*,*) "SCF failed"
                    stop
                end if

                call hf_fock_matrix(s, mat_ch, mat_c, mat_f)

                ! F' <- Xt*F*X
                call matrix_mult_transpose_normal(mat_x, mat_f, s%num_basis, s%num_basis, s%num_basis, mat_tmp)
                call matrix_mult_normal_normal(mat_tmp, mat_x, s%num_basis, s%num_basis, s%num_basis, mat_fd)

                call hf_assert_matrix_symetric(s%num_basis, mat_fd)

                ! F'C' = C'e
                call matrix_sym_diagonalize(mat_fd, s%num_basis, mat_e, mat_cd)
                call hf_assert_hfr_answer(s%num_basis, mat_fd, mat_cd, mat_e)

                ! C <- XC'
                call matrix_mult_normal_normal(mat_x, mat_cd, s%num_basis, s%num_basis, s%num_basis, mat_c)

                call hf_assert_hf_answer(s%num_basis, mat_f, mat_c, mat_s, mat_e)

                call hf_density_matrix(s, mat_c, mat_pnext)

                if ( dump ) then
                    call hf_dump_iteration(s, iterations, mat_f, mat_c, mat_e)
                end if

                write (*,'(I5,"-th iteration, energy=", E20.8)') iterations, energy

                if ( hf_converged(s, mat_p, mat_pnext) ) then
                    exit
                end if

                mat_p = mat_pnext
                energy = hf_electron_energy(s, mat_ch, mat_f, mat_p)
            end do

            res%nucleus_potential_energy = hf_potential_energy(s)
            res%electron_energy = energy
            res%molecular_energy = res%electron_energy + res%nucleus_potential_energy

            res%num_orbitals = s%num_basis

            allocate(res%orbital_coefficients(s%num_basis, res%num_orbitals))
            res%orbital_coefficients = mat_c

            allocate(res%orbital_energies(res%num_orbitals))
            do i=1,res%num_orbitals
                res%orbital_energies(i) = mat_e(i,i)
            end do

            deallocate(mat_s)
            deallocate(mat_sd)
            deallocate(mat_u)
            deallocate(mat_p)
            deallocate(mat_pnext)
            deallocate(mat_ch)
            deallocate(mat_f)
            deallocate(mat_tmp)
            deallocate(mat_x)
            deallocate(mat_fd)
            deallocate(mat_e)
            deallocate(mat_cd)
            deallocate(mat_c)
        end subroutine

        pure function hf_converged(s, mat_p1, mat_p2) result(converged)!{{{
            logical, intent(out) :: converged
            type(situation), intent(in) :: s
            real(8), intent(in), dimension(:,:) :: mat_p1, mat_p2

            real(8) :: sd
            integer :: r, c

            sd = 0d0
            do r=1,s%num_basis
            do c=1,s%num_basis
                sd = sd + (mat_p1(r,c) - mat_p2(r,c))**2d0
            end do
            end do

            sd = sqrt(sd * s%num_basis**-2d0)

            converged = sd .lt. 1d-4
        end function!}}}

        subroutine hf_dump_iteration(s, idx, mat_f, mat_c, mat_e)!{{{
            type(situation), intent(in) :: s
            integer, intent(in) :: idx
            real(8), dimension(:,:) :: mat_f, mat_c, mat_e

            integer :: r, c
            character :: filename*128

            write (filename, '("iterations/mfock_", I5.5, ".dat")') idx
            open(18, file=filename, status='replace')

            do r=1,s%num_basis
            do c=1,s%num_basis
                write (18,*) r, c, mat_f(r,c)
            end do
            write (18,*) ""
            end do

            close(18)

            write (filename, '("iterations/mcoef_", I5.5, ".dat")') idx
            open(18, file=filename, status='replace')

            do r=1,s%num_basis
            do c=1,s%num_basis
                write (18,*) r, c, mat_c(r,c)
            end do
            write (18,*) ""
            end do

            close(18)

            write (filename, '("iterations/energy_", I5.5, ".dat")') idx
            open(18, file=filename, status='replace')

            do r=1,s%num_basis
                write (18,*) r, mat_e(r,r)
            end do

            close(18)
        end subroutine!}}}

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
                if ( abs(mat_left(r,c) - mat_right(r,c)) .gt. 1e-2 ) then
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
                if ( abs(mat_left(r,c) - mat_right(r,c)) .gt. 1e-2 ) then
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

        SZABO subroutine hf_density_matrix(s, mat_c, mat)!{{{
            type(situation), intent(in) :: s
            real(8), intent(in), dimension(:,:) :: mat_c
            real(8), intent(out), dimension(:,:) :: mat

            integer :: r, c, i, j
            real(8) :: val

            mat = 0d0
            do concurrent(r=1:s%num_basis, c=1:s%num_basis) local(val, j)
                val = 0d0
                do j=1,s%num_electrons/2
                    val = val + mat_c(c,j)*mat_c(r,j)
                end do

                mat(r,c) = 2d0 * val
            end do
        end subroutine!}}}

        SZABO pure function hf_electron_energy(s, mat_ch, mat_f, mat_p) result(energy)!{{{
            real(8), intent(out) :: energy
            type(situation), intent(in) :: s
            real(8), intent(in), dimension(:,:) :: mat_ch, mat_f, mat_p

            integer :: r, c
            real(8) :: val

            val = 0d0
            do r=1,s%num_basis
            do c=1,s%num_basis
                val = val + mat_p(r,c)*(mat_ch(r,c) + mat_f(r,c))
            end do
            end do

            energy = 0.5d0 * val
        end function!}}}

        CHECKED function hf_potential_energy(s) result(energy)!{{{
            real(8), intent(out) :: energy
            type(situation), intent(in) :: s

            integer :: i1, i2
            real(8) :: val, dist

            val = 0d0

            do i1=1,s%num_nucleuses
                do i2=i1+1,s%num_nucleuses
                dist = sqrt((s%nucleuses(i1)%cx-s%nucleuses(i2)%cx)**2 + (s%nucleuses(i1)%cy-s%nucleuses(i2)%cy)**2 &
                    + (s%nucleuses(i1)%cz-s%nucleuses(i2)%cz)**2)
                val = val + s%nucleuses(i1)%charge * s%nucleuses(i2)%charge / dist
                end do
            end do

            energy = val
        end function!}}}

        subroutine hf_overlap_matrix(s, mat)!{{{
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
        end subroutine!}}}

        subroutine hf_core_hamiltonian_matrix(s, mat)!{{{
            type(situation), intent(in) :: s
            real(8), intent(out), dimension(:,:) :: mat

            integer :: r, c, i, ni
            real(8) :: val

            do r=1,s%num_basis
            do c=1,s%num_basis
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
            end do
        end subroutine!}}}

        SZABO subroutine hf_fock_matrix(s, ch, matC, mat)!{{{
            type(situation), intent(inout) :: s
            real(8), intent(in), dimension(:,:) :: ch, matC
            real(8), intent(out), dimension(:,:) :: mat

            integer :: k, l, j, m, n, a
            real(8) :: val

            ! Szabo (3.154)
            do m=1,s%num_basis
            do n=1,s%num_basis
                val = ch(m,n)

                do a=1,s%num_electrons/2
                    do k=1,s%num_basis
                    do l=1,s%num_basis
                        val = val + matC(k,a)*matC(l,a)*(2d0*_hf_eri_cache(s,m,n,k,l) - _hf_eri_cache(s,m,l,n,k))
                    end do
                    end do
                end do

                mat(m,n) = val
            end do
            end do
        end subroutine!}}}

#define BF(x) s%basis_functions(x)
        subroutine hf_test_eri(s)
            type(situation), intent(inout) :: s

            real(8) :: val

            val = stong_eri(BF(3), BF(2), BF(6), BF(1))
            write (*,*) 3261, val
        end subroutine

        subroutine _hf_create_eri_table(s)!{{{
            type(situation), intent(inout) :: s
            integer :: e1, e2, e3, e4, idx, n, i
            real(8) :: val
            logical, allocatable, dimension(:) :: calculated

            n = s%num_basis
            allocate(calculated(n**6))
            calculated = .false.

            do e1=1,n
            do e2=1,n
            do e3=1,n
            do e4=1,n
                idx = _hf_eri_index(s, e1, e2, e3, e4)
                if ( calculated(idx) ) then
                    cycle
                end if

                val = stong_eri(BF(e1), BF(e2), BF(e3), BF(e4))

                s%eri_table(idx) = val

                idx = _hf_eri_index(s, e2, e1, e3, e4)
                s%eri_table(idx) = val
                calculated(idx) = .true.
                idx = _hf_eri_index(s, e2, e1, e4, e3)
                s%eri_table(idx) = val
                calculated(idx) = .true.
                idx = _hf_eri_index(s, e3, e4, e1, e2)
                s%eri_table(idx) = val
                calculated(idx) = .true.
                idx = _hf_eri_index(s, e3, e4, e2, e1)
                s%eri_table(idx) = val
                calculated(idx) = .true.
                idx = _hf_eri_index(s, e4, e3, e1, e2)
                s%eri_table(idx) = val
                calculated(idx) = .true.
                idx = _hf_eri_index(s, e4, e3, e2, e1)
                s%eri_table(idx) = val
                calculated(idx) = .true.
            end do
            end do
            end do
            end do

            do e1=1,n
            do e2=1,n
            do e3=1,n
            do e4=1,n
                idx = _hf_eri_index(s, e1, e2, e3, e4)
                write (*,'(5(I5),E20.8)') idx, e1, e2, e3, e4, s%eri_table(idx)
            end do
            end do
            end do
            end do

            deallocate(calculated)
        end subroutine!}}}

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
