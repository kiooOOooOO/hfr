module matrix
    implicit none

    real(8), parameter :: PI = 3.1415926
    real(8), parameter :: MATRIX_ZEROFY_EPS = 1e-10

    contains

        ! matA ... i0行j0列 行列
        ! matB ... j0行k0列 行列
        ! matR ... i0行k0列 行列
        ! matR <- matA x matB
        subroutine matrix_mult_normal_normal(matA, matB, i0, j0, k0, matR)
            real(8), dimension(:,:) :: matA, matB, matR
            integer, intent(in) :: i0, j0, k0
            call dgemm('n', 'n', i0, k0, j0, 1d0, matA, i0, matB, j0, 0d0, matR, i0)
        end subroutine

        ! matA ... i0行j0列 行列
        ! matB ... i0行k0列 行列
        ! matR ... j0行k0列 行列
        ! matR <- matA_T x matB
        subroutine matrix_mult_transpose_normal(matA, matB, i0, j0, k0, matR)
            real(8), dimension(:,:) :: matA, matB, matR
            integer, intent(in) :: i0, j0, k0
            call dgemm('t', 'n', j0, k0, i0, 1d0, matA, i0, matB, i0, 0d0, matR, j0)
        end subroutine

        ! matA ... i0行j0列 行列
        ! matB ... k0行j0列 行列
        ! matR ... i0行k0列 行列
        ! matR <- matA_T x matB
        subroutine matrix_mult_normal_transpose(matA, matB, i0, j0, k0, matR)
            real(8), dimension(:,:) :: matA, matB, matR
            integer, intent(in) :: i0, j0, k0
            call dgemm('n', 't', i0, k0, j0, 1d0, matA, i0, matB, k0, 0d0, matR, i0)
        end subroutine

        ! 行列 A に単位行列を代入する．
        ! 行列 A は nlen 次製法行列である．
        subroutine matrix_set_elemental(matA, nlen)
            real(8), dimension(:,:) :: matA
            integer :: nlen

            integer :: i

            matA = 0
            do i=1,nlen
                matA(i,i) = 1d0
            end do
        end subroutine

        ! 行列 A の逆行列を B に代入する．
        ! 行列 A, B はそれぞれ nlen 次正方行列である．
        subroutine matrix_a_inverse(matA, matB, nlen)
            real(8), dimension(:,:) :: matA, matB
            integer :: nlen

            real(8), allocatable, dimension(:,:) :: matL, matU, matLinv, matUinv
            allocate(matL(nlen,nlen))
            allocate(matU(nlen,nlen))
            allocate(matLinv(nlen,nlen))
            allocate(matUinv(nlen,nlen))

            call matrix_lu_disasm(matA, matL, matU, nlen)
            call matrix_set_elemental(matLinv, nlen)
            call matrix_set_elemental(matUinv, nlen)

            call dtrsm("l", "l", "n", "n", nlen, nlen, 1d0, matL, nlen, matLinv, nlen)
            call dtrsm("l", "u", "n", "n", nlen, nlen, 1d0, matU, nlen, matUinv, nlen)

            call dgemm("n", "n", nlen, nlen, nlen, 1d0, matUinv, nlen, matLinv, nlen, 0d0, matB, nlen)

            deallocate(matL)
            deallocate(matU)
            deallocate(matLinv)
            deallocate(matUinv)
        end subroutine

        ! 行列 A を対角化する
        ! 行列 A は nlen 次正方対象行列
        ! 行列 R は対角化された行列
        ! 行列 U は対角化行列
        subroutine matrix_sym_diagonalize(matA, nlen, matR, matU)
            real(8), dimension(:,:) :: matA, matR, matU
            integer :: nlen

            real(8), allocatable, dimension(:) :: evals
            real(8), allocatable, dimension(:,:) :: evecs
            integer, allocatable, dimension(:) :: order
            integer :: i, j, smallest
            real(8), allocatable, dimension(:) :: tmparr
            real(8) :: tmp

            allocate(tmparr(nlen))
            allocate(evals(nlen))
            allocate(evecs(nlen,nlen))

            call matrix_sym_eigenvalue(matA, nlen, evals, evecs)

            do i=1,nlen-1
                smallest = i
                do j=i+1,nlen
                    if ( evals(j) .lt. evals(smallest) ) then
                        smallest = j
                    end if
                end do
                tmp = evals(i)
                evals(i) = evals(smallest)
                evals(smallest) = tmp

                tmparr = evecs(i,1:nlen)
                evecs(i,1:nlen) = evecs(smallest,1:nlen)
                evecs(smallest,1:nlen) = tmparr
            end do

            matR = 0
            do i=1,nlen
               matR(i,i) = evals(i)
            do j=1,nlen
               matU(i,j) = evecs(j,i)
!               matR(nlen-i+1,nlen-i+1) = evals(i)
!               matU(i,nlen-j+1) = evecs(j,i)
            end do
            end do

            deallocate(tmparr)
            deallocate(evals)
            deallocate(evecs)
        end subroutine

        ! 行列 A の固有値を計算する
        ! 行列 A は nlen 次正方対象行列
        ! 固有値 evals(i) に対応する固有ベクトルは evecs(i, 1:nlen)
        subroutine matrix_sym_eigenvalue(matA, nlen, evals, evecs)
            real(8), dimension(:,:) :: matA
            integer :: nlen
            real(8), dimension(:) :: evals
            real(8), dimension(:,:) :: evecs

            real(8), allocatable, dimension(:,:) :: matT, matE

            integer :: i, j, k
            integer :: p, q
            real(8) :: maxValue, theta, t1, t2

            allocate(matT(nlen,nlen))
            allocate(matE(nlen,nlen))
            matE = 0d0
            do i=1,nlen
            do j=1,nlen
                matT(i,j) = matA(i,j)
            end do
            matE(i,i) = 1d0
            end do

            maxValue = 1
            do while ( maxValue .gt. 1e-5 )
                p = 1
                q = 1
                maxValue = 0
                do i=1,nlen
                do j=i+1,nlen
                    if ( abs(matT(i,j)) .gt. maxValue ) then
                        maxValue = abs(matT(i,j))
                        p = i
                        q = j
                    end if
                end do
                end do

                theta = 0d0
                if ( abs(matT(p,p) - matT(q,q)) .lt. 1e-5 ) then
                    theta = PI / 4.0d0
                    if ( matT(p,p) .lt. 0 ) then
                        theta = -theta
                    end if
                else
                    theta = 0.5 * atan(2.0*matT(p,q) / (matT(p,p)-matT(q,q)))
                end if

                do i=1,nlen
                    t1 = matT(p,i)*cos(theta) + matT(q,i)*sin(theta)
                    t2 = -matT(p,i)*sin(theta) + matT(q,i)*cos(theta)
                    matT(p,i) = t1
                    matT(q,i) = t2

                    t1 = matE(p,i)*cos(theta) + matE(q,i)*sin(theta)
                    t2 = -matE(p,i)*sin(theta) + matE(q,i)*cos(theta)
                    matE(p,i) = t1
                    matE(q,i) = t2
                end do

                do i=1,nlen
                    t1 = matT(i,p)*cos(theta) + matT(i,q)*sin(theta)
                    t2 = -matT(i,p)*sin(theta) + matT(i,q)*cos(theta)
                    matT(i,p) = t1
                    matT(i,q) = t2

                end do
            end do

!            write (*,*) ">", matT
!            write (*,*) ">>", matE

            do i=1,nlen
                evals(i) = matT(i,i)
            end do

            do i=1,nlen
            do j=1,nlen
                evecs(i,j) = matE(i,j)
            end do
            end do

            deallocate(matT)
            deallocate(matE)
        end subroutine

        ! 行列 A に対して LU 分解を行う． L が下三角行列となる．
        ! 行列 A, L, U はそれぞれ nlen 次正方行列である．
        subroutine matrix_lu_disasm(matA, matL, matU, nlen)
            real(8), dimension(:,:) :: matA, matL, matU
            integer :: nlen

            real(8) :: r
            integer :: i, j, k

            matL = 0
            matU = 0

            do j=1,nlen
                do i=1,j
                    if ( i .eq. 1 ) then
                        matU(i,j) = matA(i,j)
                    else
                        r = matA(i,j)
                        do k=1,i-1
                            r = r - matL(i,k)*matU(k,j)
                        end do
                        matU(i,j) = r
                    end if
                end do
                do i=j+1,nlen
                    r = matA(i,j)
                    do k=1,j-1
                        r = r - matL(i,k)*matU(k,j)
                    end do
                    if ( is_zero(matU(j,j)) ) then
                        write (*,*) "div zero found"
                    end if
                    matL(i,j) = r/matU(j,j)
                end do

                matL(j,j) = 1
            end do
        end subroutine

        function is_zero(r)
            logical :: is_zero
            real(8) :: r

            is_zero = abs(r) .lt. 1e-10
        end function

        ! 行列 A のうち 0.0 であろう項を本当に 0.0 にする
        ! 行列 A は nlen 次正方行列である．
        subroutine matrix_zerofy(matA, nlen, eps)
            real(8), dimension(:,:) :: matA
            integer :: nlen
            real(8) :: eps

            integer :: i, j

            do i=1,nlen
            do j=1,nlen
                if ( abs(matA(i,j)) .lt. eps ) then
                    matA(i,j) = 0d0
                end if
            end do
            end do
        end subroutine

end module
