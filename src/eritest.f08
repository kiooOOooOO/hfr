program eritest_main
    use, intrinsic :: iso_fortran_env
    use pgto
    use pgto2
    implicit none

    integer, parameter :: N = 10

    type timer_instance
        integer(int32) :: _tbegin, _tend, _cps, _cm
    end type

    integer, parameter :: CASE_COUNT = 7
    real(8), parameter, dimension(4*7+1,CASE_COUNT) :: C = reshape( (/ &
        ! expo, x,y,z, nx,ny,nz

        ! s orbital only
        0.1d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.2d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.3d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.4d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.4661054519858802d0, &
        0.1d0, 4d0,0d0,0d0, 0d0,0d0,0d0, &
        0.2d0, 4d0,0d0,0d0, 0d0,0d0,0d0, &
        0.3d0, 4d0,0d0,0d0, 0d0,0d0,0d0, &
        0.4d0, 4d0,0d0,0d0, 0d0,0d0,0d0, &
        0.4661054519858802d0, &
        0.1d0, 0d0,4d0,0d0, 0d0,0d0,0d0, &
        0.2d0, 0d0,4d0,0d0, 0d0,0d0,0d0, &
        0.3d0, 0d0,4d0,0d0, 0d0,0d0,0d0, &
        0.4d0, 0d0,4d0,0d0, 0d0,0d0,0d0, &
        0.4661054519858802d0, &
        0.1d0, 0d0,0d0,4d0, 0d0,0d0,0d0, &
        0.2d0, 0d0,0d0,4d0, 0d0,0d0,0d0, &
        0.3d0, 0d0,0d0,4d0, 0d0,0d0,0d0, &
        0.4d0, 0d0,0d0,4d0, 0d0,0d0,0d0, &
        0.4661054519858802d0, &

        ! (px, s, s, s)
        0.1d0, 2d0,0d0,0d0, 1d0,0d0,0d0, &
        0.2d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.3d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.4d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        -0.5324180896198726d0, &

        ! (py, s, s, s)
        0.1d0, 0d0,2d0,0d0, 0d0,1d0,0d0, &
        0.2d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.3d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.4d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        -0.5324180896198726d0, &

        ! (pz, s, s, s)
        0.1d0, 0d0,0d0,2d0, 0d0,0d0,1d0, &
        0.2d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.3d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        0.4d0, 0d0,0d0,0d0, 0d0,0d0,0d0, &
        -0.5324180896198726d0, &
        2525d0 /), (/29,CASE_COUNT/) )
    real(8), parameter :: EPS = 1d-4

    integer :: i, j, k, base
    type(pgto), dimension(4) :: g
    type(timer_instance) :: ti1, ti2
    real(8) :: val, msec, actual, expected

    real(8) :: cx, cy, cz, expo
    integer :: nx, ny, nz

    do i=1,CASE_COUNT
        do k=1,4
            base = (k-1)*7

            expo = C(base+1,i)
            cx = C(base+2,i)
            cy = C(base+3,i)
            cz = C(base+4,i)
            nx = C(base+5,i)
            ny = C(base+6,i)
            nz = C(base+7,i)

!            write (*,'(4(E15.5), 3(I3))') expo, cx, cy, cz, nx, ny, nz

            g(k) = pgto_new(expo, cx,cy,cz, nx,ny,nz)
        end do

        actual = pgto2_eri(g(1), g(2), g(3), g(4))
        expected = C(29,i)

        call timer_start(ti1)
        do j=1,N
            val = pgto2_eri(g(1), g(2), g(3), g(4))
        end do
        call timer_stop(ti1, msec)
        msec = 1000*msec/N

        write (*,'("Case #", I0, F15.12, F15.3, " msec")') i, actual, msec
        if ( abs(expected-actual) .gt. EPS ) then
            write (*,*) "wrong value", expected
            val = pgto_eri(g(1), g(2), g(3), g(4))
            write (*,*) val
        end if
    end do

    val = _pgto_fm(7d0/75d0, 0)
    write (*,*) "F0(7/75)", val
    val = _pgto_fm(7d0/75d0, 1)
    write (*,*) "F1(7/75)", val

!    g1 = pgto_new(0.1611957475d2, 0d0,10d0,0d0, 1,0,0)
!    g2 = pgto_new(3.425250914d0,  -1.3228d0,10d0,0d0, 1,0,0)
!    g3 = pgto_new(13.425250914d0,  -1.3228d0,10d0,0d0, 1,0,0)
!    g4 = pgto_new(53.425250914d0,  -1.3228d0,10d0,0d0, 1,0,0)
!
!    call timer_start(ti1)
!    do i=1,N
!        val = pgto2_eri(g1,g2,g3,g4)
!    end do
!    call timer_stop(ti1, sec)
!    write (*,*) val
!    write (*,*) 1000*sec/N, "msec"
!
!    g1 = pgto_new(0.1611957475d2, 0d0,10d0,0d0, 0,0,1)
!    g2 = pgto_new(3.425250914d0,  -1.3228d0,10d0,0d0, 0,0,1)
!    g3 = pgto_new(13.425250914d0,  -1.3228d0,10d0,0d0, 0,0,1)
!    g4 = pgto_new(53.425250914d0,  -1.3228d0,10d0,0d0, 0,0,1)
!
!    call timer_start(ti2)
!    do i=1,N
!        val = pgto2_eri(g1,g2,g3,g4)
!    end do
!    call timer_stop(ti2, sec)
!    write (*,*) val
!    write (*,*) 1000*sec/N, "msec"

    contains

        subroutine timer_start(ti)
            type(timer_instance) :: ti
            call system_clock(ti%_tbegin, ti%_cps, ti%_cm)
        end subroutine

        subroutine timer_stop(ti, r)
            type(timer_instance) :: ti
            real(8), intent(out) :: r
            call system_clock(ti%_tend)
            r = real(ti%_tend-ti%_tbegin)/ti%_cps
        end subroutine
end program
