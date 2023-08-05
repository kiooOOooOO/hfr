module timer
    use, intrinsic :: iso_Fortran_env
    implicit none

    type timer_instance
        integer(int32) :: _tbegin, _tend, _cps, _cm
    end type

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
end module
