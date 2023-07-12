module pgto
    implicit none

    real(8), parameter :: PI = 3.141592653d0

    type pgto
        real(8) :: cx, cy, cz
        real(8) :: expo
        real(8) :: norm
        integer :: nx, ny, nz
    end type

    contains

        pure function pgto_new(expo, cx, cy, cz, nx, ny, nz)
            type(pgto), intent(out) :: pgto_new
            real(8), intent(in) :: cx, cy, cz
            real(8), intent(in) :: expo
            integer, intent(in) :: nx, ny, nz

            pgto_new%expo = expo
            pgto_new%cx = cx
            pgto_new%cy = cy
            pgto_new%cz = cz
            pgto_new%nx = nx
            pgto_new%ny = ny
            pgto_new%nz = nz

            pgto_new%norm = _pgto_norm(expo, nx, ny, nz)
        end function

        pure function _pgto_norm(expo, nx, ny, nz) result(norm)
            real(8), intent(out) :: norm
            real(8), intent(in) :: expo
            integer, intent(in) :: nx, ny, nz

            norm = (2d0*expo/PI)**0.75d0 * (4d0*expo)**(0.5d0*(nx+ny+nz)) &
                * _pgto_odd_factorial(2*nx-1) * _pgto_odd_factorial(2*ny-1) * _pgto_odd_factorial(2*nz-1)
        end function

        pure recursive function _pgto_odd_factorial(i) result(ret)
            real(8), intent(out) :: ret
            integer, intent(in) :: i

            if ( i .le. 1 ) then
                ret = 1d0
            else
                ret = ret * _pgto_odd_factorial(i-2)
            end if
        end function

        pure function _pgto_all_n_zero(g) result(ret)
            logical, intent(out) :: ret
            type(pgto), intent(in) :: g

            ret = (g%nx .eq. 0) .and. (g%ny .eq. 0) .and. (g%nz .eq. 0)
        end function

        function _pgto_clone(g, dnx, dny, dnz) result(ret)
            type(pgto) :: ret
            type(pgto) :: g
            integer :: dnx, dny, dnz

            ret%expo = g%expo
            ret%cx = g%cx
            ret%cy = g%cy
            ret%cz = g%cz
            ret%nx = g%nx + dnx
            ret%ny = g%ny + dny
            ret%nz = g%nz + dnz
        end function

        recursive function pgto_overlap(ga, gb) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb

            real(8) :: ab, am1b, abm1, zeta, px, py, pz

            if ( _pgto_all_n_zero(ga) .and. _pgto_all_n_zero(gb) ) then
                ret = _pgto_overlap0(ga, gb)
            else
                zeta = ga%expo + gb%expo
                ab = 0
                am1b = 0
                abm1 = 0

                call _pgto_internal_division_point(px, py, pz, ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)

                if ( ga%nx .gt. 0 ) then
                    ab = pgto_overlap(_pgto_clone(ga, -1, 0, 0), gb)
                    if ( ga%nx .gt. 1 ) then
                        am1b = pgto_overlap(_pgto_clone(ga, -2, 0, 0), gb)
                    end if
                    if ( gb%nx .gt. 0 ) then
                        abm1 = pgto_overlap(_pgto_clone(ga, -1, 0, 0), _pgto_clone(gb, -1, 0, 0))
                    end if
                    ret = (px-ga%cx)*ab + 0.5d0*(ga%nx-1)*am1b/zeta + 0.5d0*gb%nx*abm1/zeta
                else if ( ga%ny .gt. 0 ) then
                    ab = pgto_overlap(_pgto_clone(ga, 0, -1, 0), gb)
                    if ( ga%ny .gt. 1 ) then
                        am1b = pgto_overlap(_pgto_clone(ga, 0, -2, 0), gb)
                    end if
                    if ( gb%ny .gt. 0 ) then
                        abm1 = pgto_overlap(_pgto_clone(ga, 0, -1, 0), _pgto_clone(gb, 0, -1, 0))
                    end if
                    ret = (py-ga%cy)*ab + 0.5d0*(ga%ny-1)*am1b/zeta + 0.5d0*gb%ny*abm1/zeta
                else if ( ga%nz .gt. 0 ) then
                    ab = pgto_overlap(_pgto_clone(ga, 0, 0, -1), gb)
                    if ( ga%nz .gt. 1 ) then
                        am1b = pgto_overlap(_pgto_clone(ga, 0, 0, -2), gb)
                    end if
                    if ( gb%nz .gt. 0 ) then
                        abm1 = pgto_overlap(_pgto_clone(ga, 0, 0, -1), _pgto_clone(gb, 0, 0, -1))
                    end if
                    ret = (pz-ga%cz)*ab + 0.5d0*(ga%nz-1)*am1b/zeta + 0.5d0*gb%nz*abm1/zeta
                else
                    ret = pgto_overlap(gb, ga)
                end if
            end if
        end function

        function _pgto_overlap0(ga, gb) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb

            real(8) :: gzi, zeta, r2

            r2 = (ga%cx-gb%cx)**2 + (ga%cy-gb%cy)**2 + (ga%cz-gb%cz)**2
            zeta = ga%expo + gb%expo
            gzi = ga%expo*gb%expo/zeta

            ret = ga%norm * gb%norm * (PI/zeta)**1.5d0 * exp(-gzi*r2)
        end function

        subroutine _pgto_internal_division_point(ox, oy, oz, a, ax, ay, az, b, bx, by, bz)
            real(8), intent(out) :: ox, oy, oz
            real(8), intent(in) :: a, ax, ay, az, b, bx, by, bz

            ox = (a*ax + b*bx)/(a+b)
            oy = (a*ay + b*by)/(a+b)
            oz = (a*az + b*bz)/(a+b)
        end subroutine
end module
