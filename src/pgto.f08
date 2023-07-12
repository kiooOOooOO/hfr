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

        subroutine _pgto_internal_division_point(ox, oy, oz, a, ax, ay, az, b, bx, by, bz)
            real(8), intent(out) :: ox, oy, oz
            real(8), intent(in) :: a, ax, ay, az, b, bx, by, bz

            ox = (a*ax + b*bx)/(a+b)
            oy = (a*ay + b*by)/(a+b)
            oz = (a*az + b*bz)/(a+b)
        end subroutine

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

        function pgto_eri(ga, gb, gc, gd) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb, gc, gd

            ret = ga%norm * gb%norm * gc%norm * gd%norm * _pgto_eri_internal(ga, gb, gc, gd, 0)
        end function

        recursive function _pgto_eri_internal(ga, gb, gc, gd, m) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb, gc, gd
            integer, intent(in) :: m

            real(8) :: am0m = 0, am0mp1 = 0
            real(8) :: am1m = 0, am1mp1 = 0
            real(8) :: bm1m = 0, bm1mp1 = 0
            real(8) :: cm1mp1 = 0
            real(8) :: dm1mp1 = 0

            real(8) :: px = 0, py = 0, pz = 0
            real(8) :: qx = 0, qy = 0, qz = 0
            real(8) :: wx = 0, wy = 0, wz = 0

            real(8) :: zeta = 0, eta = 0, rho = 0

            zeta = ga%expo + gb%expo
            eta  = gc%expo + gd%expo
            rho = zeta*eta/(zeta+eta)

            qx = (gc%expo*gc%cx + gd%expo*gd%cx)/eta
            qy = (gc%expo*gc%cy + gd%expo*gd%cy)/eta
            qz = (gc%expo*gc%cz + gd%expo*gd%cz)/eta

            call gto_internal_division_point(px, py, pz, ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)
            call gto_internal_division_point(wx, wy, wz, zeta, px, py, pz, eta, qx, qy, qz)

            if ( _pgto_all_n_zero(ga) .and. _pgto_all_n_zero(gb) .and. &
                 _pgto_all_n_zero(gc) .and. _pgto_all_n_zero(gd) ) then

                 ret = _pgto_eri_internal0(ga, gb, gc, gd, m)
             else if ( ga%nx .gt. 0 )  then
                 am0m = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, gc, gd, m)
                 am0mp1 = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, gc, gd, m+1)
                 if ( ga%nx .gt. 1 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, -2, 0, 0), gb, gc, gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, -2, 0, 0), gb, gc, gd, m+1)
                 end if
                 if ( gb%nx .gt. 0 ) then
                    bm1m = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), _pgto_clone(gb, -1, 0, 0), gc, gd, m)
                    bm1mp1 = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), _pgto_clone(gb, -1, 0, 0), gc, gd, m+1)
                 end if
                 if ( gc%nx .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, _pgto_clone(gc, -1, 0, 0), gd, m+1)
                 end if
                 if ( gd%nx .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, gc, _pgto_clone(gd, -1, 0, 0), m+1)
                 end if

                 ret = (px-ga%cx)*am0m + (wx-px)*am0mp1 + &
                     (0.5d0*(ga%nx-1)/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%nx/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%nx*cm1mp1 + gd%nx*dm1mp1)/(2d0*(zeta+eta))
             else if ( ga%ny .gt. 0 ) then
                 am0m = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, gc, gd, m)
                 am0mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, gc, gd, m+1)
                 if ( ga%ny .gt. 1 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, -2, 0), gb, gc, gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -2, 0), gb, gc, gd, m+1)
                 end if
                 if ( gb%ny .gt. 0 ) then
                    bm1m = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), _pgto_clone(gb, 0, -1, 0), gc, gd, m)
                    bm1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), _pgto_clone(gb, 0, -1, 0), gc, gd, m+1)
                 end if
                 if ( gc%ny .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, _pgto_clone(gc, 0, -1, 0), gd, m+1)
                 end if
                 if ( gd%ny .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, gc, _pgto_clone(gd, 0, -1, 0), m+1)
                 end if

                 ret = (py-ga%cy)*am0m + (wy-py)*am0mp1 + &
                     (0.5d0*(ga%ny-1)/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%ny/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%ny*cm1mp1 + gd%ny*dm1mp1)/(2d0*(zeta+eta))
             else if ( ga%nz .gt. 0 ) then
                 am0m = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, gc, gd, m)
                 am0mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, gc, gd, m+1)
                 if ( ga%ny .gt. 1 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -2), gb, gc, gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -2), gb, gc, gd, m+1)
                 end if
                 if ( gb%ny .gt. 0 ) then
                    bm1m = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), _pgto_clone(gb, 0, 0, -1), gc, gd, m)
                    bm1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), _pgto_clone(gb, 0, 0, -1), gc, gd, m+1)
                 end if
                 if ( gc%ny .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, _pgto_clone(gc, 0, 0, -1), gd, m+1)
                 end if
                 if ( gd%ny .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, gc, _pgto_clone(gd, 0, 0, -1), m+1)
                 end if

                 ret = (pz-ga%cz)*am0m + (wz-pz)*am0mp1 + &
                     (0.5d0*(ga%nz-1)/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%nz/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%nz*cm1mp1 + gd%nz*dm1mp1)/(2d0*(zeta+eta))
             else if ( gb%nx .gt. 0 )  then
                 am0m = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), gc, gd, m)
                 am0mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), gc, gd, m+1)
                 if ( ga%nx .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), _pgto_clone(gb, -1, 0, 0), gc, gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), _pgto_clone(gb, -1, 0, 0), gc, gd, m+1)
                 end if
                 if ( gb%nx .gt. 1 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, -2, 0, 0), gc, gd, m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, -2, 0, 0), gc, gd, m+1)
                 end if
                 if ( gc%nx .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), _pgto_clone(gc, -1, 0, 0), gd, m+1)
                 end if
                 if ( gd%nx .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), gc, _pgto_clone(gd, -1, 0, 0), m+1)
                 end if

                 ret = (px-gb%cx)*am0m + (wx-px)*am0mp1 + &
                     (0.5d0*ga%nx/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*(gb%nx-1)/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%nx*cm1mp1 + gd%nx*dm1mp1)/(2d0*(zeta+eta))
             else if ( gb%ny .gt. 0 ) then
                 am0m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), gc, gd, m)
                 am0mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), gc, gd, m+1)
                 if ( ga%ny .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), _pgto_clone(gb, 0, -1, 0), gc, gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), _pgto_clone(gb, 0, -1, 0), gc, gd, m+1)
                 end if
                 if ( gb%ny .gt. 1 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -2, 0), gc, gd, m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -2, 0), gc, gd, m+1)
                 end if
                 if ( gc%ny .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), _pgto_clone(gc, 0, -1, 0), gd, m+1)
                 end if
                 if ( gd%ny .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), gc, _pgto_clone(gd, 0, -1, 0), m+1)
                 end if

                 ret = (py-gb%cy)*am0m + (wy-py)*am0mp1 + &
                     (0.5d0*ga%ny/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*(gb%ny-1)/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%ny*cm1mp1 + gd%ny*dm1mp1)/(2d0*(zeta+eta))
             else if ( gb%nz .gt. 0 ) then
                 am0m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), gc, gd, m)
                 am0mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), gc, gd, m+1)
                 if ( ga%nz .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), _pgto_clone(gb, 0, 0, -1), gc, gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), _pgto_clone(gb, 0, 0, -1), gc, gd, m+1)
                 end if
                 if ( gb%nz .gt. 1 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -2), gc, gd, m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -2), gc, gd, m+1)
                 end if
                 if ( gc%nz .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), _pgto_clone(gc, 0, 0, -1), gd, m+1)
                 end if
                 if ( gd%nz .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), gc, _pgto_clone(gd, 0, 0, -1), m+1)
                 end if

                 ret = (pz-gb%cz)*am0m + (wz-pz)*am0mp1 + &
                     (0.5d0*ga%nz/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*(gb%nz-1)/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%nz*cm1mp1 + gd%nz*dm1mp1)/(2d0*(zeta+eta))
             else if ( gc%nx .gt. 0 )  then
                 am0m = _pgto_eri_internal(ga, gb, _pgto_clone(gc, -1, 0, 0), gd, m)
                 am0mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, -1, 0, 0), gd, m+1)
                 if ( ga%nx .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, _pgto_clone(gc, -1, 0, 0), gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, _pgto_clone(gc, -1, 0, 0), gd, m+1)
                 end if
                 if ( gb%nx .gt. 0 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), _pgto_clone(gc, -1, 0, 0), gd, m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), _pgto_clone(gc, -1, 0, 0), gd, m+1)
                 end if
                 if ( gc%nx .gt. 1 ) then
                     cm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, -2, 0, 0), gd, m+1)
                 end if
                 if ( gd%nx .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, -1, 0, 0), _pgto_clone(gd, -1, 0, 0), m+1)
                 end if

                 ret = (px-gc%cx)*am0m + (wx-px)*am0mp1 + &
                     (0.5d0*ga%nx/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%nx/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     ((gc%nx-1)*cm1mp1 + gd%nx*dm1mp1)/(2d0*(zeta+eta))
             else if ( gc%ny .gt. 0 ) then
                 am0m = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, -1, 0), gd, m)
                 am0mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, -1, 0), gd, m+1)
                 if ( ga%ny .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, _pgto_clone(gc, 0, -1, 0), gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, _pgto_clone(gc, 0, -1, 0), gd, m+1)
                 end if
                 if ( gb%ny .gt. 0 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), _pgto_clone(gc, 0, -1, 0), gd, m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), _pgto_clone(gc, 0, -1, 0), gd, m+1)
                 end if
                 if ( gc%ny .gt. 1 ) then
                     cm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, -2, 0), gd, m+1)
                 end if
                 if ( gd%ny .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, -1, 0), _pgto_clone(gd, 0, -1, 0), m+1)
                 end if

                 ret = (py-gc%cy)*am0m + (wy-py)*am0mp1 + &
                     (0.5d0*ga%ny/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%ny/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     ((gc%ny-1)*cm1mp1 + gd%ny*dm1mp1)/(2d0*(zeta+eta))
             else if ( gc%nz .gt. 0 ) then
                 am0m = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, 0, -1), gd, m)
                 am0mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, 0, -1), gd, m+1)
                 if ( ga%nz .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, _pgto_clone(gc, 0, 0, -1), gd, m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, _pgto_clone(gc, 0, 0, -1), gd, m+1)
                 end if
                 if ( gb%nz .gt. 0 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), _pgto_clone(gc, 0, 0, -1), gd, m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), _pgto_clone(gc, 0, 0, -1), gd, m+1)
                 end if
                 if ( gc%nz .gt. 1 ) then
                     cm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, 0, -2), gd, m+1)
                 end if
                 if ( gd%nz .gt. 0 ) then
                     dm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, 0, -1), _pgto_clone(gd, 0, 0, -1), m+1)
                 end if

                 ret = (pz-gc%cz)*am0m + (wz-pz)*am0mp1 + &
                     (0.5d0*ga%nz/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%nz/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     ((gc%nz-1)*cm1mp1 + gd%nz*dm1mp1)/(2d0*(zeta+eta))
             else if ( gd%nx .gt. 0 )  then
                 am0m = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, -1, 0, 0), m)
                 am0mp1 = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, -1, 0, 0), m+1)
                 if ( ga%nx .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, gc, _pgto_clone(gd, -1, 0, 0), m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, -1, 0, 0), gb, gc, _pgto_clone(gd, -1, 0, 0), m+1)
                 end if
                 if ( gb%nx .gt. 0 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), gc, _pgto_clone(gd, -1, 0, 0), m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, -1, 0, 0), gc, _pgto_clone(gd, -1, 0, 0), m+1)
                 end if
                 if ( gc%nx .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, -1, 0, 0), _pgto_clone(gd, -1, 0, 0), m+1)
                 end if
                 if ( gd%nx .gt. 1 ) then
                     dm1mp1 = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, -2, 0, 0), m+1)
                 end if

                 ret = (px-gd%cx)*am0m + (wx-px)*am0mp1 + &
                     (0.5d0*ga%nx/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%nx/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%nx*cm1mp1 + (gd%nx-1)*dm1mp1)/(2d0*(zeta+eta))
             else if ( gd%ny .gt. 0 ) then
                 am0m = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, 0, -1, 0), m)
                 am0mp1 = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, 0, -1, 0), m+1)
                 if ( ga%ny .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, gc, _pgto_clone(gd, 0, -1, 0), m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, -1, 0), gb, gc, _pgto_clone(gd, 0, -1, 0), m+1)

                    
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), gc, _pgto_clone(gd, 0, -1, 0), m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, -1, 0), gc, _pgto_clone(gd, 0, -1, 0), m+1)
                 end if
                 if ( gc%ny .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, -1, 0), _pgto_clone(gd, 0, -1, 0), m+1)
                 end if
                 if ( gd%ny .gt. 1 ) then
                     dm1mp1 = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, 0, -2, 0), m+1)
                 end if

                 ret = (py-gd%cy)*am0m + (wy-py)*am0mp1 + &
                     (0.5d0*ga%ny/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%ny/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%ny*cm1mp1 + (gd%ny-1)*dm1mp1)/(2d0*(zeta+eta))
             else if ( gc%nz .gt. 0 ) then
                 am0m = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, 0, 0, -1), m)
                 am0mp1 = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, 0, 0, -1), m+1)
                 if ( ga%nz .gt. 0 ) then
                    am1m = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, gc, _pgto_clone(gd, 0, 0, -1), m)
                    am1mp1 = _pgto_eri_internal(_pgto_clone(ga, 0, 0, -1), gb, gc, _pgto_clone(gd, 0, 0, -1), m+1)
                 end if
                 if ( gb%nz .gt. 0 ) then
                    bm1m = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), gc, _pgto_clone(gd, 0, 0, -1), m)
                    bm1mp1 = _pgto_eri_internal(ga, _pgto_clone(gb, 0, 0, -1), gc, _pgto_clone(gd, 0, 0, -1), m+1)
                 end if
                 if ( gc%nz .gt. 0 ) then
                     cm1mp1 = _pgto_eri_internal(ga, gb, _pgto_clone(gc, 0, 0, -1), _pgto_clone(gd, 0, 0, -1), m+1)
                 end if
                 if ( gd%nz .gt. 1 ) then
                     dm1mp1 = _pgto_eri_internal(ga, gb, gc, _pgto_clone(gd, 0, 0, -2), m+1)
                 end if

                 ret = (pz-gd%cz)*am0m + (wz-pz)*am0mp1 + &
                     (0.5d0*ga%nz/zeta)*(am1m - rho*am1mp1/zeta) + &
                     (0.5d0*gb%nz/zeta)*(bm1m - rho*bm1mp1/zeta) + &
                     (gc%nz*cm1mp1 + (gd%nz-1)*dm1mp1)/(2d0*(zeta+eta))
             else
                 write (*,*) "warning"
                 stop
             end if

        end function

        function _pgto_eri_internal0(ga, gb, gc, gd, m) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb, gc, gd
            integer, intent(in) :: m

            real(8) :: zeta = 0, eta = 0, rho = 0
            real(8) :: qx = 0, qy = 0, qz = 0
            real(8) :: px = 0, py = 0, pz = 0
            real(8) :: t = 0

            zeta = ga%expo + gb%expo
            eta  = gc%expo + gd%expo
            rho = zeta*eta/(zeta+eta)

            qx = (gc%expo*gc%cx + gd%expo*gd%cx)/eta
            qy = (gc%expo*gc%cy + gd%expo*gd%cy)/eta
            qz = (gc%expo*gc%cz + gd%expo*gd%cz)/eta


            call _pgto_internal_division_point(px, py, pz, ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)

            t = (px-qx)**2 + (py-qy)**2 + (pz-qz)**2
            t = rho * t



            ret = (zeta + eta)**-0.5d0 * _pgto_kappa(ga, gb) * _pgto_kappa(gc, gd) * _pgto_fm(t, m)
        end function

        function _pgto_kappa(ga, gb) result(kappa)
            real(8), intent(out) :: kappa
            type(pgto), intent(in) :: ga, gb

            real(8) :: d2

            d2 = (ga%cx-gb%cx)**2 + (ga%cy-gb%cy)**2 + (ga%cz-gb%cz)**2

            kappa = 2d0**0.5d0 * PI**1.25d0 * exp(-d2*ga%expo*gb%expo/(ga%expo + gb%expo)) / (ga%expo + gb%expo)
        end function

        function _pgto_fm(t, m) result(fm)
            real(8), intent(out) :: fm
            real(8), intent(in) :: t
            integer, intent(in) :: m

            real(8), parameter :: DV = 1e-5
            real(8) :: sumation = 0, v = 0

            v = 0
            sumation = 0
            do while ( v .le. 1 )
                sumation = sumation + 0.5d0 * DV * (_pgto_fm_func(t, m, v) + _pgto_fm_func(t, m, v+DV))
                v = v + DV
            end do

            fm = sumation
        end function

        function _pgto_fm_func(t, m, v) result(func)
            real(8), intent(out) :: func
            real(8), intent(in) :: t, v
            integer, intent(in) :: m

            func = v**(2*m) * exp(-t*v**2)
        end function
end module
