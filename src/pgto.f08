module pgto
    implicit none

    real(8), parameter :: PI = 3.141592653d0

    type pgto
        real(8) :: cx, cy, cz
        real(8) :: expo
        real(8) :: norm
        integer :: nx, ny, nz
    end type

    type vec3d
        real(8) :: x, y, z
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

        ! NOTE
        ! cloned pgto won't calculate normalization factor
        pure function _pgto_clone(g, dnx, dny, dnz) result(ret)
            type(pgto), intent(out) :: ret
            type(pgto), intent(in) :: g
            integer, intent(in) :: dnx, dny, dnz

            ret%expo = g%expo
            ret%cx = g%cx
            ret%cy = g%cy
            ret%cz = g%cz
            ret%nx = g%nx + dnx
            ret%ny = g%ny + dny
            ret%nz = g%nz + dnz
            ret%norm = g%norm
        end function

        pure function _pgto_internal_division_point(a, ax, ay, az, b, bx, by, bz) result(ret)
            type(vec3d), intent(out) :: ret
            real(8), intent(in) :: a, ax, ay, az, b, bx, by, bz

            ret%x = (a*ax + b*bx)/(a+b)
            ret%y = (a*ay + b*by)/(a+b)
            ret%z = (a*az + b*bz)/(a+b)
        end function

        recursive function pgto_overlap(ga, gb) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb

            real(8) :: ab, am1b, abm1, zeta, px, py, pz
            type(vec3d) :: point

            if ( _pgto_all_n_zero(ga) .and. _pgto_all_n_zero(gb) ) then
                ret = _pgto_overlap0(ga, gb)
            else
                zeta = ga%expo + gb%expo
                ab = 0
                am1b = 0
                abm1 = 0

                point = _pgto_internal_division_point(ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)
                px = point%x
                py = point%y
                pz = point%z

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

        ! ∬ dr1 dr2 ga(1)gb(1)gc(2)gd(2)/r12
        pure function pgto_eri(ga, gb, gc, gd) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb, gc, gd

            ret = ga%norm * gb%norm * gc%norm * gd%norm * _pgto_eri_internal(ga, gb, gc, gd, 0)
        end function

        recursive pure function _pgto_eri_internal(pga, pgb, pgc, pgd, m) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: pga, pgb, pgc, pgd
            integer, intent(in) :: m

            type(pgto) :: ga, gb, gc, gd
            real(8) :: px, py, pz
            real(8) :: qx, qy, qz
            real(8) :: wx, wy, wz
            real(8) :: zeta, eta, rho
            real(8) :: am, amp1, am1m, am1mp1, bm1m, bm1mp1, cm1mp1, dm1mp1
            logical :: za, zb, zc, zd

            type(vec3d) :: point

            ga = _pgto_clone(pga, 0,0,0)
            gb = _pgto_clone(pgb, 0,0,0)
            gc = _pgto_clone(pgc, 0,0,0)
            gd = _pgto_clone(pgd, 0,0,0)
            za = _pgto_all_n_zero(ga)
            zb = _pgto_all_n_zero(gb)
            zc = _pgto_all_n_zero(gc)
            zd = _pgto_all_n_zero(gd)

            if ( za .and. zb .and. zc .and. zd ) then
                ret = _pgto_eri_internal0(ga, gb, gc, gd, m)
            else if ( za .and. ( .not. zb ) ) then
                ret = _pgto_eri_internal(gb, ga, gc, gd, m)
            else if ( (za .and. zb) .and. ((.not. zc) .or. (.not. zd)) ) then
                ret = _pgto_eri_internal(gc, gd, ga, gb, m)
            else
                zeta = ga%expo + gb%expo
                eta  = gc%expo + gd%expo
                rho = zeta*eta/(zeta+eta)

                qx = (gc%expo*gc%cx + gd%expo*gd%cx)/eta
                qy = (gc%expo*gc%cy + gd%expo*gd%cy)/eta
                qz = (gc%expo*gc%cz + gd%expo*gd%cz)/eta

                point = _pgto_internal_division_point(ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)
                px = point%x
                py = point%y
                pz = point%z
                point = _pgto_internal_division_point(zeta, px, py, pz, eta, qx, qy, qz)
                wx = point%x
                wy = point%y
                wz = point%z

                if ( ga%nx .gt. 0 ) then
                    ga%nx = ga%nx - 1
                    am = _pgto_eri_internal(ga, gb, gc, gd, m)
                    amp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)

                    if ( ga%nx .gt. 0 ) then
                        ga%nx = ga%nx - 1
                        am1m = _pgto_eri_internal(ga, gb, gc, gd, m)
                        am1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        ga%nx = ga%nx + 1
                    else
                        am1m = 0d0
                        am1mp1 = 0d0
                    end if

                    if ( gb%nx .gt. 0 ) then
                        gb%nx = gb%nx - 1
                        bm1m = _pgto_eri_internal(ga, gb, gc, gd, m)
                        bm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gb%nx = gb%nx + 1
                    else
                        bm1m = 0d0
                        bm1mp1 = 0d0
                    end if

                    if ( gc%nx .gt. 0 ) then
                        gc%nx = gc%nx - 1
                        cm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gc%nx = gc%nx + 1
                    else
                        cm1mp1 = 0d0
                    end if

                    if ( gd%nx .gt. 0 ) then
                        gd%nx = gd%nx - 1
                        dm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gd%nx = gd%nx + 1
                    else
                        dm1mp1 = 0d0
                    end if

                    ret = (px-ga%cx)*am + (wx-px)*amp1 + &
                        0.5d0 * (am1m - rho*am1mp1/zeta) * ga%nx/zeta + &
                        0.5d0 * (bm1m - rho*bm1mp1/zeta) * gb%nx/zeta + &
                        0.5d0 * cm1mp1/(zeta+eta) * gc%nx + &
                        0.5d0 * dm1mp1/(zeta+eta) * gd%nx
                    ga%nx = ga%nx + 1
                else if ( ga%ny .gt. 0) then
                    ga%ny = ga%ny - 1
                    am = _pgto_eri_internal(ga, gb, gc, gd, m)
                    amp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)

                    if ( ga%ny .gt. 0 ) then
                        ga%ny = ga%ny - 1
                        am1m = _pgto_eri_internal(ga, gb, gc, gd, m)
                        am1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        ga%ny = ga%ny + 1
                    else
                        am1m = 0d0
                        am1mp1 = 0d0
                    end if

                    if ( gb%ny .gt. 0 ) then
                        gb%ny = gb%ny - 1
                        bm1m = _pgto_eri_internal(ga, gb, gc, gd, m)
                        bm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gb%ny = gb%ny + 1
                    else
                        bm1m = 0d0
                        bm1mp1 = 0d0
                    end if

                    if ( gc%ny .gt. 0 ) then
                        gc%ny = gc%ny - 1
                        cm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gc%ny = gc%ny + 1
                    else
                        cm1mp1 = 0d0
                    end if

                    if ( gd%ny .gt. 0 ) then
                        gd%ny = gd%ny - 1
                        dm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gd%ny = gd%ny + 1
                    else
                        dm1mp1 = 0d0
                    end if

                    ret = (py-ga%cy)*am + (wy-py)*amp1 + &
                        0.5d0 * (am1m - rho*am1mp1/zeta) * ga%ny/zeta + &
                        0.5d0 * (bm1m - rho*bm1mp1/zeta) * gb%ny/zeta + &
                        0.5d0 * cm1mp1/(zeta+eta) * gc%ny + &
                        0.5d0 * dm1mp1/(zeta+eta) * gd%ny
                    ga%ny = ga%ny + 1
                else if ( ga%nz .gt. 0) then
                    ga%nz = ga%nz - 1
                    am = _pgto_eri_internal(ga, gb, gc, gd, m)
                    amp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)

                    if ( ga%nz .gt. 0 ) then
                        ga%nz = ga%nz - 1
                        am1m = _pgto_eri_internal(ga, gb, gc, gd, m)
                        am1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        ga%nz = ga%nz + 1
                    else
                        am1m = 0d0
                        am1mp1 = 0d0
                    end if

                    if ( gb%nz .gt. 0 ) then
                        gb%nz = gb%nz - 1
                        bm1m = _pgto_eri_internal(ga, gb, gc, gd, m)
                        bm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gb%nz = gb%nz + 1
                    else
                        bm1m = 0d0
                        bm1mp1 = 0d0
                    end if

                    if ( gc%nz .gt. 0 ) then
                        gc%nz = gc%nz - 1
                        cm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gc%nz = gc%nz + 1
                    else
                        cm1mp1 = 0d0
                    end if

                    if ( gd%nz .gt. 0 ) then
                        gd%nz = gd%nz - 1
                        dm1mp1 = _pgto_eri_internal(ga, gb, gc, gd, m+1)
                        gd%nz = gd%nz + 1
                    else
                        dm1mp1 = 0d0
                    end if

                    ret = (pz-ga%cz)*am + (wz-pz)*amp1 + &
                        0.5d0 * (am1m - rho*am1mp1/zeta) * ga%nz/zeta + &
                        0.5d0 * (bm1m - rho*bm1mp1/zeta) * gb%nz/zeta + &
                        0.5d0 * cm1mp1/(zeta+eta) * gc%nz + &
                        0.5d0 * dm1mp1/(zeta+eta) * gd%nz
                    ga%nz = ga%nz + 1
                else
                    ret = 25252525d0
                end if
            end if
        end function

        pure function _pgto_eri_internal0(ga, gb, gc, gd, m) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb, gc, gd
            integer, intent(in) :: m

            type(vec3d) :: point

            real(8) :: zeta, eta, rho
            real(8) :: qx, qy, qz
            real(8) :: px, py, pz
            real(8) :: t

            zeta = ga%expo + gb%expo
            eta  = gc%expo + gd%expo
            rho = zeta*eta/(zeta+eta)

            qx = (gc%expo*gc%cx + gd%expo*gd%cx)/eta
            qy = (gc%expo*gc%cy + gd%expo*gd%cy)/eta
            qz = (gc%expo*gc%cz + gd%expo*gd%cz)/eta


            point = _pgto_internal_division_point(ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)
            px = point%x
            py = point%y
            pz = point%z

            t = (px-qx)**2 + (py-qy)**2 + (pz-qz)**2
            t = rho * t

            ret = (zeta + eta)**-0.5d0 * _pgto_kappa(ga, gb) * _pgto_kappa(gc, gd) * _pgto_fm(t, m)
        end function

        pure function _pgto_kappa(ga, gb) result(kappa)
            real(8), intent(out) :: kappa
            type(pgto), intent(in) :: ga, gb

            real(8) :: d2

            d2 = (ga%cx-gb%cx)**2 + (ga%cy-gb%cy)**2 + (ga%cz-gb%cz)**2

            kappa = 2d0**0.5d0 * PI**1.25d0 * exp(-d2*ga%expo*gb%expo/(ga%expo + gb%expo)) / (ga%expo + gb%expo)
        end function

        pure function _pgto_fm(t, m) result(fm)
            real(8), intent(out) :: fm
            real(8), intent(in) :: t
            integer, intent(in) :: m

            real(8), parameter :: DV = 1e-5
            real(8) :: sumation = 0, v = 0

            integer, parameter :: MAX_M = 16
            integer, parameter :: MAX_T = 100000
            real(8), dimension(MAX_T, MAX_M) :: buffer = 0d0
            integer :: t_idx, m_idx

            t_idx = int(t*MAX_T) + 1
            m_idx = m+1

            if ( (t_idx .gt. MAX_T) .or. (m .gt. MAX_M) ) then
                fm = _pgto_fm_do_calc(t, m)
            else
                if ( buffer(t_idx, m_idx) .eq. 0d0 ) then
                    buffer(t_idx, m_idx) = _pgto_fm_do_calc(t, m)
                end if

                fm = buffer(t_idx, m_idx)
            end if
        end function

        pure function _pgto_fm_do_calc(t, m) result(fm)
            real(8), intent(out) :: fm
            real(8), intent(in) :: t
            integer, intent(in) :: m

            real(8), parameter :: DV = 1e-5
            real(8) :: sumation = 0, v = 0

            v = 0d0
            sumation = 0d0
            do while ( v .le. 1 )
                sumation = sumation + 0.5d0 * DV * (_pgto_fm_func(t, m, v) + _pgto_fm_func(t, m, v+DV))
                v = v + DV
            end do

            fm = sumation
        end function

        pure function _pgto_fm_func(t, m, v) result(func)
            real(8), intent(out) :: func
            real(8), intent(in) :: t, v
            integer, intent(in) :: m

            func = v**(2*m) * exp(-t*v**2)
        end function

        recursive function pgto_kinetic_energy(ga, gb) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb

            type(vec3d) :: point

            real(8) :: ab, am1b, abm1, zeta, gzi, px, py, pz
            real(8) :: oap1b, oam1b

            if ( _pgto_all_n_zero(ga) .and. _pgto_all_n_zero(gb) ) then
                ret = _pgto_kinetic_energy0(ga, gb)
            else
                zeta = ga%expo + gb%expo
                gzi = ga%expo*gb%expo/zeta
                ab = 0d0
                am1b = 0d0
                abm1 = 0d0
                oap1b = 0d0
                oam1b = 0d0

                point = _pgto_internal_division_point(ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)
                px = point%x
                py = point%y
                pz = point%z

                oap1b = pgto_overlap(ga, gb)

                if ( ga%nx .gt. 0 ) then
                    ab = pgto_kinetic_energy(_pgto_clone(ga, -1, 0, 0), gb)
                    if ( ga%nx .gt. 1 ) then
                        am1b = pgto_kinetic_energy(_pgto_clone(ga, -2, 0, 0), gb)
                        oam1b = pgto_overlap(_pgto_clone(ga, -2, 0, 0), gb)
                    end if
                    if ( gb%nx .gt. 0 ) then
                        abm1 = pgto_kinetic_energy(_pgto_clone(ga, -1, 0, 0), _pgto_clone(gb, -1, 0, 0))
                    end if
                    ret = (px-ga%cx)*ab + 0.5d0*(ga%nx-1)*am1b/zeta + 0.5d0*gb%nx*abm1/zeta + &
                        2*gzi * (oap1b - 0.5d0*(ga%nx-1)*oam1b/ga%expo)
                else if ( ga%ny .gt. 0 ) then
                    ab = pgto_kinetic_energy(_pgto_clone(ga, 0, -1, 0), gb)
                    if ( ga%ny .gt. 1 ) then
                        am1b = pgto_kinetic_energy(_pgto_clone(ga, 0, -2, 0), gb)
                        oam1b = pgto_overlap(_pgto_clone(ga, 0, -2, 0), gb)
                    end if
                    if ( gb%ny .gt. 0 ) then
                        abm1 = pgto_kinetic_energy(_pgto_clone(ga, 0, -1, 0), _pgto_clone(gb, 0, -1, 0))
                    end if
                    ret = (py-ga%cy)*ab + 0.5d0*(ga%ny-1)*am1b/zeta + 0.5d0*gb%ny*abm1/zeta + &
                        2*gzi * (oap1b - 0.5d0*(ga%ny-1)*oam1b/ga%expo)
                else if ( ga%nz .gt. 0 ) then
                    ab = pgto_kinetic_energy(_pgto_clone(ga, 0, 0, -1), gb)
                    if ( ga%nz .gt. 1 ) then
                        am1b = pgto_kinetic_energy(_pgto_clone(ga, 0, 0, -2), gb)
                        oam1b = pgto_overlap(_pgto_clone(ga, 0, 0, -2), gb)
                    end if
                    if ( gb%nz .gt. 0 ) then
                        abm1 = pgto_kinetic_energy(_pgto_clone(ga, 0, 0, -1), _pgto_clone(gb, 0, 0, -1))
                    end if
                    ret = (pz-ga%cz)*ab + 0.5d0*(ga%nz-1)*am1b/zeta + 0.5d0*gb%nz*abm1/zeta + &
                        2*gzi * (oap1b - 0.5d0*(ga%nz-1)*oam1b/ga%expo)
                else
                    ret = pgto_kinetic_energy(gb, ga)
                end if
            end if
        end function

        function _pgto_kinetic_energy0(ga, gb) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb

            real(8) :: zeta, gzi, r2

            zeta = ga%expo + gb%expo
            gzi = ga%expo*gb%expo/zeta

            r2 = (ga%cx-gb%cx)**2 + (ga%cy-gb%cy)**2 + (ga%cz-gb%cz)**2
            ret = gzi * (3d0 - 2d0*gzi*r2) * _pgto_overlap0(ga, gb) ! norm multiplied in overlap0
        end function

        ! FIXME m should not be passed from caller
        recursive function pgto_nuclear_attr(ga, gb, cx, cy, cz, m) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb
            real(8), intent(in) :: cx, cy, cz
            integer, intent(in) :: m

            type(vec3d) :: point

            real(8) :: zeta
            real(8) :: abm, abmp1, am1bm, am1bmp1, abm1m, abm1mp1
            real(8) :: px, py, pz

            am1bm = 0d0
            am1bmp1 = 0d0
            abm1m = 0d0
            abm1mp1 = 0d0

            zeta = ga%expo + gb%expo
            point = _pgto_internal_division_point(ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)
            px = point%x
            py = point%y
            pz = point%z

            if ( _pgto_all_n_zero(ga) .and. _pgto_all_n_zero(gb) ) then
                ret = _pgto_nuclear_attr0(ga, gb, cx, cy, cz, m)
            else if ( ga%nx .gt. 0 ) then
                abm = pgto_nuclear_attr(_pgto_clone(ga, -1,0,0), gb, cx, cy, cz, m)
                abmp1 = pgto_nuclear_attr(_pgto_clone(ga, -1,0,0), gb, cx, cy, cz, m+1)

                if ( ga%nx .gt. 1 ) then
                    am1bm = pgto_nuclear_attr(_pgto_clone(ga, -2,0,0), gb, cx, cy, cz, m)
                    am1bmp1 = pgto_nuclear_attr(_pgto_clone(ga, -2,0,0), gb, cx, cy, cz, m+1)
                end if

                if ( gb%nx .gt. 0 ) then
                    abm1m = pgto_nuclear_attr(_pgto_clone(ga, -1,0,0), _pgto_clone(gb, -1,0,0), cx, cy, cz, m)
                    abm1mp1 = pgto_nuclear_attr(_pgto_clone(ga, -1,0,0), _pgto_clone(gb, -1,0,0), cx, cy, cz, m+1)
                end if

                ret = (px-ga%cx)*abm - (px-cx)*abmp1 + 0.5d0*(ga%nx-1)*(am1bm - am1bmp1)/zeta &
                    + 0.5d0*gb%nx*(abm1m - abm1mp1)/zeta
            else if ( ga%ny .gt. 0 ) then
                abm = pgto_nuclear_attr(_pgto_clone(ga, 0,-1,0), gb, cx, cy, cz, m)
                abmp1 = pgto_nuclear_attr(_pgto_clone(ga, 0,-1,0), gb, cx, cy, cz, m+1)

                if ( ga%ny .gt. 1 ) then
                    am1bm = pgto_nuclear_attr(_pgto_clone(ga, 0,-2,0), gb, cx, cy, cz, m)
                    am1bmp1 = pgto_nuclear_attr(_pgto_clone(ga, 0,-2,0), gb, cx, cy, cz, m+1)
                end if

                if ( gb%ny .gt. 0 ) then
                    abm1m = pgto_nuclear_attr(_pgto_clone(ga, 0,-1,0), _pgto_clone(gb, 0,-1,0), cx, cy, cz, m)
                    abm1mp1 = pgto_nuclear_attr(_pgto_clone(ga, 0,-1,0), _pgto_clone(gb, 0,-1,0), cx, cy, cz, m+1)
                end if

                ret = (py-ga%cy)*abm - (py-cy)*abmp1 + 0.5d0*(ga%ny-1)*(am1bm - am1bmp1)/zeta &
                    + 0.5d0*gb%ny*(abm1m - abm1mp1)/zeta
            else if ( ga%nz .gt. 0 ) then
                abm = pgto_nuclear_attr(_pgto_clone(ga, 0,0,-1), gb, cx, cy, cz, m)
                abmp1 = pgto_nuclear_attr(_pgto_clone(ga, 0,0,-1), gb, cx, cy, cz, m+1)

                if ( ga%nz .gt. 1 ) then
                    am1bm = pgto_nuclear_attr(_pgto_clone(ga, 0,0,-2), gb, cx, cy, cz, m)
                    am1bmp1 = pgto_nuclear_attr(_pgto_clone(ga, 0,0,-2), gb, cx, cy, cz, m+1)
                end if

                if ( gb%nz .gt. 0 ) then
                    abm1m = pgto_nuclear_attr(_pgto_clone(ga, 0,0,-1), _pgto_clone(gb, 0,0,-1), cx, cy, cz, m)
                    abm1mp1 = pgto_nuclear_attr(_pgto_clone(ga, 0,0,-1), _pgto_clone(gb, 0,0,-1), cx, cy, cz, m+1)
                end if

                ret = (pz-ga%cz)*abm - (pz-cz)*abmp1 + 0.5d0*(ga%nz-1)*(am1bm - am1bmp1)/zeta &
                    + 0.5d0*gb%nz*(abm1m - abm1mp1)/zeta
            else
                ret = pgto_nuclear_attr(gb, ga, cx, cy, cz, m)
            end if
        end function

        function _pgto_nuclear_attr0(ga, gb, cx, cy, cz, m) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb
            real(8), intent(in) :: cx, cy, cz
            integer, intent(in) :: m

            type(vec3d) :: point

            real(8) :: zeta, u, px, py, pz

            point = _pgto_internal_division_point(ga%expo, ga%cx, ga%cy, ga%cz, gb%expo, gb%cx, gb%cy, gb%cz)
            px = point%x
            py = point%y
            pz = point%z

            zeta = ga%expo + gb%expo
            u = zeta * ((px-cx)**2 + (py-cy)**2 + (pz-cz)**2)
            
            ret = 2d0*(zeta/PI)**0.5d0 * _pgto_overlap0(ga, gb) * _pgto_fm(u, m) ! norm multiplied in overlap0
        end function
end module
