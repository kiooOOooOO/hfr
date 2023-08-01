module pgto2
    use pgto
    implicit none

    contains

        ! ∬ dr1 dr2 ga(1)gb(1)gc(2)gd(2)/r12
        pure function pgto2_eri(ga, gb, gc, gd) result(ret)
            real(8), intent(out) :: ret
            type(pgto), intent(in) :: ga, gb, gc, gd

            ret = ga%norm * gb%norm * gc%norm * gd%norm * _pgto2_eri_internal(ga, gb, gc, gd, 0)
        end function

        recursive pure function _pgto2_eri_internal(pga, pgb, pgc, pgd, m) result(ret)
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

!                write (*,*) "initial"
!                    write (*,'("A(",I0,",",I0,",",I0,")")') ga%nx, ga%ny, ga%nz
!                    write (*,'("B(",I0,",",I0,",",I0,")")') gb%nx, gb%ny, gb%nz
!                    write (*,'("C(",I0,",",I0,",",I0,")")') gc%nx, gc%ny, gc%nz
!                    write (*,'("D(",I0,",",I0,",",I0,")")') gd%nx, gd%ny, gd%nz
            if ( za .and. zb .and. zc .and. zd ) then
                ret = _pgto_eri_internal0(ga, gb, gc, gd, m)
            else if ( za .and. ( .not. zb ) ) then
                ret = _pgto2_eri_internal(gb, ga, gc, gd, m)
            else if ( (za .and. zb) .and. ((.not. zc) .or. (.not. zd)) ) then
                ret = _pgto2_eri_internal(gc, gd, ga, gb, m)
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
                    am = _pgto2_eri_internal(ga, gb, gc, gd, m)
                    amp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)

                    if ( ga%nx .gt. 0 ) then
                        ga%nx = ga%nx - 1
                        am1m = _pgto2_eri_internal(ga, gb, gc, gd, m)
                        am1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        ga%nx = ga%nx + 1
                    else
                        am1m = 0d0
                        am1mp1 = 0d0
                    end if

                    if ( gb%nx .gt. 0 ) then
                        gb%nx = gb%nx - 1
                        bm1m = _pgto2_eri_internal(ga, gb, gc, gd, m)
                        bm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        gb%nx = gb%nx + 1
                    else
                        bm1m = 0d0
                        bm1mp1 = 0d0
                    end if

                    if ( gc%nx .gt. 0 ) then
                        gc%nx = gc%nx - 1
                        cm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        gc%nx = gc%nx + 1
                    else
                        cm1mp1 = 0d0
                    end if

                    if ( gd%nx .gt. 0 ) then
                        gd%nx = gd%nx - 1
                        dm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
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
                    am = _pgto2_eri_internal(ga, gb, gc, gd, m)
                    amp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)

                    if ( ga%ny .gt. 0 ) then
                        ga%ny = ga%ny - 1
                        am1m = _pgto2_eri_internal(ga, gb, gc, gd, m)
                        am1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        ga%ny = ga%ny + 1
                    else
                        am1m = 0d0
                        am1mp1 = 0d0
                    end if

                    if ( gb%ny .gt. 0 ) then
                        gb%ny = gb%ny - 1
                        bm1m = _pgto2_eri_internal(ga, gb, gc, gd, m)
                        bm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        gb%ny = gb%ny + 1
                    else
                        bm1m = 0d0
                        bm1mp1 = 0d0
                    end if

                    if ( gc%ny .gt. 0 ) then
                        gc%ny = gc%ny - 1
                        cm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        gc%ny = gc%ny + 1
                    else
                        cm1mp1 = 0d0
                    end if

                    if ( gd%ny .gt. 0 ) then
                        gd%ny = gd%ny - 1
                        dm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
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
                    am = _pgto2_eri_internal(ga, gb, gc, gd, m)
                    amp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)

                    if ( ga%nz .gt. 0 ) then
                        ga%nz = ga%nz - 1
                        am1m = _pgto2_eri_internal(ga, gb, gc, gd, m)
                        am1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        ga%nz = ga%nz + 1
                    else
                        am1m = 0d0
                        am1mp1 = 0d0
                    end if

                    if ( gb%nz .gt. 0 ) then
                        gb%nz = gb%nz - 1
                        bm1m = _pgto2_eri_internal(ga, gb, gc, gd, m)
                        bm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        gb%nz = gb%nz + 1
                    else
                        bm1m = 0d0
                        bm1mp1 = 0d0
                    end if

                    if ( gc%nz .gt. 0 ) then
                        gc%nz = gc%nz - 1
                        cm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
                        gc%nz = gc%nz + 1
                    else
                        cm1mp1 = 0d0
                    end if

                    if ( gd%nz .gt. 0 ) then
                        gd%nz = gd%nz - 1
                        dm1mp1 = _pgto2_eri_internal(ga, gb, gc, gd, m+1)
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

end module
