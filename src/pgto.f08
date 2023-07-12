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

end module
