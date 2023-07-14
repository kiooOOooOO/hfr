module sto_ng
    use pgto
    implicit none

    type sto_ng
        real(8), allocatable, dimension(:) :: coefs
        type(pgto), allocatable, dimension(:) :: pgtos
        integer :: n
        real(8) :: norm
    end type

    contains

        function stong_new(n, pgtos, coefs)
            type(sto_ng), intent(out) :: stong_new
            integer, intent(in) :: n
            type(pgto), intent(in), dimension(:) :: pgtos
            real(8), intent(in), dimension(:) :: coefs

            stong_new%n = n
            allocate(stong_new%coefs(n))
            allocate(stong_new%pgtos(n))

            stong_new%coefs = coefs(1:n)
            stong_new%pgtos = pgtos(1:n)

            stong_new%norm = _stong_normalization_factor(n, pgtos, coefs)
        end function

        function _stong_normalization_factor(n, pgtos,coefs) result(ret)
            real(8), intent(out) :: ret
            integer, intent(in) :: n
            type(pgto), intent(in), dimension(:) :: pgtos
            real(8), intent(in), dimension(:) :: coefs

            real(8) :: r
            integer :: i1, i2

            r = 0d0
            do i1=1,n
            do i2=1,n
                r = r + coefs(i1)*coefs(i2)*pgto_overlap(pgtos(i1),pgtos(i2))
            end do
            end do

            ret = r**-0.5d0
        end function

        function stong_overlap(sa, sb)
            real(8), intent(out) :: stong_overlap
            type(sto_ng), intent(in) :: sa, sb

            real(8) :: r
            integer :: a1, b1

            r = 0d0
            do a1=1,sa%n
            do b1=1,sb%n
                r = r + sa%coefs(a1) * sb%coefs(b1) * pgto_overlap(sa%pgtos(a1), sb%pgtos(b1))
            end do
            end do

            stong_overlap = sa%norm * sb%norm * r
        end function

        pure function stong_eri(sa, sb, sc, sd)
            real(8), intent(out) :: stong_eri
            type(sto_ng), intent(in) :: sa, sb, sc, sd

            real(8) :: r
            integer :: a1, b1, c1, d1

            r =0d0
            do a1=1,sa%n
            do b1=1,sb%n
            do c1=1,sc%n
            do d1=1,sd%n
                r = r + sa%coefs(a1) * sb%coefs(b1) * sc%coefs(c1) * sd%coefs(d1) &
                    * pgto_eri(sa%pgtos(a1), sb%pgtos(b1), sc%pgtos(c1), sd%pgtos(d1))
            end do
            end do
            end do
            end do

            stong_eri = sa%norm * sb%norm * sc%norm * sd%norm * r
        end function

        function stong_kinetic_energy(sa, sb)
            real(8), intent(out) :: stong_kinetic_energy
            type(sto_ng), intent(in) :: sa, sb

            real(8) :: r
            integer :: a1, b1

            r = 0d0
            do a1=1,sa%n
            do b1=1,sb%n
                r = r + sa%coefs(a1) * sb%coefs(b1) * pgto_kinetic_energy(sa%pgtos(a1), sb%pgtos(b1))
            end do
            end do

            stong_kinetic_energy = sa%norm * sb%norm * r
        end function

        function stong_nuclear_attr(sa, sb, cx, cy, cz)
            real(8), intent(out) :: stong_nuclear_attr
            type(sto_ng), intent(in) :: sa, sb
            real(8), intent(in) :: cx, cy, cz

            real(8) :: r
            integer :: a1, b1

            r = 0d0
            do a1=1,sa%n
            do b1=1,sb%n
                r = r + sa%coefs(a1) * sb%coefs(b1) * pgto_nuclear_attr(sa%pgtos(a1), sb%pgtos(b1), cx, cy, cz, 0)
            end do
            end do

            stong_nuclear_attr = sa%norm * sb%norm * r
        end function

end module
