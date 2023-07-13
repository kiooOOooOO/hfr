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

end module
