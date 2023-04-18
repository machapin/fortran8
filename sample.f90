module fenep
    implicit none
    real(8), parameter :: Lp = 55
    real(8), parameter :: tau = 1.0
    real(8), parameter :: dt = 1.0d-3

contains
    subroutine Eigenvalues(a, re1, re2)
        real(8), intent(in) :: a(3)
        real(8), intent(out) :: re1, Re2
        real(8) D

        D = (a(1)-a(3))**2 + 4.0d0 * a(2)**2
        re1 = (a(1) + a(3) + sqrt(D))*0.5d0
        re2 = (a(1) + a(3) - sqrt(D))*0.5d0

    end subroutine Eigenvalues
end module fenep

program main
    use fenep
    implicit none
    real(8) C(3)
    real(8) re1, re2

    C(:) = (/8.0, 2.0, 5.0/)
    call Eigenvalues(C, re1, re2)
    write(*, *) re1, re2
    
end program main

