module fenep
    implicit none
contains
    subroutine tdma(a, b, c, d, x, size)
        implicit none
        integer, intent(in) :: size
        complex(8), intent(in) :: a(size), b(size), c(size), d(size)
        complex(8), intent(out) :: x(size)
        complex(8) :: P(size), Q(size)
        integer :: i
      
        P(1) = -b(1)/a(1)
        Q(1) = d(1)/a(1)
      
        do i = 2, size
          P(i) = -b(i)/(a(i)+c(i)*P(i-1))
          Q(i) = (d(i)-c(i)*Q(i-1))/(a(i)+c(i)*P(i-1))
        enddo

        x(size) = Q(size)
      
        do i = size-1, 1, -1
          x(i) = P(i)*x(i+1) + Q(i)
        enddo
      
      end subroutine tdma
end module fenep


program main
    use fenep
    implicit none
    integer, parameter :: n = 4 ! 連立方程式のサイズ
    integer :: i
    complex(8) :: a(n), b(n), c(n), d(n), x(n) ! 複素数型の配列
    
    ! 係数の設定
    ! a = [cmplx(1.0, 0.0), cmplx(1.0, 0.0), cmplx(-1.0, 0.0), cmplx(1.0, 0.0)]
    ! b = [cmplx(1.0, 0.0), cmplx(1.0, 0.0), cmplx(1.0, 0.0), cmplx(0.0, 0.0)]
    ! c = [cmplx(0.0, 0.0), cmplx(1.0, 0.0), cmplx(1.0, 0.0), cmplx(1.0, 0.0)]
    ! d = [cmplx(0.0, 0.0), cmplx(4.0, 2.0), cmplx(-4.0, -2.0), cmplx(2.0, 1.0)]
    a = cmplx(2.0, 1.0)
    b = cmplx(1.0, 2.0)
    c = cmplx(1.0, -2.0)
    d = cmplx(1.0, 1.0)
    

    ! 三重対角行列の連立方程式を解く
    call tdma(a, b, c, d, x, n)

    ! 結果の出力
    do i = 1, n
        write(*,*) "x(",i,") = ", x(i)
    enddo
    
end program main

! gfortran sample6.f90 && ./a.out
