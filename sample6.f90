module fenep
    implicit none
contains
    subroutine tdma(a, b, c, d, x, size)
        ! https://mathlang.hatenablog.com/entry/2018/10/14/232741
        implicit none
        integer, intent(in) :: size
        real(8), intent(in) :: a(size), b(size), c(size), d(size)
        real(8), intent(out) :: x(size)
        real(8) :: P(size), Q(size)
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
    integer, parameter :: size = 4 ! 3x3行列を考える
    real(8) :: a(size), b(size), c(size), d(size), x(size)
    integer :: i
    
    ! 行列Aを設定する
    a = [2.0d0, 2.0d0, 2.0d0, 2.0d0]
    b = [1.0d0, 1.0d0, 1.0d0, 0.0d0]
    c = [0.0d0, 1.0d0, 1.0d0, 1.0d0]
    d = [4.0d0, 8.0d0, 12.0d0, 11.0d0]

    ! tdmaサブルーチンで解く
    call tdma(a, b, c, d, x, size)
    
    ! 結果を出力する
    do i = 1, size
        print *, 'x(', i, ') = ', x(i)  ! 1, 2, 3, 4
    enddo

    
end program main


! gfortran sample6.f90 && ./a.out
