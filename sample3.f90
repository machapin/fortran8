module fenep
    implicit none
    real(8), parameter :: Lp = 55
    real(8), parameter :: tau = 1.0
    real(8), parameter :: dt = 1.0d-3

contains
    function f(C) result(y)  ! FENE-Pモデルの核
        real(8), intent(in) :: C(6) 
        real(8) y  ! y=f(C)
        y = (Lp**2 - 3)/(Lp**2 - C(1)-C(4)-C(6))

    end function f


    subroutine Lyapunov(Cx, U, V, W, C)
        real(8), intent(in) :: Cx(6), U(3, 3, 3), V(3, 3, 3), W(3, 3, 3)
        real(8), intent(out) :: C(6)
        integer i, j
        real(8) A0(3, 3), b(6)
        real(8) x(6), x1(6), x2(6), l(6), l1(6), l2(6), Jac(6, 6), dx
        integer itr, itrmax
        real(8) r(6), eps
        dx = 1.0d-5  ! 3, 4, 5あたりがいい。
        itrmax = 100
        eps = 1.0d-10

        ! 定数部分を計算
        A0(1, 1) = U(1, 1, 1)+V(1, 1, 1)+W(1, 1, 1)
        ! A0(1, 1) = -dt*(U(i+1, j, k)-U(i, j, k))/dX  ! 要素1つずつ計算していくしかない。スタッガード格子に注意
        ! A0(1, :) = (/1.0, 0.0, 0.0/)
        ! A0(2, :) = (/0.0, 1.0, 0.0/)
        ! A0(3, :) = (/0.0, 0.0, 1.0/)
        call random_number(A0)
        A0 = transpose(A0)
        b(:) = Cx(:)
        b(1) = b(1) + dt/tau
        b(4) = b(4) + dt/tau
        b(6) = b(6) + dt/tau
        
        ! Newton法
        x(:) = Cx(:)  ! 探索開始点、なんとなくCxが近そう。
        do itr = 1, itrmax
            ! ヤコビアンの計算
            do j = 1, 6
                x1(:) = x(:)
                x2(:) = x(:)
                x1(j) = x(j) + dx
                x2(j) = x(j) - dx
                call Lyapunov_func(A0, b, x, l)
                call Lyapunov_func(A0, b, x, l)
                do i = 1, 6
                    Jac(i, j) = (l1(i) - l2(i))/(2*dx)
                enddo
            enddo

            call Lyapunov_func(A0, b, x, l)  ! f(x_i)を計算

            call LU_solve(Jac, l, r)
            x(:) = x(:) - r(:)

            ! write(*, '(6e12.4)') r(:)
            write(*, '(6e12.4)') sum(r**2), sum(x**2), sqrt(sum(r**2)/sum(x**2))
            ! write(*, *) ''

            if (sqrt(sum(r**2)/sum(x**2)) < eps) exit
        enddo

        C(:) = x(:)
        

    end subroutine Lyapunov


    subroutine Lyapunov_func(A0, b, x, l)  ! l(x)=Mx-bを計算
        real(8), intent(in) :: A0(3, 3), b(6), x(6)
        real(8), intent(out) :: l(6)
        integer i
        real(8) A(3, 3), M(6, 6)

        A(:, :) = A0(:, :)
        A(1, 1) = A(1, 1) + (1 + f(x)*dt/tau)/2
        A(2, 2) = A(2, 2) + (1 + f(x)*dt/tau)/2
        A(3, 3) = A(3, 3) + (1 + f(x)*dt/tau)/2

        M(:, :) = 0.0d0
        M(:, :) = 0.0d0

        M(1, 1) = A(1, 1) + A(1, 1)
        M(1, 2) = A(1, 2)
        M(1, 3) = A(1, 3)
 
        M(2, 1) = A(2, 1)
        M(2, 2) = A(2, 2) + A(1, 1)
        M(2, 3) = A(2, 3)
        M(2, 4) = A(1, 2)
        
        M(3, 1) = A(3, 1)
        M(3, 2) = A(3, 2)
        M(3, 3) = A(3, 3) + A(1, 1)
        M(3, 5) = A(1, 2)
        M(3, 6) = A(1, 3)

        M(4, 2) = A(2, 1)
        M(4, 4) = A(2, 2) + A(2, 2)
        M(4, 5) = A(2, 3)

        M(5, 3) = A(2, 1)
        M(5, 4) = A(3, 2)
        M(5, 5) = A(3, 3) + A(2, 2)
        M(5, 6) = A(2, 3)

        M(6, 3) = A(3, 1)
        M(6, 5) = A(3, 2)
        M(6, 6) = A(3, 3) + A(3, 3)

        do i = 1, 6
            l(i) = sum(M(i, :)*x(:)) - b(i)
        enddo
    end subroutine Lyapunov_func


    subroutine LU_solve(a, b, x)  ! 外積形式ガウス法
        real(8), intent(inout) :: a(6, 6)
        real(8), intent(in) :: b(6)
        real(8), intent(out) :: x(6)
        integer i, j, k
        real(8) tmp, s, y(6)

        do k = 1, 6
            tmp = 1.0d0 / a(k, k)
            a(k+1:6, k) = a(k+1:6, k) * tmp
            do j = k + 1, 6
                tmp = a(k, j)
                a(k+1:6, j) = a(k+1:6, j) - a(k+1:6, k) * tmp
            enddo
        enddo

        ! Ax=LUx=b
        do i = 1, 6  !「Ly=b」前進代入
            s = sum((/(a(i, j)*y(j), j=1, i)/))
            y(i) = b(i) - s
        enddo
        do i = 6, 1, -1  !「Ux=y」後退代入
            s = sum((/(a(i, j)*x(j), j=i+1, 6)/))
            x(i) = (y(i) - s)/a(i, i)
        enddo
    end subroutine LU_solve

end module fenep



program main
    use fenep
    implicit none
    real(8) Cx(6), U(3, 3, 3), V(3, 3, 3), W(3, 3, 3), C(6)
    real(8) re0, re1, re2
    integer i, count, NX
    call random_number(U)
    call random_number(V)
    call random_number(W)
    
    NX = 1000000
    count = 0
    do i = 1, NX
        call random_number(Cx)
        C(:) = 0.1d0 * Cx(:)
        C(1) = 0.1d0 + Cx(1)
        C(4) = 0.1d0 + Cx(2)
        C(6) = 0.1d0 + Cx(3)
        call Cardano(C(:), re0, re1, re2)
        if (re0 < 0.0d0 .or. re1 < 0.0d0 .or. re2 < 0.0d0) count = count + 1
    enddo
    write(*, '(F6.4)') count*1.0d0/NX



    ! Cx(:) = (/1.0, 2.0, -1.0, -2.0, 2.0, 1.0/)
    ! call random_number(Cx)

    ! call Lyapunov(Cx, U, V, W, C)
    
end program main

