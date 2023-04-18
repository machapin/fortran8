module fenep
    implicit none
    integer, parameter :: N = 4
    real(8), parameter :: Lp = 55
    real(8), parameter :: tau = 1.0
    real(8), parameter :: dt = 1.0d-3

contains
    subroutine Newton(A, b, x)  ! Aが変化しない場合
        ! 参考元https://slpr.sakura.ne.jp/qp/rootfinding/#newtonndprogram
        ! call dgetrf(N, N, Jac, N, ipiv, info)  ! LAPACKをインストールする必要あり
        real(8), intent(in) :: A(N, N), b(N)
        real(8), intent(out) :: x(N)
        real(8) Jac(N, N)
        real(8) w(N), f(N)
        integer itr, itrmax
        real(8) eps
        itrmax = 10
        eps = 1.0d-5

        x(:) = 0.0d0  ! 探索開始点
        do itr = 1, itrmax
            call Jacobian(A, b, x, Jac)  ! ヤコビアンの計算

            f(:) = matmul(A, x)-b(:)  ! 定数のときのみ
            ! Jac(1, :) = (/1.0, -2.0, 3.0/)
            ! Jac(2, :) = (/2.0, -1.0, 1.0/)
            ! Jac(3, :) = (/1.0, 3.0, -5.0/)
            ! f(:) =(/5.0, 6.0, 2.0/)
            call LU_solve(Jac, f, w)
            x(:) = x(:) - w(:)

            write(*, '(6e12.4)') w(:)
            write(*, '(6e12.4)') sum(w**2)/sum(x**2)
            write(*, *) ''

            if (sqrt(sum(w**2)/sum(x**2)) < eps) exit
        enddo

    end subroutine Newton

    subroutine Jacobian(A, b, x, Jac)
        real(8), intent(in) :: A(N, N), b(N), x(N)
        real(8), intent(out) :: Jac(N, N)
        real(8), parameter :: dx = 1.0d-5
        integer i, j
        real(8) x1(N), x2(N)
        Jac(:, :) = 0.0d0

        do j = 1, N
            x1(:) = x(:)
            x2(:) = x(:)
            x1(j) = x(j) + dx
            x2(j) = x(j) - dx
            do i = 1, N
                Jac(i, j) = (sum(A(i, :)*x1(:))-sum(A(i, :)*x2(:)))/(2*dx)  ! いやこれ、傾き変化せんかったらAそのものやん。
            enddo
        enddo
        ! write(*, '(6e12.4)') Jac(:, :)
        x1(1) = b(1)  ! 使ってません

    end subroutine Jacobian


    subroutine LU_solve(a, b, x)  ! 外積形式ガウス法
        real(8), intent(inout) :: a(N, N)
        real(8), intent(in) :: b(N)
        real(8), intent(out) :: x(N)
        integer i, j, k
        real(8) s, y(N)
        y(:) = 0.0d0
        x(:) = 0.0d0

        do k = 1, N  ! Lの対角成分は1、LとUを一つの行列にまとめる
            a(k+1:N, k) = a(k+1:N, k) / a(k, k)
            do j = k + 1, N
                a(k+1:N, j) = a(k+1:N, j) - a(k+1:N, k) * a(k, j)
            enddo
        enddo

        ! Ax=LUx=b
        do i = 1, N  !「Ly=b」前進代入
            s = sum((/(a(i, j)*y(j), j=1, i)/))
            y(i) = b(i) - s
        enddo
        do i = N, 1, -1  !「Ux=y」後退代入
            s = sum((/(a(i, j)*x(j), j=i+1, N)/))
            x(i) = (y(i) - s)/a(i, i)
        enddo
    end subroutine LU_solve

    subroutine gauss(a, b, x)
        integer, parameter :: n = 3
        real(8), intent(in) :: a(n, n), b(n)
        real(8), intent(out) :: x(n)
        integer i, k, maxidx
        real(8) maxabs, temp(n+1), Ah(n, n+1)
        Ah(:, :n) = a(:, :)
        Ah(:, n+1) = b(:)

        do k = 1, n
            ! ピボット選択して並び替え
            maxidx = maxloc(abs(Ah(k:n, k)), 1) + k-1
            maxabs = abs(Ah(maxidx, k))

            if (maxabs == 0.0d0) then  ! 解けない場合 x(:) = 0.0d0
                Ah(:, n+1) = 0.0d0
                ! stop 'gauss cannot solve'
                exit
            endif

            if (k /= maxidx) then
                temp(k:n+1) = Ah(k, k:n+1)
                Ah(k, k:n+1) = Ah(maxidx, k:n+1)
                Ah(maxidx, k:n+1) = temp(k:n+1)
            endif

            Ah(k, k+1:n+1) = Ah(k, k+1:n+1) / Ah(k, k)
            Ah(k, k) = 1.0d0
            do i = 1, n
                if (i /= k) then
                    Ah(i, k+1:n+1) = Ah(i, k+1:n+1) - Ah(i, k) * Ah(k, k+1:n+1)
                    Ah(i, k) = 0.0d0
                endif
            enddo
        enddo
        x(:) = Ah(:, n+1)
    end subroutine gauss
end module fenep


program main
    use fenep
    implicit none
    ! real(8) A(N, N), b(N), x(N)
    real(8) A(3, 3), b(3), x(3)

    ! A(1, :) = (/1.0, 2.0/)
    ! A(2, :) = (/4.0, 5.0/)
    ! b(:)= (/3.0, 6.0/)  ! -1, 2

    ! A(1, :) = (/1.0, -2.0, 3.0/)
    ! A(2, :) = (/2.0, -1.0, 1.0/)
    ! A(3, :) = (/1.0, 3.0, -5.0/)
    ! b(:) = (/5.0, 6.0, 2.0/)  ! 6, 17, 11

    A(1, :) = (/1.0, 1.0, 1.0/)  !  エラー発生
    A(2, :) = (/1.0, 1.0, 1.0/)
    A(3, :) = (/1.0, -1.0, 1.0/)
    b(:) = (/1.0, 4.0, -4.0/)  ! 0, 2.5, -1.5

    ! A(1, :) = (/1.0, 1.0, 1.0, 1.0/)
    ! A(2, :) = (/1.0, 1.0, 1.0, -1.0/)
    ! A(3, :) = (/1.0, 1.0, -1.0, 1.0/)
    ! A(4, :) = (/1.0, -1.0, 1.0, 1.0/)
    ! b(:) = (/0.0, 4.0, -4.0, 2.0/)  ! 1, -1, 2, -2

    ! A(1, :) = (/1.0, 1.0, 1.0, 1.0, -1.0, 1.0/)
    ! A(2, :) = (/1.0, 1.0, 1.0, -1.0, 1.0, 1.0/)
    ! A(3, :) = (/1.0, 1.0, -1.0, 1.0, 1.0, 1.0/)
    ! A(4, :) = (/1.0, -1.0, 1.0, 1.0, 1.0, 1.0/)
    ! A(5, :) = (/-1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
    ! A(6, :) = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
    ! b(:) = (/0.0, 1.0, 2.0, 3.0, 4.0, 5.0/)  ! 0.5, 1, 1.5, 2, 2.5, -2.5
    ! call LU_solve(A, b, x)
    call gauss(A, b, x)
    write(*, *) x(:)
    
end program main

