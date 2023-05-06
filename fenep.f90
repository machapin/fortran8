module smac
    implicit none
    ! gfortran fenep.f90 -I$HOME/local/include -L$HOME/local/lib -lfftw3 && ./a.out
    ! ステップ数
    integer, parameter :: Nstep = 1
    integer, parameter :: Gstep = 1000  ! データを取得する間隔
    integer, parameter :: Estep = 1000  ! エネルギースペクトルを取得する間隔
    integer, parameter :: Dstep = 10  ! デバッグする間隔
    character(*), parameter :: dir = './debug/'
    integer, parameter :: input_step = 0  ! 0以外で初期条件をファイルから読み込む
    integer, parameter :: output_step = 100000  ! 配列を保存する間隔
    ! 手法
    integer, parameter :: method = 2  ! 0:陽解法、1:FFT、2:IBM
    real(8), parameter :: PI = acos(-1.0d0)
    ! パラメータ
    integer, parameter :: NX = 32, NY = NX, NZ = NX
    real(8), parameter :: dX = 2*PI/NX, dY = 2*PI/NY, dZ = 2*PI/NZ
    real(8), parameter :: dt = 0.01d0
    ! 無次元パラメータ
    real(8), parameter :: Re_s = 1.0d0
    real(8), parameter :: beta = 0.9d0
    real(8), parameter :: Re = Re_s*beta
    real(8), parameter :: Wi = 1.0d0
    real(8), parameter :: Lp = 55.0d0
    ! 有次元パラメータ
    real(8), parameter :: L_C = 1.0d0  ! 長さが100なら100/2*PI
    real(8), parameter :: U_C = 1.0d0  ! 本来は乱流テイラーグリーン渦の平均流の速さ
    real(8), parameter :: f0 = 1.0d0  ! 無次元化したときに1となるように
    real(8), parameter :: dX_C = dX*L_C, dY_C = dY*L_C, dZ_C = dZ*L_C
    real(8), parameter :: dt_C = dt*L_C/U_C
    real(8), parameter :: nu = L_C*U_C/Re_s
    ! グローバル変数
    integer eigen_method  ! 固有値の求め方
    real(8) counter(0:3), newton_itr
    integer poisson_itr
    ! LAPACK用変数
    ! integer, parameter :: N = 3
    ! integer :: info
    ! real(8) :: A(N, N), w0(N), work(3*N-1)
    ! integer :: lwork = 3*N-1


contains
    subroutine init(U, V, W, P, Phi, C, Fx, Fy, Fz)  ! はじめに実行
        real(8), intent(out) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: P(0:NX+1, 0:NY+1, 0:NZ+1), Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
        integer i, j, k
        real(8), parameter :: large_K = 2.0d0

        if (method == 0) then
            write(*, '(a, F8.3, F8.3, F8.3, a8)') 'Von Neumann:', beta/Re*dt/dX**2, beta/Re*dt/dY**2, beta/Re*dt/dZ**2, '< 0.167'
        endif
        write(*, '(a, F8.3, F8.3, F8.3, a8)') 'CFL :', 1.0d0*dt/dX, 1.0d0*dt/dY, 1.0d0*dt/dZ,'< 0.167'
        ! write(*, '(a5, F8.3)') 'tau =', tau
        ! write(*, '(a5, F8.3)') 'nu  =', beta/Re*U_C*Xmax  ! beta/Re*U_C*Xmax：実際のnu、beta/Re：計算上のnu
        ! write(*, '(a, F9.4)') 'U_C =', U_C
        ! write(*, '(a, F9.4)') 'Re =', Re
        ! write(*, '(a, F9.4)') 'nu_s =', nu_s
        ! write(*, '(a, F12.8)') 'dt =', dt

        ! 初期条件（Taylor-Green）
        ! U(:, :, :) = 0.0d0
        ! V(:, :, :) = 0.0d0
        ! W(:, :, :) = 0.0d0
        call random_number(U)
        call random_number(V)
        call random_number(W)
        U(:, :, :) = 0.001d0 * (U(:, :, :)-0.5d0)
        V(:, :, :) = 0.001d0 * (V(:, :, :)-0.5d0)
        W(:, :, :) = 0.001d0 * (W(:, :, :)-0.5d0)
        P(:, :, :) = 0.0d0
        Phi(:, :, :) = 0.0d0
        Fx(:, :, :) = 0.0d0
        Fy(:, :, :) = 0.0d0
        Fz(:, :, :) = 0.0d0
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    ! U(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)  ! テイラー渦
                    ! V(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)

                    ! Fx(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY) * cos((k-0.5d0)*dZ)  ! 荒木さん
                    ! Fy(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY) * cos((k-0.5d0)*dZ)
                    Fx(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)  ! 増田さん, 安房井さん
                    Fy(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)
                    ! Fx(i, j, k) = -sin(i*dX) * cos((k-0.5d0)*dZ)  ! x, z 増田さん, 安房井さん
                    ! Fz(i, j, k) = cos((i-0.5d0)*dX) * sin(k*dZ)
                    ! Fx(i, j, k) = -sin(large_K*(j-0.5d0)*dY)  ! 小井手さん
                    ! Fy(i, j, k) = sin(large_K*(i-0.5d0)*dX)
                enddo
            enddo
        enddo
        Fx(:, :, :) = f0 * Fx(:, :, :)
        Fy(:, :, :) = f0 * Fy(:, :, :)
        Fz(:, :, :) = f0 * Fz(:, :, :)

        call random_number(C)
        C(:, :, :, :) = 0.001d0 * (C(:, :, :, :)-0.5d0)
        C(1, :, :, :) = C(1, :, :, :) + 1.0d0
        C(4, :, :, :) = C(4, :, :, :) + 1.0d0
        C(6, :, :, :) = C(6, :, :, :) + 1.0d0

        ! C(:, :, :, :) = 0.0d0
        ! C(1, :, :, :) = 1.2d0
        ! C(4, :, :, :) = 1.1d0
        ! C(6, :, :, :) = 1.0d0

        call PBM(U)
        call PBM(V)
        call PBM(W)
        call PBM(P)
        call C_PBM(C)

    end subroutine init


    subroutine CpxCnx(C, Cpx, Cnx)
        real(8), intent(in) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Cpx(6, 0:NX, NY, NZ), Cnx(6, 0:NX, NY, NZ)
        integer i, j, k, l, index
        real(8) Cptemp(6, 3), Cntemp(6, 3), Eigen(3, 6), mintemp
        real(8) temp
        temp = 0.0d0
        counter(:) = 0.0d0

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Cptemp(:, 1) = C(:, i, j, k)*3/2 - C(:, i+1, j, k)/2
                    Cntemp(:, 1) = C(:, i, j, k)/2 + C(:, i+1, j, k)/2
                    Cptemp(:, 2) = C(:, i-1, j, k)/2 + C(:, i, j, k)/2
                    Cntemp(:, 2) = -C(:, i-1, j, k)/2 + C(:, i, j, k)*3/2
                    Cptemp(:, 3) = C(:, i-1, j, k)/4 + C(:, i, j, k) - C(:, i+1, j, k)/4
                    Cntemp(:, 3) = -C(:, i-1, j, k)/4 + C(:, i, j, k) + C(:, i+1, j, k)/4

                    Eigen(:, :) = 0.0d0
                    do l = 1, 3
                        if (eigen_method == 0) call Cardano(Cptemp(:, l), Eigen(l, 1), Eigen(l, 2), Eigen(l, 3))  ! CardanoでCpnの95%ぐらい時間かかってる。
                        if (eigen_method == 0) call Cardano(Cntemp(:, l), Eigen(l, 4), Eigen(l, 5), Eigen(l, 6))
                        if (eigen_method == 1) call Sylvester(Cptemp(:, l), Eigen(l, 1))
                        if (eigen_method == 1) call Sylvester(Cntemp(:, l), Eigen(l, 4))
                    enddo

                    index = 0
                    mintemp = 0.0d0
                    do l = 1, 3  ! 全ての固有値が正、最小固有値が最大
                        if (minval(Eigen(l, :)) >= mintemp) then
                            index = l
                            mintemp = minval(Eigen(l, :))
                        endif
                    enddo

                    if (index > 0) then
                        Cpx(:, i-1, j, k) = Cptemp(:, index)
                        Cnx(:, i, j, k) = Cntemp(:, index)
                    else
                        Cpx(:, i-1, j, k) = C(:, i, j, k)
                        Cnx(:, i, j, k) = C(:, i, j, k)
                    endif
                    counter(index) = counter(index) + 1

                    ! index=3とどれだけずれているか確認
                    temp = temp + sum(Cptemp(:, index) - Cptemp(:, 3))/sum(Cptemp(:, 3))
                enddo
            enddo
        enddo
        ! write(*, '(E12.4)') temp/(NX*NY*NZ)

        Cpx(:, NX, :, :) = Cpx(:, 0, :, :)
        Cnx(:, 0, :, :) = Cnx(:, NX, :, :)
    end subroutine CpxCnx

    subroutine CpyCny(C, Cpy, Cny)
        real(8), intent(in) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Cpy(6, NX, 0:NY, NZ), Cny(6, NX, 0:NY, NZ)
        integer i, j, k, l, index
        real(8) Cptemp(6, 3), Cntemp(6, 3), Eigen(3, 6), mintemp

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Cptemp(:, 1) = C(:, i, j, k)*3/2 - C(:, i, j+1, k)/2
                    Cntemp(:, 1) = C(:, i, j, k)/2 + C(:, i, j+1, k)/2
                    Cptemp(:, 2) = C(:, i, j-1, k)/2 + C(:, i, j, k)/2
                    Cntemp(:, 2) = -C(:, i, j-1, k)/2 + C(:, i, j, k)*3/2
                    Cptemp(:, 3) = C(:, i, j-1, k)/4 + C(:, i, j, k) - C(:, i, j+1, k)/4
                    Cntemp(:, 3) = -C(:, i, j-1, k)/4 + C(:, i, j, k) + C(:, i, j+1, k)/4

                    Eigen(:, :) = 0.0d0
                    do l = 1, 3
                        if (eigen_method == 0) call Cardano(Cptemp(:, l), Eigen(l, 1), Eigen(l, 2), Eigen(l, 3))  ! CardanoでCpnの95%ぐらい時間かかってる。
                        if (eigen_method == 0) call Cardano(Cntemp(:, l), Eigen(l, 4), Eigen(l, 5), Eigen(l, 6))
                        if (eigen_method == 1) call Sylvester(Cptemp(:, l), Eigen(l, 1))
                        if (eigen_method == 1) call Sylvester(Cntemp(:, l), Eigen(l, 4))
                    enddo

                    index = 0
                    mintemp = 0.0d0
                    do l = 1, 3  ! 全ての固有値が正、最小固有値が最大
                        if (minval(Eigen(l, :)) >= mintemp) then
                            index = l
                            mintemp = minval(Eigen(l, :))
                        endif
                    enddo

                    if (index > 0) then
                        Cpy(:, i, j-1, k) = Cptemp(:, index)
                        Cny(:, i, j, k) = Cntemp(:, index)
                    else
                        Cpy(:, i, j-1, k) = C(:, i, j, k)
                        Cny(:, i, j, k) = C(:, i, j, k)
                    endif
                    counter(index) = counter(index) + 1
                enddo
            enddo
        enddo

        Cpy(:, :, NY, :) = Cpy(:, :, 0, :)
        Cny(:, :, 0, :) = Cny(:, :, NY, :)
    end subroutine CpyCny

    subroutine CpzCnz(C, Cpz, Cnz)
        real(8), intent(in) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Cpz(6, NX, NY, 0:NZ), Cnz(6, NX, NY, 0:NZ)
        integer i, j, k, l, index
        real(8) Cptemp(6, 3), Cntemp(6, 3), Eigen(3, 6), mintemp

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Cptemp(:, 1) = C(:, i, j, k)*3/2 - C(:, i, j, k+1)/2
                    Cntemp(:, 1) = C(:, i, j, k)/2 + C(:, i, j, k+1)/2
                    Cptemp(:, 2) = C(:, i, j, k-1)/2 + C(:, i, j, k)/2
                    Cntemp(:, 2) = -C(:, i, j, k-1)/2 + C(:, i, j, k)*3/2
                    Cptemp(:, 3) = C(:, i, j, k-1)/4 + C(:, i, j, k) - C(:, i, j, k+1)/4
                    Cntemp(:, 3) = -C(:, i, j, k-1)/4 + C(:, i, j, k) + C(:, i, j, k+1)/4

                    Eigen(:, :) = 0.0d0
                    do l = 1, 3
                        if (eigen_method == 0) call Cardano(Cptemp(:, l), Eigen(l, 1), Eigen(l, 2), Eigen(l, 3))  ! CardanoでCpnの95%ぐらい時間かかってる。
                        if (eigen_method == 0) call Cardano(Cntemp(:, l), Eigen(l, 4), Eigen(l, 5), Eigen(l, 6))
                        if (eigen_method == 1) call Sylvester(Cptemp(:, l), Eigen(l, 1))
                        if (eigen_method == 1) call Sylvester(Cntemp(:, l), Eigen(l, 4))
                    enddo

                    index = 0
                    mintemp = 0.0d0
                    do l = 1, 3  ! 全ての固有値が正、最小固有値が最大
                        if (minval(Eigen(l, :)) >= mintemp) then
                            index = l
                            mintemp = minval(Eigen(l, :))
                        endif
                    enddo

                    if (index > 0) then
                        Cpz(:, i, j, k-1) = Cptemp(:, index)
                        Cnz(:, i, j, k) = Cntemp(:, index)
                    else
                        Cpz(:, i, j, k-1) = C(:, i, j, k)
                        Cnz(:, i, j, k) = C(:, i, j, k)
                    endif
                    counter(index) = counter(index) + 1
                enddo
            enddo
        enddo

        Cpz(:, :, :, NZ) = Cpz(:, :, :, 0)
        Cnz(:, :, :, 0) = Cnz(:, :, :, NZ)
    end subroutine CpzCnz

    subroutine Cardano(a, re0, re1, re2)
        real(8), intent(in) :: a(6)
        real(8), intent(out) :: re0, re1, re2
        real(8) a0, a1, a2, p, q, t
        complex(8) u3, v3, w0, w1, w2, e0, e1, e2

        a2 = -(a(1) + a(4) + a(6))
        a1 = a(1)*a(4) + a(4)*a(6) + a(6)*a(1) - a(5)*a(5) - a(3)*a(3) - a(2)*a(2)
        a0 = -(a(1)*a(4)*a(6) + a(2)*a(5)*a(3) + a(2)*a(5)*a(3) - a(1)*a(5)*a(5) - a(2)*a(2)*a(6) - a(3)*a(4)*a(3))

        p = a1 - a2**2/3
        q = a0 - a1*a2/3 + 2*a2**3/27
        t = (q/2)**2 + (p/3)**3  ! 三次方程式の判別式
        ! write(*, *) t

        if (t > 0.0d0) then  ! 対称行列の固有値は全て実数で、tは必ず負になるので、このifを満たすものはない。
            ! u3 = -q/2 + cmplx(sqrt(t), 0.0d0, kind=8)
            ! v3 = -q/2 - cmplx(sqrt(t), 0.0d0, kind=8)
            ! stop 'complex eigenvalues'
            re0 = -1.0d0  ! 選ばれないようにする
            re1 = -1.0d0
            re2 = -1.0d0
        else
            u3 = -q/2 + cmplx(0.0d0, sqrt(-t), kind=8) ! t=0のときは解に虚数を含むが、倍精度の足し引きをするので、0はありえない。
            v3 = -q/2 - cmplx(0.0d0, sqrt(-t), kind=8)
        endif

        w0 = cmplx(1.0d0, 0.0d0, kind=8)
        w1 = cmplx(-0.5d0, sqrt(3.0d0)/2, kind=8)
        w2 = cmplx(-0.5d0, -sqrt(3.0d0)/2, kind=8)

        e0 = w0*u3**(1.0d0/3.0d0) + w0*v3**(1.0d0/3.0d0) - a2/3
        e1 = w1*u3**(1.0d0/3.0d0) + w2*v3**(1.0d0/3.0d0) - a2/3
        e2 = w2*u3**(1.0d0/3.0d0) + w1*v3**(1.0d0/3.0d0) - a2/3
        ! write(*, *) e0, e1, e2

        re0 = real(e0)
        re1 = real(e1)
        re2 = real(e2)
    end subroutine Cardano

    subroutine Sylvester(a, b)
        real(8), intent(in) :: a(6)
        real(8), intent(out) :: b  ! 正定値:1、それ以外:-1

        b = -1.0d0
        if (a(1)>0 .and. a(4)>0 .and. a(6)>0) then
            if (a(1)*a(4)-a(2)*a(2) > 0) then
                if (a(1)*a(4)*a(6) + a(2)*a(5)*a(3) + a(2)*a(5)*a(3) &
                  - a(1)*a(5)*a(5) - a(2)*a(2)*a(6) - a(3)*a(4)*a(3) > 0) then
                    b = 1.0d0
                    
                endif
            endif
        endif

    end subroutine Sylvester


    subroutine Cstar(Cpx, Cnx, Cpy, Cny, Cpz, Cnz, U, V, W, Cx)
        real(8), intent(in) :: Cpx(6, 0:NX, NY, NZ), Cnx(6, 0:NX, NY, NZ)
        real(8), intent(in) :: Cpy(6, NX, 0:NY, NZ), Cny(6, NX, 0:NY, NZ)
        real(8), intent(in) :: Cpz(6, NX, NY, 0:NZ), Cnz(6, NX, NY, 0:NZ)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Cx(6, NX, NY, NZ)
        integer i, j, k
        real(8) s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12
        real(8) u1, u2, v1, v2, w1, w2

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    u1 = U(i-1, j, k)
                    u2 = U(i, j, k)
                    v1 = V(i, j-1, k)
                    v2 = V(i, j, k)
                    w1 = W(i, j, k-1)
                    w2 = W(i, j, k)

                    s1 = (-u2 + abs(u2))*dt/(2*dX)
                    s2 = 1.0d0/6 + (-u2 - abs(u2))*dt/(2*dX)
                    s3 = 1.0d0/6 + (u1 - abs(u1))*dt/(2*dX)
                    s4 = (u1 + abs(u1))*dt/(2*dX)

                    s5 = (-v2 + abs(v2))*dt/(2*dY)
                    s6 = 1.0d0/6 + (-v2 - abs(v2))*dt/(2*dY)
                    s7 = 1.0d0/6 + (v1 - abs(v1))*dt/(2*dY)
                    s8 = (v1 + abs(v1))*dt/(2*dY)

                    s9 = (-w2 + abs(w2))*dt/(2*dZ)
                    s10 = 1.0d0/6 + (-w2 - abs(w2))*dt/(2*dZ)
                    s11 = 1.0d0/6 + (w1 - abs(w1))*dt/(2*dZ)
                    s12 = (w1 + abs(w1))*dt/(2*dZ)

                    Cx(:, i, j, k) = s3*Cpx(:, i-1, j, k) + s1*Cpx(:, i, j, k) + s4*Cnx(:, i-1, j, k) + s2*Cnx(:, i, j, k) &
                                    +s7*Cpy(:, i, j-1, k) + s5*Cpy(:, i, j, k) + s8*Cny(:, i, j-1, k) + s6*Cny(:, i, j, k) &
                                    +s11*Cpz(:, i, j, k-1) + s9*Cpz(:, i, j, k) + s12*Cnz(:, i, j, k-1) + s10*Cnz(:, i, j, k)
                    ! write(*, '(12e12.4)') s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12
                enddo
            enddo
        enddo

    end subroutine Cstar


    function f(C) result(y)
        real(8), intent(in) :: C(6) 
        real(8) y
        y = (Lp**2 - 3)/(Lp**2 - C(1)-C(4)-C(6))  ! FENE-Pモデル
        ! y = 1.0d0  ! Oldroyd-Bモデル
    end function f


    subroutine Lyapunov(Cx, U, V, W, C)
        real(8), intent(in) :: Cx(6, NX, NY, NZ)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        integer row, colum, itr
        real(8) A0(3, 3), b(6)
        real(8) x(6), x1(6), x2(6), l(6), l1(6), l2(6), Jac(6, 6), r(6)
        real(8), parameter :: dc = 1.0d-5  ! 3, 4, 5あたりがいい。
        integer, parameter :: itrmax = 100
        real(8), parameter :: eps = 1.0d-5
        newton_itr = 0

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    ! 定数部分を計算
                    ! i = 1; j = 1; k = 1
                    b(:) = Cx(:, i, j, k)
                    b(1) = b(1) + dt/Wi
                    b(4) = b(4) + dt/Wi
                    b(6) = b(6) + dt/Wi

                    A0(1, 1) = -dt*(-U(i-1, j, k) + U(i, j, k))/dX
                    A0(1, 2) = -dt*(-V(i-1, j-1, k) + V(i+1, j-1, k) - V(i-1, j, k) + V(i+1, j, k))/(4*dX)
                    A0(1, 3) = -dt*(-W(i-1, j, k-1) + W(i+1, j, k-1) - W(i-1, j, k) + W(i+1, j, k))/(4*dX)
                    A0(2, 1) = -dt*(-U(i-1, j-1, k) - U(i, j-1, k) + U(i-1, j+1, k) + U(i, j+1, k))/(4*dY)
                    A0(2, 2) = -dt*(-V(i, j-1, k) + V(i, j, k))/dY
                    A0(2, 3) = -dt*(-W(i, j-1, k-1) + W(i, j+1, k-1) - W(i, j-1, k) + W(i, j+1, k))/(4*dY)
                    A0(3, 1) = -dt*(-U(i-1, j, k-1) - U(i, j, k-1) + U(i-1, j, k+1) + U(i, j, k+1))/(4*dZ)
                    A0(3, 2) = -dt*(-V(i, j-1, k-1) - V(i, j, k-1) + V(i, j-1, k+1) + V(i, j, k+1))/(4*dZ)
                    A0(3, 3) = -dt*(-W(i, j, k-1) + W(i, j, k))/dZ
                    A0 = transpose(A0)  ! 転置して関数に渡す。転置しないとC11とC22の分布が入れ替わり、全体的に大きなCがでる
                    ! write(*, '(9e12.4)') A0(:, :)

                    ! ニュートン法
                    x(:) = C(:, i, j, k)  ! Cから探索開始
                    ! x(:) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
                    do itr = 1, itrmax
                        ! ヤコビアンの計算
                        do colum = 1, 6
                            x1(:) = x(:)
                            x2(:) = x(:)
                            x1(colum) = x1(colum) + dc
                            x2(colum) = x2(colum) - dc
                            call Lyapunov_func(A0, b, x1, l1)
                            call Lyapunov_func(A0, b, x2, l2)
                            do row = 1, 6
                                Jac(row, colum) = (l1(row) - l2(row))/(2.0d0*dc)
                            enddo
                        enddo

                        call Lyapunov_func(A0, b, x, l)  ! f(xi)を計算
                        ! write(*, '(6e12.4)') l(:)
                        
                        call gauss(Jac, l, r)  ! ガウスの消去法
                        x(:) = x(:) - r(:)
                        ! write(*, '(6e12.4)') sum(r**2)
                        ! write(*, *) itr, sum(r**2), sum(x**2), sqrt(sum(r**2)/sum(x**2))
                        ! write(*, *) ''

                        if (sqrt(sum(r**2)/sum(x**2)) < eps) exit
                    enddo
                    newton_itr = newton_itr + itr
                    C(:, i, j, k) = x(:)
                enddo
            enddo
        enddo
        ! write(*, '(6e12.4)') C(:, 1, 1, 1)
        ! write(*, '(6e12.4)') C(:, NX, NY, NZ)
        call C_PBM(C)
    end subroutine Lyapunov


    subroutine Lyapunov_func(A0, b, x, l)  ! l(x)=Mx-bを計算
        real(8), intent(in) :: A0(3, 3), b(6), x(6)
        real(8), intent(out) :: l(6)
        real(8) A(3, 3)

        A(:, :) = A0(:, :)
        A(1, 1) = A(1, 1) + (1.0d0 + f(x)*dt/Wi)/2.0d0
        A(2, 2) = A(2, 2) + (1.0d0 + f(x)*dt/Wi)/2.0d0
        A(3, 3) = A(3, 3) + (1.0d0 + f(x)*dt/Wi)/2.0d0

        l(1) = (A(1, 1) + A(1, 1))*x(1) + A(1, 2)*x(2) + A(1, 3)*x(3) + A(1, 2)*x(2) + A(1, 3)*x(3) - b(1)
        l(2) = A(2, 1)*x(1) + (A(2, 2) + A(1, 1))*x(2) + A(2, 3)*x(3) + A(1, 2)*x(4) + A(1, 3)*x(5) - b(2)
        l(3) = A(3, 1)*x(1) + A(3, 2)*x(2) + (A(3, 3) + A(1, 1))*x(3) + A(1, 2)*x(5) + A(1, 3)*x(6) - b(3)
        l(4) = A(2, 1)*x(2) + A(2, 1)*x(2) + (A(2, 2) + A(2, 2))*x(4) + A(2, 3)*x(5) + A(2, 3)*x(5) - b(4)
        l(5) = A(2, 1)*x(3) + A(3, 1)*x(2) + A(3, 2)*x(4) + (A(3, 3) + A(2, 2))*x(5) + A(2, 3)*x(6) - b(5)
        l(6) = A(3, 1)*x(3) + A(3, 2)*x(5) + A(3, 1)*x(3) + A(3, 2)*x(5) + (A(3, 3) + A(3, 3))*x(6) - b(6)

    end subroutine Lyapunov_func


    subroutine gauss(a, b, x)  ! ガウスの消去法（ピボット選択あり）
        integer, parameter :: n = 6
        real(8), intent(in) :: a(n, n), b(n)
        real(8), intent(out) :: x(n)
        integer i, k, maxidx
        real(8) maxabs, temp(n+1), Ah(n, n+1)
        Ah(:, :n) = a(:, :)
        Ah(:, n+1) = b(:)

        do k = 1, n
            maxidx = maxloc(abs(Ah(k:n, k)), 1) + k-1  ! ピボット選択して並び替え
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

            Ah(k, k+1:n+1) = Ah(k, k+1:n+1) / Ah(k, k)  ! 先にx1を計算
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


    subroutine C_PBM(C)
        real(8), intent(inout) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        C(:, 0, :, :) = C(:, NX, :, :)  ! 周期境界条件
        C(:, NX+1, :, :) = C(:, 1, :, :)
        C(:, :, 0, :) = C(:, :, NY, :)
        C(:, :, NY+1, :) = C(:, :, 1, :)
        C(:, :, :, 0) = C(:, :, :, NZ)
        C(:, :, :, NZ+1) = C(:, :, :, 1)
    end subroutine C_PBM

    subroutine PBM(A)
        real(8), intent(inout) :: A(0:NX+1, 0:NY+1, 0:NZ+1)
        A(0, :, :) = A(NX, :, :)
        A(NX+1, :, :) = A(1, :, :)
        A(:, 0, :) = A(:, NY, :)
        A(:, NY+1, :) = A(:, 1, :)
        A(:, :, 0) = A(:, :, NZ)
        A(:, :, NZ+1) = A(:, :, 1)
    end subroutine PBM


    subroutine polymer_stress(C, Tx, Ty, Tz)
        real(8), intent(in) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Tx(1:NX, 1:NY, 1:NZ), Ty(1:NX, 1:NY, 1:NZ), Tz(1:NX, 1:NY, 1:NZ)
        integer i, j, k
        real(8) C0(6, 0:NX+1, 0:NY+1, 0:NZ+1), E(6)

        E(:) = (/1.0, 0.0, 0.0, 1.0, 0.0, 1.0/)
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    C0(:, i, j, k) = (f(C(:, i, j, k))*C(:, i, j, k) - E(:))
                enddo
            enddo
        enddo
        call C_PBM(C0)

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Tx(i, j, k) = (-C0(1, i, j, k) + C0(1, i+1, j, k))/dX &
                                 +(-C0(2, i, j-1, k) - C0(2, i+1, j-1, k) + C0(2, i, j+1, k) + C0(2, i+1, j+1, k))/(4*dY) &
                                 +(-C0(3, i, j, k-1) - C0(3, i+1, j, k-1) + C0(3, i, j, k+1) + C0(3, i+1, j, k+1))/(4*dZ)
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Ty(i, j, k) = (-C0(2, i-1, j, k) + C0(2, i+1, j, k) - C0(2, i-1, j+1, k) + C0(2, i+1, j+1, k))/(4*dX) &
                                 +(-C0(4, i, j, k) + C0(4, i, j+1, k))/dY &
                                 +(-C0(5, i, j, k-1) - C0(5, i, j+1, k-1) + C0(5, i, j, k+1) + C0(5, i, j+1, k+1))/(4*dZ)
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Tz(i, j, k) = (-C0(3, i-1, j, k) + C0(3, i+1, j, k) - C0(3, i-1, j, k+1) + C0(3, i+1, j, k+1))/(4*dX) &
                                 +(-C0(5, i, j-1, k) + C0(5, i, j+1, k) - C0(5, i, j-1, k+1) + C0(5, i, j+1, k+1))/(4*dY) &
                                 +(-C0(6, i, j, k) + C0(6, i, j, k+1))/dZ
                enddo
            enddo
        enddo
        Tx(:, :, :) = (1.0d0-beta)/(Re*Wi)*Tx(:, :, :)
        Ty(:, :, :) = (1.0d0-beta)/(Re*Wi)*Ty(:, :, :)
        Tz(:, :, :) = (1.0d0-beta)/(Re*Wi)*Tz(:, :, :)
    end subroutine polymer_stress


    subroutine convection(U, V, W, Ax, Ay, Az)  ! 発散型
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        integer i, j, k
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Ax(i, j, k) = (-((U(i-1, j, k)+U(i, j, k))/2)**2 + ((U(i, j, k)+U(i+1, j, k))/2)**2)/dX &
                                 +(-(V(i, j-1, k)+V(i+1, j-1, k))/2*(U(i, j-1, k)+U(i, j, k))/2 &
                                  + (V(i, j, k)+V(i+1, j, k))/2*(U(i, j, k)+U(i, j+1, k))/2)/dY &
                                 +(-(W(i, j, k-1)+W(i+1, j, k-1))/2*(U(i, j, k-1)+U(i, j, k))/2 &
                                  + (W(i, j, k)+W(i+1, j, k))/2*(U(i, j, k)+U(i, j, k+1))/2)/dZ
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Ay(i, j, k) = (-(U(i-1, j, k)+U(i-1, j+1, k))/2*(V(i-1, j, k)+V(i, j, k))/2 &
                                  + (U(i, j, k)+U(i, j+1, k))/2*(V(i, j, k)+V(i+1, j, k))/2)/dX &
                                 +(-((V(i, j-1, k)+V(i, j, k))/2)**2 + ((V(i, j, k)+V(i, j+1, k))/2)**2)/dY &
                                 +(-(W(i, j, k-1)+W(i, j+1, k-1))/2*(V(i, j, k-1)+V(i, j, k))/2 &
                                  + (W(i, j, k)+W(i, j+1, k))/2*(V(i, j, k)+V(i, j, k+1))/2)/dZ
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Az(i, j, k) = (-(U(i-1, j, k)+U(i-1, j, k+1))/2*(W(i-1, j, k)+W(i, j, k))/2 &
                                  + (U(i, j, k)+U(i, j, k+1))/2*(W(i, j, k)+W(i+1, j, k))/2)/dX &
                                 +(-(V(i, j-1, k)+V(i, j-1, k+1))/2*(W(i, j-1, k)+W(i, j, k))/2 &
                                  + (V(i, j, k)+V(i, j, k+1))/2*(W(i, j, k)+W(i, j+1, k))/2)/dY &
                                 +(-((W(i, j, k-1)+W(i, j, k))/2)**2 + ((W(i, j, k)+W(i, j, k+1))/2)**2)/dZ
                enddo
            enddo
        enddo
        Ax(:, :, :) = -Ax(:, :, :)
        Ay(:, :, :) = -Ay(:, :, :)
        Az(:, :, :) = -Az(:, :, :)
    end subroutine convection

    subroutine viscous(U, V, W, Bx, By, Bz)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        integer i, j, k
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Bx(i, j, k) = (U(i-1, j, k) - 2*U(i, j, k) + U(i+1, j, k)) / dX**2 &
                                 +(U(i, j-1, k) - 2*U(i, j, k) + U(i, j+1, k)) / dY**2 &
                                 +(U(i, j, k-1) - 2*U(i, j, k) + U(i, j, k+1)) / dZ**2
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    By(i, j, k) = (V(i-1, j, k) - 2*V(i, j, k) + V(i+1, j, k)) / dX**2 &
                                 +(V(i, j-1, k) - 2*V(i, j, k) + V(i, j+1, k)) / dY**2 &
                                 +(V(i, j, k-1) - 2*V(i, j, k) + V(i, j, k+1)) / dZ**2
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Bz(i, j, k) = (W(i-1, j, k) - 2*W(i, j, k) + W(i+1, j, k)) / dX**2 &
                                 +(W(i, j-1, k) - 2*W(i, j, k) + W(i, j+1, k)) / dY**2 &
                                 +(W(i, j, k-1) - 2*W(i, j, k) + W(i, j, k+1)) / dZ**2
                enddo
            enddo
        enddo
        Bx(:, :, :) = beta/Re*Bx(:, :, :)
        By(:, :, :) = beta/Re*By(:, :, :)
        Bz(:, :, :) = beta/Re*Bz(:, :, :)
    end subroutine viscous

    subroutine navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, &
                      Bx, By, Bz, Bx0, By0, Bz0, Tx, Ty, Tz, Tx0, Ty0, Tz0, Fx, Fy, Fz, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Bx0(1:NX, 1:NY, 1:NZ), By0(1:NX, 1:NY, 1:NZ), Bz0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Tx(1:NX, 1:NY, 1:NZ), Ty(1:NX, 1:NY, 1:NZ), Tz(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Tx0(1:NX, 1:NY, 1:NZ), Ty0(1:NX, 1:NY, 1:NZ), Tz0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
        integer, intent(in) :: step
        integer i, j, k

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :, :) = Ax(:, :, :)
            Ay0(:, :, :) = Ay(:, :, :)
            Az0(:, :, :) = Az(:, :, :)
            Bx0(:, :, :) = Bx(:, :, :)
            By0(:, :, :) = By(:, :, :)
            Bz0(:, :, :) = Bz(:, :, :)
            Tx0(:, :, :) = Tx(:, :, :)
            Ty0(:, :, :) = Ty(:, :, :)
            Tz0(:, :, :) = Tz(:, :, :)
        endif

        ! NS方程式で速度の予測
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Up(i, j, k) = U(i, j, k) - dt*(-P(i, j, k)+P(i+1, j, k))/dX &
                                 +dt*(3.0d0*(Ax(i, j, k) + Bx(i, j, k)) - (Ax0(i, j, k) + Bx0(i, j, k)))/2.0d0 &
                                 +dt*(Tx(i, j, k) + Tx0(i, j, k))/2.0d0 &
                                 +dt*Fx(i, j, k)
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Vp(i, j, k) = V(i, j, k) - dt*(-P(i, j, k)+P(i, j+1, k))/dY &
                                 +dt*(3.0d0*(Ay(i, j, k) + By(i, j, k)) - (Ay0(i, j, k) + By0(i, j, k)))/2.0d0 &
                                 +dt*(Ty(i, j, k) + Ty0(i, j, k))/2.0d0 &
                                 +dt*Fy(i, j, k)
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Wp(i, j, k) = W(i, j, k) - dt*(-P(i, j, k)+P(i, j, k+1))/dZ &
                                 +dt*(3.0d0*(Az(i, j, k) + Bz(i, j, k)) - (Az0(i, j, k) + Bz0(i, j, k)))/2.0d0 &
                                 +dt*(Tz(i, j, k) + Tz0(i, j, k))/2.0d0 &
                                 +dt*Fz(i, j, k)
                enddo
            enddo
        enddo

        ! Phiを求めるときに0番目の値も必要
        call PBM(Up)
        call PBM(Vp)
        call PBM(Wp)

        ! n-1ステップ目の保存
        Ax0(:, :, :) = Ax(:, :, :)
        Ay0(:, :, :) = Ay(:, :, :)
        Az0(:, :, :) = Az(:, :, :)
        Bx0(:, :, :) = Bx(:, :, :)
        By0(:, :, :) = By(:, :, :)
        Bz0(:, :, :) = Bz(:, :, :)
        Tx0(:, :, :) = Tx(:, :, :)
        Ty0(:, :, :) = Ty(:, :, :)
        Tz0(:, :, :) = Tz(:, :, :)
    end subroutine navier

    subroutine poisson(Up, Vp, Wp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k, l, itr
        real(8) BXM, BXP, BYM, BYP, BZM, BZP, B0
        real(8) er, er0
        real(8) E(1:NX, 1:NY, 1:NZ), Q(1:NX, 1:NY, 1:NZ)
        real(8), parameter :: SOR = 2.0d0 / (1 + sin(PI*dX))
        real(8), parameter :: eps = 1.0d-5
        integer, parameter :: itrmax = 10000
        poisson_itr = 0

        ! ポアソン方程式で用いる定数の計算(今回は等間隔な格子)
        BXM = 1.0d0 / dX**2
        BXP = 1.0d0 / dX**2
        BYM = 1.0d0 / dY**2
        BYP = 1.0d0 / dY**2
        BZM = 1.0d0 / dZ**2
        BZP = 1.0d0 / dZ**2
        B0 = BXM + BXP + BYM + BYP + BZM + BZP
        ! 右辺の計算
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = ((-Up(i-1, j, k) + Up(i, j, k)) / dX &
                                 +(-Vp(i, j-1, k) + Vp(i, j, k)) / dY &
                                 +(-Wp(i, j, k-1) + Wp(i, j, k)) / dZ) / dt
                enddo
            enddo
        enddo
        er0 = sum(Q**2)
        if (er0 == 0.0d0) then  ! Qが全て0のとき計算できない
            er0 = sum(Q**2)
        endif

        ! SOR法
        ! Phi(:, :, :) = 0.0d0  ! 前回のPhiから計算したほうが早く収束する
        do itr = 1, itrmax
            do k = 1, NZ
                do j = 1, NY
                    do l = 1, 2  ! 奇数と偶数に分けて計算する
                        do i = l, NX, 2
                            E(i, j, k) = BZM*Phi(i, j, k-1) + BYM*Phi(i, j-1, k) + BXM*Phi(i-1, j, k) - B0*Phi(i, j, k) &
                                        +BXP*Phi(i+1, j, k) + BYP*Phi(i, j+1, k) + BZP*Phi(i, j, k+1) - Q(i, j, k)
                            Phi(i, j, k) = Phi(i, j, k) + SOR*E(i, j, k) / B0
                        enddo
                    enddo
                enddo
            enddo
            call PBM(Phi)

            er = sum(E**2)
            if (sqrt(er / er0) < eps) exit
        enddo
        poisson_itr = poisson_itr + itr
        ! write(*, '(16e12.4)') Phi(1:NX, 1, 1)
        ! write(*, *) 'SOR', SOR, 'itr =', itr, 'er', er, 'er0', er0
    end subroutine poisson

    subroutine march(Up, Vp, Wp, U, V, W, Phi, P)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    U(i, j, k) = Up(i, j, k) - dt*(-Phi(i, j, k)+Phi(i+1, j, k))/dX
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    V(i, j, k) = Vp(i, j, k) - dt*(-Phi(i, j, k)+Phi(i, j+1, k))/dY
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    W(i, j, k) = Wp(i, j, k) - dt*(-Phi(i, j, k)+Phi(i, j, k+1))/dZ
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    P(i, j, k) = P(i, j, k) + Phi(i, j, k)
                enddo
            enddo
        enddo

        call PBM(U)
        call PBM(V)
        call PBM(W)
        call PBM(P)
    end subroutine march

    subroutine taylor_debag(U, V, W, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k
        real(8) err
        real(8) U0(0:NX+1, 0:NY+1, 0:NZ+1), V0(0:NX+1, 0:NY+1, 0:NZ+1), W0(0:NX+1, 0:NY+1, 0:NZ+1)
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    U0(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)
                    V0(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)
                enddo
            enddo
        enddo
        W0(:, :, :) = 0.0d0
        U0(:, :, :) = U0(:, :, :) * exp(-2.0d0/Re*step*dt)  ! 解析解
        V0(:, :, :) = V0(:, :, :) * exp(-2.0d0/Re*step*dt)
        W0(:, :, :) = W0(:, :, :) * exp(-2.0d0/Re*step*dt)

        i = NX/4
        j = NY/4
        k = NZ/4
        err = (sum((U0(1:NX, 1:NY, 1:NZ)-U(1:NX, 1:NY, 1:NZ))**2) &
              +sum((V0(1:NX, 1:NY, 1:NZ)-V(1:NX, 1:NY, 1:NZ))**2) &
              +sum((W0(1:NX, 1:NY, 1:NZ)-W(1:NX, 1:NY, 1:NZ))**2))/(NX*NY*NZ)
        open(10, file=dir//'debag2.d', position='append')
        write(10, *) step, U0(i, j, k), U(i, j, k), U0(i, j, k) - U(i, j, k), sqrt(err)
        close(10)
    end subroutine taylor_debag


    subroutine get_data(U, V, W, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k, l
        real(8) D(3, 3), Omega(3), S(3), Qti
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        ! open(10, file='data/y_'//str//'.d')
        open(10, file=dir//'all_'//str//'.d')

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    D(1, 1) = (U(i, j, k) - U(i-1, j, k))/dX
                    D(1, 2) = (U(i, j+1, k) - U(i, j-1, k) + U(i-1, j+1, k) - U(i-1, j-1, k))/(4*dY)
                    D(1, 3) = (U(i, j, k+1) - U(i, j, k-1) + U(i-1, j, k+1) - U(i-1, j, k-1))/(4*dZ)
                    D(2, 1) = (V(i+1, j, k) - V(i-1, j, k) + V(i+1, j-1, k) - V(i-1, j-1, k))/(4*dX)
                    D(2, 2) = (V(i, j, k) - V(i, j-1, k))/dY
                    D(2, 3) = (V(i, j, k+1) - V(i, j, k-1) + V(i, j-1, k+1) - V(i, j-1, k-1))/(4*dZ)
                    D(3, 1) = (W(i+1, j, k) - W(i-1, j, k) + W(i+1, j, k-1) - W(i-1, j, k-1))/(4*dX)
                    D(3, 2) = (W(i, j+1, k) - W(i, j-1, k) + W(i, j+1, k-1) - W(i, j-1, k-1))/(4*dY)
                    D(3, 3) = (W(i, j, k) - W(i, j, k-1))/dZ
                    D(:, :) = D(:, :)*U_C

                    Omega(1) = (D(3, 2) - D(2, 3))
                    Omega(2) = (D(1, 3) - D(3, 1))
                    Omega(3) = (D(2, 1) - D(1, 2))
                    S(1) = (D(3, 2) + D(2, 3))
                    S(2) = (D(1, 3) + D(3, 1))
                    S(3) = (D(2, 1) + D(1, 2))
                    Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                    write(10, '(11e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
                                        (U(i-1, j, k)+U(i, j, k))/2*U_C, &
                                        (V(i, j-1, k)+V(i, j, k))/2*U_C, &
                                        (W(i, j, k-1)+W(i, j, k))/2*U_C, &
                                        Omega(3), sum(Omega**2)/2, S(3), sum(S**2)/2, Qti
                enddo
            enddo
        enddo
        close(10)
    end subroutine get_data

    subroutine get_data_xy(U, V, W, C, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k
        real(8) D(3, 3), omega, sendan, Qti, trC
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換
        open(10, file=dir//'z_'//str//'.d')

        k = NZ/2
        do j = 1, NY
            do i = 1, NX
                D(1, 1) = (U(i, j, k) - U(i-1, j, k))/dX
                D(1, 2) = (U(i, j+1, k) - U(i, j-1, k) + U(i-1, j+1, k) - U(i-1, j-1, k))/(4*dY)
                D(1, 3) = (U(i, j, k+1) - U(i, j, k-1) + U(i-1, j, k+1) - U(i-1, j, k-1))/(4*dZ)
                D(2, 1) = (V(i+1, j, k) - V(i-1, j, k) + V(i+1, j-1, k) - V(i-1, j-1, k))/(4*dX)
                D(2, 2) = (V(i, j, k) - V(i, j-1, k))/dY
                D(2, 3) = (V(i, j, k+1) - V(i, j, k-1) + V(i, j-1, k+1) - V(i, j-1, k-1))/(4*dZ)
                D(3, 1) = (W(i+1, j, k) - W(i-1, j, k) + W(i+1, j, k-1) - W(i-1, j, k-1))/(4*dX)
                D(3, 2) = (W(i, j+1, k) - W(i, j-1, k) + W(i, j+1, k-1) - W(i, j-1, k-1))/(4*dY)
                D(3, 3) = (W(i, j, k) - W(i, j, k-1))/dZ
                D(:, :) = D(:, :)*U_C

                omega = (D(2, 1) - D(1, 2))
                sendan = (D(2, 1) + D(1, 2))
                Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                trC = C(1, i, j, k) + C(4, i, j, k) + C(6, i, j, k)

                write(10, '(18e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
                                       (U(i-1, j, k)+U(i, j, k))/2*U_C, &
                                       (V(i, j-1, k)+V(i, j, k))/2*U_C, &
                                       (W(i, j, k-1)+W(i, j, k))/2*U_C, &
                                       omega, D(1, 1), sendan, D(2, 2), Qti, trC, &
                                       C(1, i, j, k), C(2, i, j, k), C(3, i, j, k), C(4, i, j, k), C(5, i, j, k), C(6, i, j, k)
            enddo
        enddo
        close(10)
    end subroutine get_data_xy

    subroutine get_data_xz(U, V, W, C, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k, l
        real(8) D(3, 3), Omega(3), S(3), Qti, trC
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換
        open(10, file=dir//'y_'//str//'.d')

        j = NY/2
        do k = 1, NZ
            do i = 1, NX
                D(1, 1) = (U(i, j, k) - U(i-1, j, k))/dX
                D(1, 2) = (U(i, j+1, k) - U(i, j-1, k) + U(i-1, j+1, k) - U(i-1, j-1, k))/(4*dY)
                D(1, 3) = (U(i, j, k+1) - U(i, j, k-1) + U(i-1, j, k+1) - U(i-1, j, k-1))/(4*dZ)
                D(2, 1) = (V(i+1, j, k) - V(i-1, j, k) + V(i+1, j-1, k) - V(i-1, j-1, k))/(4*dX)
                D(2, 2) = (V(i, j, k) - V(i, j-1, k))/dY
                D(2, 3) = (V(i, j, k+1) - V(i, j, k-1) + V(i, j-1, k+1) - V(i, j-1, k-1))/(4*dZ)
                D(3, 1) = (W(i+1, j, k) - W(i-1, j, k) + W(i+1, j, k-1) - W(i-1, j, k-1))/(4*dX)
                D(3, 2) = (W(i, j+1, k) - W(i, j-1, k) + W(i, j+1, k-1) - W(i, j-1, k-1))/(4*dY)
                D(3, 3) = (W(i, j, k) - W(i, j, k-1))/dZ
                D(:, :) = D(:, :)*U_C

                Omega(1) = (D(3, 2) - D(2, 3))
                Omega(2) = (D(1, 3) - D(3, 1))
                Omega(3) = (D(2, 1) - D(1, 2))

                trC = C(1, i, j, k) + C(4, i, j, k) + C(6, i, j, k)

                write(10, '(18e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
                                        (U(i-1, j, k)+U(i, j, k))/2*U_C, &
                                        (V(i, j-1, k)+V(i, j, k))/2*U_C, &
                                        (W(i, j, k-1)+W(i, j, k))/2*U_C, &
                                        Omega(2), 0.0d0, 0.0d0, 0.0d0, 0.0d0, trC, &
                                        C(1, i, j, k), C(2, i, j, k), C(3, i, j, k), C(4, i, j, k), C(5, i, j, k), C(6, i, j, k)
            enddo
        enddo
        close(10)
        l = 1  ! 使ってません対策
    end subroutine get_data_xz

    subroutine logging(U, V, W, C, step, t_start)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        real(8), intent(in) :: t_start
        integer i, j, k, count
        real(8) trC(NX, NY, NZ), re0, re1, re2
        real(8) t_temp
        real(8) K_energy
        integer(8) t_pre, hour, min, sec
        real(8) U_tmp(0:NX+1, 0:NY+1, 0:NZ+1), V_tmp(0:NX+1, 0:NY+1, 0:NZ+1), W_tmp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) U_grad(NX, NY, NZ), V_grad(NX, NY, NZ), W_grad(NX, NY, NZ)
        real(8) U_rms, lamda, Re_lamda, tmp, D(3, 3), S(3, 3), epsilon, eta

        if (mod(step, Dstep) == 0) then
            ! テイラー長レイノルズ数
            U_tmp(:, :, :) = U(:, :, :) - sum(U(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
            V_tmp(:, :, :) = V(:, :, :) - sum(V(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
            W_tmp(:, :, :) = W(:, :, :) - sum(W(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
            U_tmp(:, :, :) = U_tmp(:, :, :) * U_C
            V_tmp(:, :, :) = V_tmp(:, :, :) * U_C
            W_tmp(:, :, :) = W_tmp(:, :, :) * U_C
            do k = 1, NZ
                do j = 1, NY
                    do i = 1, NX
                        U_grad(i, j, k) = (U_tmp(i, j, k) - U_tmp(i-1, j, k))/dX
                        V_grad(i, j, k) = (V_tmp(i, j, k) - V_tmp(i, j-1, k))/dY
                        W_grad(i, j, k) = (W_tmp(i, j, k) - W_tmp(i, j, k-1))/dZ
                    enddo
                enddo
            enddo
            U_rms = sqrt((sum(U_tmp(1:NX, 1:NY, 1:NZ)**2) &
                        + sum(V_tmp(1:NX, 1:NY, 1:NZ)**2) &
                        + sum(W_tmp(1:NX, 1:NY, 1:NZ)**2))/(3.0d0*NX*NY*NZ))
            lamda = sqrt(sum(U_tmp(1:NX, 1:NY, 1:NZ)**2)/sum(U_grad**2) &  ! 空間平均同士で割るので格子点数で割らなくてよい
                    + sum(V_tmp(1:NX, 1:NY, 1:NZ)**2)/sum(V_grad**2) &
                    + sum(W_tmp(1:NX, 1:NY, 1:NZ)**2)/sum(W_grad**2))
            Re_lamda = U_rms*lamda/nu

            ! コルモゴロフ長
            tmp = 0.0d0
            do k = 1, NZ
                do j = 1, NY
                    do i = 1, NX
                        D(1, 1) = (U(i, j, k) - U(i-1, j, k))/dX
                        D(1, 2) = (U(i, j+1, k) - U(i, j-1, k) + U(i-1, j+1, k) - U(i-1, j-1, k))/(4*dY)
                        D(1, 3) = (U(i, j, k+1) - U(i, j, k-1) + U(i-1, j, k+1) - U(i-1, j, k-1))/(4*dZ)
                        D(2, 1) = (V(i+1, j, k) - V(i-1, j, k) + V(i+1, j-1, k) - V(i-1, j-1, k))/(4*dX)
                        D(2, 2) = (V(i, j, k) - V(i, j-1, k))/dY
                        D(2, 3) = (V(i, j, k+1) - V(i, j, k-1) + V(i, j-1, k+1) - V(i, j-1, k-1))/(4*dZ)
                        D(3, 1) = (W(i+1, j, k) - W(i-1, j, k) + W(i+1, j, k-1) - W(i-1, j, k-1))/(4*dX)
                        D(3, 2) = (W(i, j+1, k) - W(i, j-1, k) + W(i, j+1, k-1) - W(i, j-1, k-1))/(4*dY)
                        D(3, 3) = (W(i, j, k) - W(i, j, k-1))/dZ
                        D(:, :) = D(:, :)*U_C
                        S(1, 1) = (D(1, 1) + D(1, 1))
                        S(1, 2) = (D(1, 2) + D(2, 1))
                        S(1, 3) = (D(1, 3) + D(3, 1))
                        S(2, 1) = (D(2, 1) + D(1, 2))
                        S(2, 2) = (D(2, 2) + D(2, 2))
                        S(2, 3) = (D(2, 3) + D(3, 2))
                        S(3, 1) = (D(3, 1) + D(1, 3))
                        S(3, 2) = (D(3, 2) + D(2, 3))
                        S(3, 3) = (D(3, 3) + D(3, 3))
                        tmp = tmp + sum(S**2)
                    enddo
                enddo
            enddo
            epsilon = 0.5d0*nu*tmp/(NX*NY*NZ)
            eta = (nu**3/epsilon)**0.25d0
            ! write(*, '(12e12.4)') Re_lamda, epsilon, eta, dX
        endif

        if (mod(step, Dstep) == 0) then
            write(*, '(a, I6)', advance='no') 'step:', step

            ! 運動エネルギーとエネルギーカスケードの総和
            K_energy = sum(U(1:NX, 1:NY, 1:NZ)**2) + sum(V(1:NX, 1:NY, 1:NZ)**2) + sum(W(1:NX, 1:NY, 1:NZ)**2)
            K_energy = K_energy*U_C**2/2
            K_energy = K_energy/(NX*NY*NZ)
            write(*, '(a, e12.4)', advance='no') '  | K_energy:', K_energy

            write(*, '(a, F6.3)', advance='no') '  | index 0:', counter(0)*1.0d0/sum(counter)
            write(*, '(a, F6.3)', advance='no') '  1:', counter(1)*1.0d0/sum(counter)
            write(*, '(a, F6.3)', advance='no') '  2:', counter(2)*1.0d0/sum(counter)
            write(*, '(a, F6.3)', advance='no') '  3:', counter(3)*1.0d0/sum(counter)
            ! write(*, '(a, F6.2)', advance='no') '  | newton_itr:', newton_itr*1.0d0/(NX*NY*NZ)
            ! write(*, '(a, I5)', advance='no') '  | poisson_itr:', poisson_itr

            trC(:, :, :) = C(1, 1:NX, 1:NY, 1:NZ)+C(4, 1:NX, 1:NY, 1:NZ)+C(6, 1:NX, 1:NY, 1:NZ)
            write(*, '(a, F7.3, a, F8.3)', advance='no') '  |', minval(trC(:, :, :)), '< trC <', maxval(trC(:, :, :))
            ! write(*, '(a, F7.3, a, F7.3)', advance='no') '  |', minval(C(1, :, :, :)), '< Cxx <', maxval(C(1, :, :, :))
            ! write(*, '(a, F7.3, a, F7.3)', advance='no') '  |', minval(C(4, :, :, :)), '< Cyy <', maxval(C(4, :, :, :))
            ! write(*, '(a, F7.3, a, F7.3)', advance='no') '  |', minval(C(6, :, :, :)), '< Czz <', maxval(C(6, :, :, :))
            ! write(*, *) ''

            count = 0
            do k = 1, NZ
                do j = 1, NY
                    do i = 1, NX
                        ! call Cardano(C(:, i, j, k), re0, re1, re2)
                        ! if (re0 < 0.0d0 .or. re1 < 0.0d0 .or. re2 < 0.0d0) count = count + 1
                        call Sylvester(C(:, i, j, k), re0)
                        if (re0 > 0.0d0) count = count + 1
                    enddo
                enddo
            enddo
            write(*, '(a, F7.3)', advance='no') '  | SPD:', count*1.0d0/(NX*NY*NZ)

            ! write(*, '(a, e12.4)', advance='no') '  | Re_lamda:', Re_lamda
            ! write(*, '(a, e12.4)', advance='no') '  | epsilon:', epsilon
            ! write(*, '(a, e12.4)', advance='no') '  | eta:', eta

            call cpu_time(t_temp)
            t_pre = int((t_temp-t_start)*(Nstep-step)/step)
            hour = t_pre / 3600
            min = mod(t_pre, 3600) / 60
            sec = mod(t_pre, 60)
            write(*, '(a, I3, a, I2, a, I2)', advance='no') '  | time_left:', hour, ':', min, ':', sec

            write(*, *) ''
        endif
    end subroutine logging


    subroutine input(U, V, W, P, C)  ! initを実行した後に書く
        real(8), intent(out) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: P(0:NX+1, 0:NY+1, 0:NZ+1), C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        character(8) str
        real(8) K_energy
        write(str, '(I8.8)') input_step  ! 数値を文字列に変換
        open(10, file=dir//str//'.d')
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    read(10, '(10e12.4)') U(i, j, k), V(i, j, k), W(i, j, k), P(i, j, k), C(:, i, j, k)
                enddo
            enddo
        enddo
        close(10)
        call PBM(U)
        call PBM(V)
        call PBM(W)
        call PBM(P)

        K_energy = sum(U(1:NX, 1:NY, 1:NZ)**2) + sum(V(1:NX, 1:NY, 1:NZ)**2) + sum(W(1:NX, 1:NY, 1:NZ)**2)
        K_energy = K_energy*U_C**2/2
        K_energy = K_energy/(NX*NY*NZ)
        write(*, '(a, a, e12.4)') 'input_file='//dir//str//'.d', '  | K_energy:', K_energy

    end subroutine input

    subroutine output(U, V, W, P, C, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1), C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換
        open(10, file=dir//str//'.d')

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    write(10, '(10e12.4)') U(i, j, k), V(i, j, k), W(i, j, k), P(i, j, k), C(:, i, j, k)
                enddo
            enddo
        enddo
        close(10)
    end subroutine output

    subroutine mk_dir(outdir)
        implicit none
        character(*), intent(in) :: outdir
        character(256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dir

end module smac

module fft
    use smac
    implicit none
    include 'fftw3.f'
contains
    subroutine fft_solve(Q, LHS, Phi)
        real(8), intent(in) :: Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)
        real(8), intent(out) :: Phi(1:NX, 1:NY, 1:NZ)
        complex(8) Q_hat(1:NX/2+1, 1:NY, 1:NZ), Phi_hat(1:NX/2+1, 1:NY, 1:NZ)
        integer i, j, k
        integer(8) plan

        ! 左辺をフーリエ変換
        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, Q, Q_hat, FFTW_ESTIMATE)
        call dfftw_execute(plan, Q, Q_hat)
        call dfftw_destroy_plan(plan)
        
        Q_hat(:, :, :) = Q_hat(:, :, :)/(NX*NY*NZ)

        ! 定数で割ることで左辺を求める
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX/2+1
                    if (abs(LHS(i, j, k)) < 1.0d-16) then  ! マシン零以下
                        Phi_hat(i, j, k) = (0.0d0, 0.0d0)
                    else
                        Phi_hat(i, j, k) = Q_hat(i, j, k) / LHS(i, j, k)
                    endif
                enddo
            enddo
        enddo

        ! 右辺を逆フーリエ変換
        call dfftw_plan_dft_c2r_3d(plan, NX, NY, NZ, Phi_hat, Phi, FFTW_ESTIMATE)
        call dfftw_execute(plan, Phi_hat, Phi)
        call dfftw_destroy_plan(plan)

    end subroutine fft_solve

    subroutine fft_navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, &
                      Bx, By, Bz, Tx, Ty, Tz, Tx0, Ty0, Tz0, Fx, Fy, Fz, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Tx(1:NX, 1:NY, 1:NZ), Ty(1:NX, 1:NY, 1:NZ), Tz(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Tx0(1:NX, 1:NY, 1:NZ), Ty0(1:NX, 1:NY, 1:NZ), Tz0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
        integer, intent(in) :: step
        integer i, j, k
        real(8) Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :, :) = Ax(:, :, :)
            Ay0(:, :, :) = Ay(:, :, :)
            Az0(:, :, :) = Az(:, :, :)
            Tx0(:, :, :) = Tx(:, :, :)
            Ty0(:, :, :) = Ty(:, :, :)
            Tz0(:, :, :) = Tz(:, :, :)
        endif

        ! 左辺の定数部分の計算
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    LHS(i, j, k) = 1.0d0 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                                           + 2*(1-cos((j-1)*dY))/dY**2 &
                                                           + 2*(1-cos((k-1)*dZ))/dZ**2)
                enddo
            enddo
        enddo

        ! 速度場Upを予測
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = U(i, j, k) &
                                -dt*(-P(i, j, k)+P(i+1, j, k))/dX &
                                +dt*(3.0d0*Ax(i, j, k) - Ax0(i, j, k))/2.0d0 &
                                +dt*Bx(i, j, k)/2.0d0 &
                                +dt*(Tx(i, j, k) + Tx0(i, j, k))/2.0d0 &
                                +dt*Fx(i, j, k)
                enddo
            enddo
        enddo

        call fft_solve(Q, LHS, Up(1:NX, 1:NY, 1:NZ))
        

        ! 同様にVpについても求める
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = V(i, j, k) &
                                -dt*(-P(i, j, k)+P(i, j+1, k))/dY &
                                +dt*(3.0d0*Ay(i, j, k) - Ay0(i, j, k))/2.0d0 &
                                +dt*By(i, j, k)/2.0d0 &
                                +dt*(Ty(i, j, k) + Ty0(i, j, k))/2.0d0 &
                                +dt*Fy(i, j, k)
                enddo
            enddo
        enddo

        call fft_solve(Q, LHS, Vp(1:NX, 1:NY, 1:NZ))
        

        ! 同様にWpについても求める
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = W(i, j, k) &
                                -dt*(-P(i, j, k)+P(i, j, k+1))/dZ &
                                +dt*(3.0d0*Az(i, j, k) - Az0(i, j, k))/2.0d0 &
                                +dt*Bz(i, j, k)/2.0d0 &
                                +dt*(Tz(i, j, k) + Tz0(i, j, k))/2.0d0 &
                                +dt*Fz(i, j, k)
                enddo
            enddo
        enddo

        call fft_solve(Q, LHS, Wp(1:NX, 1:NY, 1:NZ))
        
        ! Phiを求めるときに0番目の値も必要
        call PBM(Up)
        call PBM(Vp)
        call PBM(Wp)

        ! n-1ステップ目の保存
        Ax0(:, :, :) = Ax(:, :, :)
        Ay0(:, :, :) = Ay(:, :, :)
        Az0(:, :, :) = Az(:, :, :)
        Tx0(:, :, :) = Tx(:, :, :)
        Ty0(:, :, :) = Ty(:, :, :)
        Tz0(:, :, :) = Tz(:, :, :)
    end subroutine fft_navier

    subroutine fft_poisson(Up, Vp, Wp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        real(8) Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)
        poisson_itr = -1

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    LHS(i, j, k) = -(2*(1-cos((i-1)*dX))/dX**2 + &
                                     2*(1-cos((j-1)*dY))/dY**2 + &
                                     2*(1-cos((k-1)*dZ))/dZ**2)
                enddo
            enddo
        enddo
        
        ! 右辺の計算
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = ((-Up(i-1, j, k) + Up(i, j, k)) / dX &
                                 +(-Vp(i, j-1, k) + Vp(i, j, k)) / dY &
                                 +(-Wp(i, j, k-1) + Wp(i, j, k)) / dZ) / dt
                enddo
            enddo
        enddo

        call fft_solve(Q, LHS, Phi(1:NX, 1:NY, 1:NZ))
        call PBM(Phi)
        
    end subroutine fft_poisson

    subroutine fft_march(Up, Vp, Wp, U, V, W, Phi, P)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    U(i, j, k) = Up(i, j, k) - dt*(-Phi(i, j, k)+Phi(i+1, j, k))/dX
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    V(i, j, k) = Vp(i, j, k) - dt*(-Phi(i, j, k)+Phi(i, j+1, k))/dY
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    W(i, j, k) = Wp(i, j, k) - dt*(-Phi(i, j, k)+Phi(i, j, k+1))/dZ
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    P(i, j, k) = P(i, j, k) + Phi(i, j, k) &
                                -dt*beta/Re/2.0d0*((Phi(i-1, j, k)-2*Phi(i, j, k)+Phi(i+1, j, k))/dX**2 &
                                                  +(Phi(i, j-1, k)-2*Phi(i, j, k)+Phi(i, j+1, k))/dY**2 &
                                                  +(Phi(i, j, k-1)-2*Phi(i, j, k)+Phi(i, j, k+1))/dZ**2)
                enddo
            enddo
        enddo

        call PBM(U)
        call PBM(V)
        call PBM(W)
        call PBM(P)
    end subroutine fft_march

    subroutine energy_single(U, V, W, step)  ! エネルギースペクトル
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        real(8) U_tmp(0:NX+1, 0:NY+1, 0:NZ+1), V_tmp(0:NX+1, 0:NY+1, 0:NZ+1), W_tmp(0:NX+1, 0:NY+1, 0:NZ+1)
        complex(8) U_hat(1:NX/2+1, 1:NY, 1:NZ), V_hat(1:NX/2+1, 1:NY, 1:NZ), W_hat(1:NX/2+1, 1:NY, 1:NZ)
        real(8) E_tmp(1:NX/2+1, 1:NY, 1:NZ)
        real(8) K_abs(1:NX/2+1, 1:NY, 1:NZ)
        real(8) Energy(0:NX)  ! 格子数程度あれば十分
        real(8) K_energy
        integer i, j, k, index
        integer(8) plan
        character(8) str
        real(8) kx, ky, kz
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        ! そのまま用いると，Energy(0)=0ではないが，そのまま求めたK_energyと一致する．
        U_tmp(:, :, :) = U(:, :, :) * U_C
        V_tmp(:, :, :) = V(:, :, :) * U_C
        W_tmp(:, :, :) = W(:, :, :) * U_C
        ! 変動成分にすると，Energy(0)=0になり，変動成分を使ったK_energyと一致する．
        ! U_tmp(:, :, :) = U(:, :, :) - sum(U(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
        ! V_tmp(:, :, :) = V(:, :, :) - sum(V(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
        ! W_tmp(:, :, :) = W(:, :, :) - sum(W(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
        ! U_tmp(:, :, :) = U_tmp(:, :, :) * U_C
        ! V_tmp(:, :, :) = V_tmp(:, :, :) * U_C
        ! W_tmp(:, :, :) = W_tmp(:, :, :) * U_C


        ! 速度場をフーリエ変換
        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, U_tmp(1:NX, 1:NY, 1:NZ), U_hat, FFTW_ESTIMATE)  ! 分けて書かないとバグった
        call dfftw_execute(plan, U_tmp(1:NX, 1:NY, 1:NZ), U_hat)
        call dfftw_destroy_plan(plan)

        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, V_tmp(1:NX, 1:NY, 1:NZ), V_hat, FFTW_ESTIMATE)
        call dfftw_execute(plan, V_tmp(1:NX, 1:NY, 1:NZ), V_hat)
        call dfftw_destroy_plan(plan)

        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, W_tmp(1:NX, 1:NY, 1:NZ), W_hat, FFTW_ESTIMATE)
        call dfftw_execute(plan, W_tmp(1:NX, 1:NY, 1:NZ), W_hat)
        call dfftw_destroy_plan(plan)

        ! 格子点数で割る
        U_hat(:, :, :) = U_hat(:, :, :)/dble(NX*NY*NZ)
        V_hat(:, :, :) = V_hat(:, :, :)/dble(NX*NY*NZ)
        W_hat(:, :, :) = W_hat(:, :, :)/dble(NX*NY*NZ)


        ! 配列要素に対するエネルギー
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX/2+1
                    E_tmp(i, j, k) = abs(U_hat(i, j, k))**2.0d0 + abs(V_hat(i, j, k))**2.0d0 + abs(W_hat(i, j, k))**2.0d0
                enddo
            enddo
        enddo

        ! 配列要素に対する波数
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX/2+1
                    ! 0, 1, 2, ..., 63, 64, -63, -62, ..., -3, -2, -1
                    kx = sign(1.0d0, dble(NX + 3)/2.0d0 - dble(i)) * (- dble(abs(NX/2 + 1 - i)) + dble(NX)/2.0d0)  ! 増田さん参考
                    ky = sign(1.0d0, dble(NY + 3)/2.0d0 - dble(j)) * (- dble(abs(NY/2 + 1 - j)) + dble(NY)/2.0d0)
                    kz = sign(1.0d0, dble(NZ + 3)/2.0d0 - dble(k)) * (- dble(abs(NZ/2 + 1 - k)) + dble(NZ)/2.0d0)
                    K_abs(i, j, k) = sqrt(kx**2.0d0 + ky**2.0d0 + kz**2.0d0)
                enddo
            enddo
        enddo

        ! 波数を四捨五入し、対応する整数の波数にエネルギーを足し合わせる
        Energy(:) = 0.0d0
        ! kx = 2 ~ NX/2 は，絶対値が同じ波数が2つあるので2倍する．
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX/2+1
                    index = nint(K_abs(i, j, k))
                    if (i==1 .or. i==NX/2+1) then
                        Energy(index) = Energy(index) + E_tmp(i, j, k)/2.d0
                    else
                        Energy(index) = Energy(index) + E_tmp(i, j, k)/2.d0*2.0d0
                    endif
                enddo
            enddo
        enddo

        ! 出力
        ! open(10, file=dir//str//'.d')
        ! do i = 0, NX
        !     write(10, '(I4, e12.4)') i, Energy(i)
        ! enddo
        ! close(10)


        ! 運動エネルギーとエネルギーカスケードの総和
        K_energy = sum(U_tmp(1:NX, 1:NY, 1:NZ)**2.0d0) + sum(V_tmp(1:NX, 1:NY, 1:NZ)**2.0d0) + sum(W_tmp(1:NX, 1:NY, 1:NZ)**2.0d0)
        K_energy = K_energy/2.d0
        K_energy = K_energy/(NX*NY*NZ)
        write(*, '(I6, 5e12.4)') step, K_energy, sum(Energy(:)), K_energy-sum(Energy(:)), Energy(0), K_energy/sum(Energy(:))   ! 一致するはず
        ! open(30, file = 'fft/debag_energy.d', position='append')
        ! write(30, '(I6, 2e12.4)') step, K_energy, sum(Energy(:))
        ! close(30)

    end subroutine energy_single
end module fft

module ibm
    use smac
    use fft
    implicit none
    ! IBM用パラメータ
    real(8), parameter :: DC = 2*PI/8.0d0  ! 円柱直径
    integer, parameter :: NC = 4  ! 円柱の個数
    integer, parameter :: NS = 1  ! 反復回数
    integer, parameter :: NL = nint(PI*DC/dX)  ! 円周方向の分割数
    real(8), parameter :: dV = PI*DC/NL * dX * dX
contains
    subroutine ibm_init(X, Y, Z, Xc, Yc, Zc, Uc, Vc, Wc)
        real(8), intent(out) :: X(0:NX+1, 0:NY+1, 0:NZ+1), Y(0:NX+1, 0:NY+1, 0:NZ+1), Z(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Xc(1:NC, 1:NL, 1:NZ), Yc(1:NC, 1:NL, 1:NZ), Zc(1:NC, 1:NL, 1:NZ)
        real(8), intent(out) :: Uc(1:NC, 1:NL, 1:NZ), Vc(1:NC, 1:NL, 1:NZ), Wc(1:NC, 1:NL, 1:NZ)
        integer i, j, k

        if (dX /= dY .or. dX /= dZ) then
            stop 'dX, dY and dZ are not equal'
        endif

        write(*, '(a, I5)') 'NL   :', NL
        write(*, '(a, E12.4)') 'dX**3:', dX**3
        write(*, '(a, E12.4)') 'dV   :', dV

        ! 格子点中心座標
        do i = 0, NX+1
            X(i, :, :) = (i-0.5d0)*dX
        enddo
        do j = 0, NY+1
            Y(:, j, :) = (j-0.5d0)*dY
        enddo
        do k = 0, NZ+1
            Z(:, :, k) = (k-0.5d0)*dZ
        enddo

        ! 円柱上の座標
        if (NC == 4) then
            do j = 1, NL
                Xc(1, j, :) = 2*PI* 5/16 + DC/2 * cos(2*PI/NL*j)
                Yc(1, j, :) = 2*PI* 5/16 + DC/2 * sin(2*PI/NL*j)
                Xc(2, j, :) = 2*PI*11/16 + DC/2 * cos(2*PI/NL*j)
                Yc(2, j, :) = 2*PI* 5/16 + DC/2 * sin(2*PI/NL*j)
                Xc(3, j, :) = 2*PI* 5/16 + DC/2 * cos(2*PI/NL*j)
                Yc(3, j, :) = 2*PI*11/16 + DC/2 * sin(2*PI/NL*j)
                Xc(4, j, :) = 2*PI*11/16 + DC/2 * cos(2*PI/NL*j)
                Yc(4, j, :) = 2*PI*11/16 + DC/2 * sin(2*PI/NL*j)
            enddo
        else if (NC == 1) then
            do j = 1, NL
                Xc(1, j, :) = 2*PI/4 + DC/2 * cos(2*PI/NL*j)
                Yc(1, j, :) = 2*PI/2 + DC/2 * sin(2*PI/NL*j)
            enddo
        endif
        do k = 1, NZ
            Zc(:, :, k) = (k-0.25d0)*dZ
        enddo
        Xc(:, :, :) = Xc(:, :, :) + 10d-10*dX  ! 速度定義点から少しずらして外力点を定義
        Yc(:, :, :) = Yc(:, :, :) + 10d-10*dY
        Zc(:, :, :) = Zc(:, :, :) + 10d-10*dZ

        ! 円柱上の座標での速度
        Uc(:, :, :) = 0.0d0
        Vc(:, :, :) = 0.0d0
        Wc(:, :, :) = 0.0d0
        if (NC == 4) then
            do j = 1, NL
                Uc(1, j, :) = -sin(2*PI/NL*j)  ! 代表速度は出力時にかける
                Vc(1, j, :) = cos(2*PI/NL*j)
                Uc(2, j, :) = sin(2*PI/NL*j)
                Vc(2, j, :) = -cos(2*PI/NL*j)
                Uc(3, j, :) = sin(2*PI/NL*j)
                Vc(3, j, :) = -cos(2*PI/NL*j)
                Uc(4, j, :) = -sin(2*PI/NL*j)
                Vc(4, j, :) = cos(2*PI/NL*j)
            enddo
        endif
    end subroutine ibm_init

    subroutine ibm_get(Xc, Yc)
        real(8), intent(in) :: Xc(1:NC, 1:NL, 1:NZ), Yc(1:NC, 1:NL, 1:NZ)
        integer i, j, k
        character(8) str

        k = NZ/2
        do i = 1, NC
            write(str, '(I8.8)') i
            open(10, file=dir//str//'.d')
            do j = 1, NL
                write(10, '(2e12.4)') Xc(i, j, k)*L_C, Yc(i, j, k)*L_C
            enddo
            write(10, '(2e12.4)') Xc(i, 1, k)*L_C, Yc(i, 1, k)*L_C  ! 1周させる
            close(10)
        enddo
    end subroutine ibm_get

    subroutine ibm_preliminary(U, V, W, P, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Tx, Ty, Tz, Tx0, Ty0, Tz0, &
                               Ua, Va, Wa, Fxc, Fyc, Fzc, Up, Vp, Wp, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Tx(1:NX, 1:NY, 1:NZ), Ty(1:NX, 1:NY, 1:NZ), Tz(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Tx0(1:NX, 1:NY, 1:NZ), Ty0(1:NX, 1:NY, 1:NZ), Tz0(1:NX, 1:NY, 1:NZ)
        real(8), intent(out) :: Ua(0:NX+1, 0:NY+1, 0:NZ+1), Va(0:NX+1, 0:NY+1, 0:NZ+1), Wa(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: Fxc(1:NC, 1:NL, 1:NZ), Fyc(1:NC, 1:NL, 1:NZ), Fzc(1:NC, 1:NL, 1:NZ)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k
        
        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :, :) = Ax(:, :, :)
            Ay0(:, :, :) = Ay(:, :, :)
            Az0(:, :, :) = Az(:, :, :)
            Tx0(:, :, :) = Tx(:, :, :)
            Ty0(:, :, :) = Ty(:, :, :)
            Tz0(:, :, :) = Tz(:, :, :)
        endif

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Ua(i, j, k) = U(i, j, k) - dt*(-P(i, j, k)+P(i+1, j, k))/dX &
                                 +dt*(3.0d0*Ax(i, j, k) - Ax0(i, j, k))/2.0d0 &
                                 +dt*Bx(i, j, k) &
                                 +dt*(Tx(i, j, k) + Tx0(i, j, k))/2.0d0

                    Va(i, j, k) = V(i, j, k) - dt*(-P(i, j, k)+P(i, j+1, k))/dY &
                                 +dt*(3.0d0*Ay(i, j, k) - Ay0(i, j, k))/2.0d0 &
                                 +dt*By(i, j, k) &
                                 +dt*(Ty(i, j, k) + Ty0(i, j, k))/2.0d0

                    Wa(i, j, k) = W(i, j, k) - dt*(-P(i, j, k)+P(i, j, k+1))/dZ &
                                 +dt*(3.0d0*Az(i, j, k) - Az0(i, j, k))/2.0d0 &
                                 +dt*Bz(i, j, k) &
                                 +dt*(Tz(i, j, k) + Tz0(i, j, k))/2.0d0
                enddo
            enddo
        enddo
        ! 周期境界
        call PBM(Ua)
        call PBM(Va)
        call PBM(Wa)

        ! n-1ステップ目の保存
        Ax0(:, :, :) = Ax(:, :, :)
        Ay0(:, :, :) = Ay(:, :, :)
        Az0(:, :, :) = Az(:, :, :)
        Tx0(:, :, :) = Tx(:, :, :)
        Ty0(:, :, :) = Ty(:, :, :)
        Tz0(:, :, :) = Tz(:, :, :)

        ! 円柱上での外力初期化
        Fxc(:, :, :) = 0.0d0
        Fyc(:, :, :) = 0.0d0
        Fzc(:, :, :) = 0.0d0

        ! 反復回数によらずHelmholtzの式を使うため
        Up(:, :, :) = Ua(:, :, :)
        Vp(:, :, :) = Va(:, :, :)
        Wp(:, :, :) = Wa(:, :, :)

    end subroutine ibm_preliminary

    function delta(x, y, z, xc, yc, zc) result(out)
        real(8), intent(in) :: x, y, z, xc, yc, zc
        real(8) out
        out = weight((x-xc)/dX) * weight((y-yc)/dY) * weight((z-zc)/dZ)/(dX*dY*dZ)
    end function delta

    function weight(r) result(y)
        real(8), intent(in) :: r
        real(8) y
        if (abs(r)<1.5d0 .and. abs(r)>=0.5d0) then
            y = (5.0d0 - 3.0d0*abs(r) - sqrt(-3.0d0*(1-abs(r))**2 + 1.0d0))/6.0d0
        else if (abs(r)<0.5d0) then
            y = (1.0d0 + sqrt(-3.0d0 * r ** 2 + 1.0d0))/3.0d0
        else
            y = 0.0d0
        endif
    end function weight

    subroutine ibm_Helmholtz(Up, Vp, Wp, X, Y, Z, Xc, Yc, Zc, Uc, Vc, Wc, Fxc, Fyc, Fzc, fxint, fyint, fzint)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: X(0:NX+1, 0:NY+1, 0:NZ+1), Y(0:NX+1, 0:NY+1, 0:NZ+1), Z(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Xc(1:NC, 1:NL, 1:NZ), Yc(1:NC, 1:NL, 1:NZ), Zc(1:NC, 1:NL, 1:NZ)
        real(8), intent(in) :: Uc(1:NC, 1:NL, 1:NZ), Vc(1:NC, 1:NL, 1:NZ), Wc(1:NC, 1:NL, 1:NZ)
        real(8), intent(inout) :: Fxc(1:NC, 1:NL, 1:NZ), Fyc(1:NC, 1:NL, 1:NZ), Fzc(1:NC, 1:NL, 1:NZ)
        real(8), intent(out) :: fxint(1:NX, 1:NY, 1:NZ), fyint(1:NX, 1:NY, 1:NZ), fzint(1:NX, 1:NY, 1:NZ)
        integer i, j, k, l, m, n
        real(8) Ub(1:NC, 1:NL, 1:NZ), Vb(1:NC, 1:NL, 1:NZ), Wb(1:NC, 1:NL, 1:NZ)
        real(8) fxtmp(1:NX, 1:NY, 0:NZ+1), fytmp(1:NX, 1:NY, 0:NZ+1), fztmp(1:NX, 1:NY, 0:NZ+1)
        real(8) er
        real(8) count3d(0:NX+1, 0:NY+1, 0:NZ+1), count3d_debug(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        count3d(:, :, :) = 0.0d0
        count3d_debug(:, :, :, :) = 0.0d0
        
        Ub(:, :, :) = 0.0d0
        Vb(:, :, :) = 0.0d0
        Wb(:, :, :) = 0.0d0

        do n = 1, NZ
            do m = 1, NL
                do l = 1, NC
                    ! 外力点周囲の3*3*3点の和をとる
                    do k = n-1, n+1
                        do j = int(Yc(l, m, n)/dY), int(Yc(l, m, n)/dY) + 2
                            do i = int(Xc(l, m, n)/dX - 0.5d0), int(Xc(l, m, n)/dX - 0.5d0) + 2
                                Ub(l, m, n) = Ub(l, m, n) &
                                            + Up(i, j, k) * delta(X(i, j, k)+0.5d0*dX, Y(i, j, k), Z(i, j, k), &
                                                                  Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc(l, m, n)/dY - 0.5d0), int(Yc(l, m, n)/dY - 0.5d0) + 2
                            do i = int(Xc(l, m, n)/dX), int(Xc(l, m, n)/dX) + 2
                                Vb(l, m, n) = Vb(l, m, n) &
                                            + Vp(i, j, k) * delta(X(i, j, k), Y(i, j, k)+0.5d0*dY, Z(i, j, k), &
                                                                  Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc(l, m, n)/dY) , int(Yc(l, m, n)/dY) + 2
                            do i = int(Xc(l, m, n)/dX), int(Xc(l, m, n)/dX) + 2
                                Wb(l, m, n) = Wb(l, m, n) &
                                            + Wp(i, j, k) * delta(X(i, j, k), Y(i, j, k), Z(i, j, k)+0.5d0*dZ, &
                                                                  Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do n = 1, NZ
            do m = 1, NL
                do l = 1, NC
                    Fxc(l, m, n) = Fxc(l, m, n) + (Uc(l, m, n) - Ub(l, m, n))/dt
                    Fyc(l, m, n) = Fyc(l, m, n) + (Vc(l, m, n) - Vb(l, m, n))/dt
                    Fzc(l, m, n) = Fzc(l, m, n) + (Wc(l, m, n) - Wb(l, m, n))/dt
                enddo
            enddo
        enddo

        fxtmp(:, :, :) = 0.0d0
        fytmp(:, :, :) = 0.0d0
        fztmp(:, :, :) = 0.0d0
        do n = 1, NZ
            do m = 1, NL
                do l = 1, NC
                    ! 格子点周辺の外力点の和をとる、つまり外力点周辺の格子点へ和をのこす
                    do k = n-1, n+1
                        do j = int(Yc(l, m, n)/dY), int(Yc(l, m, n)/dY) + 2
                            do i = int(Xc(l, m, n)/dX - 0.5d0), int(Xc(l, m, n)/dX - 0.5d0) + 2
                                fxtmp(i, j, k) = fxtmp(i, j, k) &
                                               + Fxc(l, m, n) * delta(X(i, j, k)+0.5d0*dX, Y(i, j, k), Z(i, j, k), &
                                                                      Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dV
                                count3d(i, j, k) = count3d(i, j, k) + 1
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc(l, m, n)/dY - 0.5d0), int(Yc(l, m, n)/dY - 0.5d0) + 2
                            do i = int(Xc(l, m, n)/dX), int(Xc(l, m, n)/dX) + 2
                                fytmp(i, j, k) = fytmp(i, j, k) &
                                               + Fyc(l, m, n) * delta(X(i, j, k), Y(i, j, k)+0.5d0*dY, Z(i, j, k), &
                                                                      Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dV
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc(l, m, n)/dY) , int(Yc(l, m, n)/dY) + 2
                            do i = int(Xc(l, m, n)/dX), int(Xc(l, m, n)/dX) + 2
                                fztmp(i, j, k) = fztmp(i, j, k) &
                                               + Fzc(l, m, n) * delta(X(i, j, k), Y(i, j, k), Z(i, j, k)+0.5d0*dZ, &
                                                                      Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dV
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! do k = 0, NZ+1
        !     do j = 0, NY+1
        !         write(*, '(100I3)') count3d(:, j, k)
        !     enddo
        ! enddo
        ! do j = 0, NY+1
        !     do i = 0, NX+1
        !         write(*, '(100I6)', advance='no') count3d(i, j, 1)
        !     enddo
        !     write(*, *) ''
        ! enddo
        ! write(*, *) ''
        ! do k = 0, NZ+1
        !     write(*, '(100I6)', advance='no') sum(count3d(:, :, k))
        ! enddo
        ! write(*, *) ''

        count3d_debug(1, :, :, :) = count3d(:, :, :)

        call get_data_xy(Up, Vp, Wp, count3d_debug, 0)

        fxint(:, :, :) = fxtmp(:, :, 1:NZ)
        fyint(:, :, :) = fytmp(:, :, 1:NZ)
        fzint(:, :, :) = fztmp(:, :, 1:NZ)
        fxint(:, :, 1) = fxint(:, :, 1) + fxtmp(:, :, NZ+1)
        fyint(:, :, 1) = fyint(:, :, 1) + fytmp(:, :, NZ+1)
        fzint(:, :, 1) = fzint(:, :, 1) + fztmp(:, :, NZ+1)
        fxint(:, :, NZ) = fxint(:, :, NZ) + fxtmp(:, :, 0)
        fyint(:, :, NZ) = fyint(:, :, NZ) + fytmp(:, :, 0)
        fzint(:, :, NZ) = fzint(:, :, NZ) + fztmp(:, :, 0)

    end subroutine ibm_Helmholtz

    subroutine ibm_predict(Ua, Va, Wa, fxint, fyint, fzint, Bx, By, Bz, Up, Vp, Wp)
        real(8), intent(in) :: Ua(0:NX+1, 0:NY+1, 0:NZ+1), Va(0:NX+1, 0:NY+1, 0:NZ+1), Wa(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: fxint(1:NX, 1:NY, 1:NZ), fyint(1:NX, 1:NY, 1:NZ), fzint(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        real(8) Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)

        ! 左辺の定数部分の計算
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    LHS(i, j, k) = 1 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                                       + 2*(1-cos((j-1)*dY))/dY**2 &
                                                       + 2*(1-cos((k-1)*dZ))/dZ**2)
                enddo
            enddo
        enddo

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = Ua(i, j, k) + dt*fxint(i, j, k) - dt*Bx(i, j, k)/2
                enddo
            enddo
        enddo
        call fft_solve(Q, LHS, Up(1:NX, 1:NY, 1:NZ))

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = Va(i, j, k) + dt*fyint(i, j, k) - dt*By(i, j, k)/2
                enddo
            enddo
        enddo
        call fft_solve(Q, LHS, Vp(1:NX, 1:NY, 1:NZ))

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = Wa(i, j, k) + dt*fzint(i, j, k) - dt*Bz(i, j, k)/2
                enddo
            enddo
        enddo
        call fft_solve(Q, LHS, Wp(1:NX, 1:NY, 1:NZ))

        ! Phiを求めるときに0番目の値も必要
        call PBM(Up)
        call PBM(Vp)
        call PBM(Wp)
    end subroutine ibm_predict
end module ibm


program main
    use smac
    use fft
    use ibm
    implicit none
    real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) Cx(6, NX, NY, NZ)
    real(8) Cpx(6, 0:NX, NY, NZ), Cnx(6, 0:NX, NY, NZ)
    real(8) Cpy(6, NX, 0:NY, NZ), Cny(6, NX, 0:NY, NZ)
    real(8) Cpz(6, NX, NY, 0:NZ), Cnz(6, NX, NY, 0:NZ)
    real(8) Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
    real(8) Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
    real(8) Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
    real(8) Bx0(1:NX, 1:NY, 1:NZ), By0(1:NX, 1:NY, 1:NZ), Bz0(1:NX, 1:NY, 1:NZ)
    real(8) Tx(1:NX, 1:NY, 1:NZ), Ty(1:NX, 1:NY, 1:NZ), Tz(1:NX, 1:NY, 1:NZ)
    real(8) Tx0(1:NX, 1:NY, 1:NZ), Ty0(1:NX, 1:NY, 1:NZ), Tz0(1:NX, 1:NY, 1:NZ)
    real(8) Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) P(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) Phi(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
    ! ibmのために追加した変数
    real(8) X(0:NX+1, 0:NY+1, 0:NZ+1), Y(0:NX+1, 0:NY+1, 0:NZ+1), Z(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) Xc(1:NC, 1:NL, 1:NZ), Yc(1:NC, 1:NL, 1:NZ), Zc(1:NC, 1:NL, 1:NZ)
    real(8) Uc(1:NC, 1:NL, 1:NZ), Vc(1:NC, 1:NL, 1:NZ), Wc(1:NC, 1:NL, 1:NZ)
    real(8) Ua(0:NX+1, 0:NY+1, 0:NZ+1), Va(0:NX+1, 0:NY+1, 0:NZ+1), Wa(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) Fxc(1:NC, 1:NL, 1:NZ), Fyc(1:NC, 1:NL, 1:NZ), Fzc(1:NC, 1:NL, 1:NZ)
    real(8) fxint(1:NX, 1:NY, 1:NZ), fyint(1:NX, 1:NY, 1:NZ), fzint(1:NX, 1:NY, 1:NZ)
    integer s
    integer step
    real(8) t1, t2, t3, t4, t5, t12, t23, t34, t45
    real(8) total_time, t_start, t_end, others_time
    ! real(8) Ep(1:NX, 1:NY, 1:NZ)
    ! integer i, j, k
    call mk_dir(dir)
    call cpu_time(t_start)
    t12 = 0.0d0
    t23 = 0.0d0
    t34 = 0.0d0
    t45 = 0.0d0

    call init(U, V, W, P, Phi, C, Fx, Fy, Fz)
    if (input_step > 0) call input(U, V, W, P, C)
    if (method == 2) call ibm_init(X, Y, Z, Xc, Yc, Zc, Uc, Vc, Wc)
    if (method == 2) call ibm_get(Xc, Yc) 
    ! call get_data(U, V, W, 0)
    ! call get_data_xy(U, V, W, C, 0)
    ! call get_data_xz(U, V, W, C, 50000)
    ! do k = 1, NZ
    !     do j = 1, NY
    !         do i = 1, NX
    !             Ep(i, j, k) = 0.5d0*(1.0d0-beta)/Re/Wi*(Lp**2.0d0-3.0d0)*log(f(C(:, i, j, k)))
    !         enddo
    !     enddo
    ! enddo
    ! write(*, *) sum(Ep)/(NX*NY*NZ)
    
    do step = 1, Nstep
        eigen_method = 0
        ! if (mod(step, 100)==0) eigen_method = 0  ! 0:カルダノ、1:シルベスター
        call cpu_time(t1)
        call CpxCnx(C, Cpx, Cnx)
        call CpyCny(C, Cpy, Cny)
        call CpzCnz(C, Cpz, Cnz)

        call cpu_time(t2)
        call Cstar(Cpx, Cnx, Cpy, Cny, Cpz, Cnz, U, V, W, Cx)
        call Lyapunov(Cx, U, V, W, C)
        call cpu_time(t3)

        call polymer_stress(C, Tx, Ty, Tz)
        call convection(U, V, W, Ax, Ay, Az)
        call viscous(U, V, W, Bx, By, Bz)

        if (method == 0) then
            call navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, &
                        Bx, By, Bz, Bx0, By0, Bz0, Tx, Ty, Tz, Tx0, Ty0, Tz0, Fx, Fy, Fz, step)
            call cpu_time(t4)
            call poisson(Up, Vp, Wp, Phi)
            call march(Up, Vp, Wp, U, V, W, Phi, P)
        endif

        if (method == 1) then
            call fft_navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, &
                            Bx, By, Bz, Tx, Ty, Tz, Tx0, Ty0, Tz0, Fx, Fy, Fz, step)
            call cpu_time(t4)
            call fft_poisson(Up, Vp, Wp, Phi)
            call fft_march(Up, Vp, Wp, U, V, W, Phi, P)
        endif

        if (method == 2) then
            call ibm_preliminary(U, V, W, P, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Tx, Ty, Tz, Tx0, Ty0, Tz0, &
                                 Ua, Va, Wa, Fxc, Fyc, Fzc, Up, Vp, Wp, step)
            do s = 1, NS
                call ibm_Helmholtz(Up, Vp, Wp, X, Y, Z, Xc, Yc, Zc, Uc, Vc, Wc, Fxc, Fyc, Fzc, fxint, fyint, fzint)
                call ibm_predict(Ua, Va, Wa, fxint, fyint, fzint, Bx, By, Bz, Up, Vp, Wp)
            enddo
            call cpu_time(t4)
            call fft_poisson(Up, Vp, Wp, Phi)
            call fft_march(Up, Vp, Wp, U, V, W, Phi, P)
        endif

        ! if (step > 50000) then
        !     Fx(:, :, :) = 0.0d0
        !     Fy(:, :, :) = 0.0d0
        !     Fz(:, :, :) = 0.0d0
        !     if (mod(step, 50) == 0) call get_data(U, V, W, C, step)
        ! endif

        call cpu_time(t5)
        
        call logging(U, V, W, C, step, t_start)
        ! if (mod(step, Gstep)==0) call get_data(U, V, W, C, step)
        ! if (mod(step, Gstep)==0) call get_data_xz(U, V, W, C, step)
        if (mod(step, Gstep)==0) call taylor_debag(U, V, W, step)
        ! if (mod(step, Estep)==0) call energy_single(U, V, W, step)
        if (mod(step, output_step)==0) call output(U, V, W, P, C, step)
        if (sum(U**2)*0.0d0 /= 0.0d0) stop 'NaN value'  ! NaNの判定
        t12 = t12 + t2-t1
        t23 = t23 + t3-t2
        t34 = t34 + t4-t3
        t45 = t45 + t5-t4
        ! write(*, '(12e12.4)') U(1:12, 1, 1)
    enddo

    call cpu_time(t_end)
    total_time = t_end-t_start
    others_time = total_time - (t12+t23+t34+t45)
    write(*, '(a9, F10.3, a3)') 'Total   :', total_time, '[s]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'C+-     :', t12, '[s]', t12/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Newton  :', t23, '[s]', t23/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Navier  :', t34, '[s]', t34/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Poisson :', t45, '[s]', t45/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Others  :', others_time, '[s]', others_time/total_time*100, '[%]'

    ! open(10, file = 'debag3.d', position='append')
    ! write(10, *) NX, total_time
    ! close(10)
    ! write(*, *) NX, total_time
    
end program main


! 名大スパコン
! module load fftw
! frtpx fenep.f90 -lfftw3
! pjsub single.sh