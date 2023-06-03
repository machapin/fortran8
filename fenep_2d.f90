module smac
    implicit none
    ! gfortran fenep_2d.f90 -I$HOME/local/include -L$HOME/local/lib -lfftw3 && ./a.out
    ! ステップ数
    integer, parameter :: Nstep = 1
    integer, parameter :: Gstep = 10  ! データを取得する間隔
    integer, parameter :: Estep = 1000000  ! エネルギースペクトルを取得する間隔
    integer, parameter :: Dstep = 200  ! デバッグする間隔
    integer, parameter :: input_step = 0  ! 0以外で初期条件をファイルから読み込む
    integer, parameter :: output_step = 1000000  ! 配列を保存する間隔
    real(8), parameter :: PI = acos(-1.0d0)
    ! パラメータ
    integer, parameter :: NX = 128, NY = NX
    real(8), parameter :: dX = 2*PI/NX, dY = 2*PI/NY  ! 規格化長さは2*PI
    real(8), parameter :: dt = 0.005d0
    ! 手法
    integer, parameter :: method = 2  ! 0:陽解法、1:FFT、2:IBM、3:固体壁
    integer, parameter :: ibm_type = 21  ! 1:壁表面、2:壁の中全て、3:壁の中半分、11:円柱表面1つ、12:円柱表面4つ、21:円柱内部1つ、22:円柱内部4つ
    real(8), parameter :: DC = 2*PI / 4.0d0  ! 円柱の直径、円柱の中心点を確認すること！
    integer, parameter :: flow_type = 0  ! 0:外力なし(f=0)、1:クエット流れ(Uwall=1)、2:ポアズイユ流れ(fx=1)、3:テイラーグリーン外力、4:テイラーグリーン渦の減衰
    integer, parameter :: eigen_method = 0  ! 0:固有値計算、1:SPDならindex3
    ! 無次元パラメータ
    ! real(8), parameter :: Re = 1.0d0
    real(8), parameter :: beta = 1.0d0
    real(8), parameter :: Wi = 1.0d0
    real(8), parameter :: Lp = 55.0d0
    ! 有次元パラメータ
    ! real(8), parameter :: L_C = 1.0d0  ! 長さが100なら100/2*PI
    ! real(8), parameter :: U_C = 1.0d0  ! 本来は乱流テイラーグリーン渦の平均流の速さ
    ! real(8), parameter :: nu = L_C*U_C/Re
    ! 実験と同様のパラメータ
    integer, parameter :: f_C = 2  ! 円柱回転速度[rps]
    real(8), parameter :: D_C = 0.03d0  ! 円柱直径[m]
    real(8), parameter :: U_C = f_C * PI * D_C  ! 代表速度
    real(8), parameter :: L_C = D_C / DC  ! 代表長さ
    real(8), parameter :: nu = 1.0d-6  ! 動粘性係数
    real(8), parameter :: Re = U_C*L_C/nu
    ! その他のパラメータ
    real(8), parameter :: dX_C = dX*L_C, dY_C = dY*L_C
    real(8), parameter :: dt_C = dt*L_C/U_C
    ! グローバル変数
    character(64) :: dir = './tmp/'
    real(8) counter(0:3)

contains
    subroutine init(U, V, P, Phi, C, Fx, Fy)  ! はじめに実行
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(out) :: P(0:NX+1, 0:NY+1), Phi(0:NX+1, 0:NY+1)
        real(8), intent(out) :: C(3, 0:NX+1, 0:NY+1)
        real(8), intent(out) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        integer i, j
        character(2) str
        write(str, '(I2.2)') ibm_type

        if (method == 1) dir = trim(dir)//'fft'
        if (method == 2) dir = trim(dir)//'ibm'//str
        if (method == 3) dir = trim(dir)//'kabe'
        if (flow_type == 1) dir = trim(dir)//'_couette'
        if (flow_type == 2) dir = trim(dir)//'_poiseuille'
        if (flow_type == 3) dir = trim(dir)//'_taylor'
        if (flow_type == 4) dir = trim(dir)//'_taylor_decay'
        if (beta == 1.0d0) dir = trim(dir)//'_newton'
        dir = trim(dir)//'/'
        write(*, '(a, a)') 'dir:  ', dir
        call mk_dir(dir)

        ! write(*, '(a, F8.3, F8.3, a8)') 'Von Neumann:', 1.0d0/Re*dt/dX**2, 1.0d0/Re*dt/dY**2, '< 0.167'
        write(*, '(a, F8.3, F8.3, a8)') 'CFL:', 1.0d0*dt/dX*6.0d0, 1.0d0*dt/dY*6.0d0,'< 1.0'
        write(*, '(a, F8.3)') 'L_C =', L_C
        write(*, '(a, F8.3)') 'U_C =', U_C
        write(*, '(a, E12.4)') 'nu =', nu
        write(*, '(a, F10.3)') 'Re  =', Re
        write(*, '(a, E12.4)') 'dX_C =', dX_C
        write(*, '(a, E12.4)') 'dt_C =', dt_C
        write(*, '(a, F8.3)') 'step /[s]=', 1.0d0/dt_C

        ! 初期条件（Taylor-Green）
        ! U(:, :) = 0.0d0
        ! V(:, :) = 0.0d0
        call random_number(U)
        call random_number(V)
        U(:, :) = 0.001d0 * (U(:, :)-0.5d0)
        V(:, :) = 0.001d0 * (V(:, :)-0.5d0)
        P(:, :) = 0.0d0
        Phi(:, :) = 0.0d0
        Fx(:, :) = 0.0d0
        Fy(:, :) = 0.0d0
        if (flow_type == 2) then  ! ポアズイユ流れ
            Fx(:, :) = 1.0d0
        endif
        if (flow_type == 3) then
            do j = 1, NY
                do i = 1, NX
                    Fx(i, j) = -sin(i*dX) * cos((j-0.5d0)*dY)  ! 増田さん, 安房井さん
                    Fy(i, j) = cos((i-0.5d0)*dX) * sin(j*dY)
                    ! Fx(i, j) = -sin(2.0d0*(j-0.5d0)*dY)  ! 小井手さん
                    ! Fy(i, j) = sin(2.0d0*(i-0.5d0)*dX)
                enddo
            enddo
        endif
        if (flow_type == 4) then
            do j = 1, NY
                do i = 1, NX
                    U(i, j) = -sin(i*dX) * cos((j-0.5d0)*dY)
                    V(i, j) = cos((i-0.5d0)*dX) * sin(j*dY)
                enddo
            enddo
        endif

        call random_number(C)
        C(:, :, :) = 0.001d0 * (C(:, :, :)-0.5d0)
        C(1, :, :) = C(1, :, :) + 1.0d0
        C(3, :, :) = C(3, :, :) + 1.0d0

        call PBM(U)
        call PBM(V)
        call PBM(P)
        call C_PBM(C)
     end subroutine init

    subroutine PBM(A)
        real(8), intent(inout) :: A(0:NX+1, 0:NY+1)
        A(0, :) = A(NX, :)
        A(NX+1, :) = A(1, :)
        A(:, 0) = A(:, NY)
        A(:, NY+1) = A(:, 1)
    end subroutine PBM

    subroutine C_PBM(C)
        real(8), intent(inout) :: C(3, 0:NX+1, 0:NY+1)
        C(:, 0, :) = C(:, NX, :)  ! 周期境界条件
        C(:, NX+1, :) = C(:, 1, :)
        C(:, :, 0) = C(:, :, NY)
        C(:, :, NY+1) = C(:, :, 1)
    end subroutine C_PBM


    subroutine CpxCnx(C, Cpx, Cnx)
        real(8), intent(in) :: C(3, 0:NX+1, 0:NY+1)
        real(8), intent(out) :: Cpx(3, 0:NX, NY), Cnx(3, 0:NX, NY)
        integer i, j, l, index
        real(8) Cptemp(3, 3), Cntemp(3, 3), Eigen(3, 4), mintemp
        counter(:) = 0.0d0

        do j = 1, NY
            do i = 1, NX
                Cptemp(:, 1) = C(:, i, j)*3/2 - C(:, i+1, j)/2
                Cntemp(:, 1) = C(:, i, j)/2 + C(:, i+1, j)/2
                Cptemp(:, 2) = C(:, i-1, j)/2 + C(:, i, j)/2
                Cntemp(:, 2) = -C(:, i-1, j)/2 + C(:, i, j)*3/2
                Cptemp(:, 3) = C(:, i-1, j)/4 + C(:, i, j) - C(:, i+1, j)/4
                Cntemp(:, 3) = -C(:, i-1, j)/4 + C(:, i, j) + C(:, i+1, j)/4

                do l = 1, 3
                    call Eigenvalues(Cptemp(:, l), Eigen(l, 1), Eigen(l, 2))
                    call Eigenvalues(Cntemp(:, l), Eigen(l, 3), Eigen(l, 4))
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
                    if (eigen_method == 1) index = 3  ! SPD条件を満たすならindexを3に固定
                    Cpx(:, i-1, j) = Cptemp(:, index)
                    Cnx(:, i, j) = Cntemp(:, index)
                else
                    Cpx(:, i-1, j) = C(:, i, j)
                    Cnx(:, i, j) = C(:, i, j)
                endif
                counter(index) = counter(index) + 1
            enddo
        enddo

        Cpx(:, NX, :) = Cpx(:, 0, :)
        Cnx(:, 0, :) = Cnx(:, NX, :)
    end subroutine CpxCnx

    subroutine CpyCny(C, Cpy, Cny)
        real(8), intent(in) :: C(3, 0:NX+1, 0:NY+1)
        real(8), intent(out) :: Cpy(3, NX, 0:NY), Cny(3, NX, 0:NY)
        integer i, j, l, index
        real(8) Cptemp(3, 3), Cntemp(3, 3), Eigen(3, 4), mintemp

        do j = 1, NY
            do i = 1, NX
                Cptemp(:, 1) = C(:, i, j)*3/2 - C(:, i, j+1)/2
                Cntemp(:, 1) = C(:, i, j)/2 + C(:, i, j+1)/2
                Cptemp(:, 2) = C(:, i, j-1)/2 + C(:, i, j)/2
                Cntemp(:, 2) = -C(:, i, j-1)/2 + C(:, i, j)*3/2
                Cptemp(:, 3) = C(:, i, j-1)/4 + C(:, i, j) - C(:, i, j+1)/4
                Cntemp(:, 3) = -C(:, i, j-1)/4 + C(:, i, j) + C(:, i, j+1)/4

                do l = 1, 3
                    call Eigenvalues(Cptemp(:, l), Eigen(l, 1), Eigen(l, 2))
                    call Eigenvalues(Cntemp(:, l), Eigen(l, 3), Eigen(l, 4))
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
                    if (eigen_method == 1) index = 3  ! SPD条件を満たすならindexを3に固定
                    Cpy(:, i, j-1) = Cptemp(:, index)
                    Cny(:, i, j) = Cntemp(:, index)
                else
                    Cpy(:, i, j-1) = C(:, i, j)
                    Cny(:, i, j) = C(:, i, j)
                endif
                counter(index) = counter(index) + 1
            enddo
        enddo

        Cpy(:, :, NY) = Cpy(:, :, 0)
        Cny(:, :, 0) = Cny(:, :, NY)
    end subroutine CpyCny

    subroutine Eigenvalues(a, re1, re2)
        real(8), intent(in) :: a(3)
        real(8), intent(out) :: re1, Re2
        real(8) D

        D = (a(1)-a(3))**2 + 4.0d0 * a(2)**2
        re1 = (a(1) + a(3) + sqrt(D))*0.5d0
        re2 = (a(1) + a(3) - sqrt(D))*0.5d0

    end subroutine Eigenvalues


    subroutine Cstar(Cpx, Cnx, Cpy, Cny, U, V, Cx)
        real(8), intent(in) :: Cpx(3, 0:NX, NY), Cnx(3, 0:NX, NY)
        real(8), intent(in) :: Cpy(3, NX, 0:NY), Cny(3, NX, 0:NY)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Cx(3, NX, NY)
        integer i, j, k
        real(8) s1, s2, s3, s4, s5, s6, s7, s8
        real(8) u1, u2, v1, v2

        do j = 1, NY
            do i = 1, NX
                u1 = U(i-1, j)
                u2 = U(i, j)
                v1 = V(i, j-1)
                v2 = V(i, j)

                s1 = (-u2 + abs(u2))*dt/(2*dX)
                s2 = 1.0d0/4 + (-u2 - abs(u2))*dt/(2*dX)
                s3 = 1.0d0/4 + (u1 - abs(u1))*dt/(2*dX)
                s4 = (u1 + abs(u1))*dt/(2*dX)

                s5 = (-v2 + abs(v2))*dt/(2*dY)
                s6 = 1.0d0/4 + (-v2 - abs(v2))*dt/(2*dY)
                s7 = 1.0d0/4 + (v1 - abs(v1))*dt/(2*dY)
                s8 = (v1 + abs(v1))*dt/(2*dY)

                Cx(:, i, j) = s3*Cpx(:, i-1, j) + s1*Cpx(:, i, j) + s4*Cnx(:, i-1, j) + s2*Cnx(:, i, j) &
                             +s7*Cpy(:, i, j-1) + s5*Cpy(:, i, j) + s8*Cny(:, i, j-1) + s6*Cny(:, i, j)
            enddo
        enddo

    end subroutine Cstar


    function f(C) result(y)
        real(8), intent(in) :: C(3) 
        real(8) y
        y = (Lp**2 - 2)/(Lp**2 - C(1)-C(3))  ! FENE-Pモデル
        ! y = 1.0d0  ! Oldroyd-Bモデル
    end function f


    subroutine Lyapunov(Cx, U, V, C)
        real(8), intent(in) :: Cx(3, NX, NY)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(out) :: C(3, 0:NX+1, 0:NY+1)
        integer i, j
        integer row, colum, itr
        real(8) A0(2, 2), b(3)
        real(8) x(3), x1(3), x2(3), l(3), l1(3), l2(3), Jac(3, 3), r(3)
        real(8), parameter :: dc = 1.0d-10  ! 3, 4, 5あたりがいい。
        integer, parameter :: itrmax = 10000
        real(8), parameter :: eps = 1.0d-10

        do j = 1, NY
            do i = 1, NX
                ! 定数部分を計算
                b(:) = Cx(:, i, j)
                b(1) = b(1) + dt/Wi
                b(3) = b(3) + dt/Wi

                A0(1, 1) = -dt*(-U(i-1, j) + U(i, j))/dX
                A0(1, 2) = -dt*(-V(i-1, j-1) + V(i+1, j-1) - V(i-1, j) + V(i+1, j))/(4*dX)
                A0(2, 1) = -dt*(-U(i-1, j-1) - U(i, j-1) + U(i-1, j+1) + U(i, j+1))/(4*dY)
                A0(2, 2) = -dt*(-V(i, j-1) + V(i, j))/dY
                A0 = transpose(A0)  ! 転置して関数に渡す。転置しないとC11とC22の分布が入れ替わり、全体的に大きなCがでる
                ! write(*, '(9e12.4)') A0(:, :)

                ! ニュートン法
                x(:) = C(:, i, j)  ! Cから探索開始
                ! x(:) = (/0.0, 0.0, 0.0/)
                do itr = 1, itrmax
                    ! ヤコビアンの計算
                    do colum = 1, 3
                        x1(:) = x(:)
                        x2(:) = x(:)
                        x1(colum) = x1(colum) + dc
                        x2(colum) = x2(colum) - dc
                        call Lyapunov_func(A0, b, x1, l1)
                        call Lyapunov_func(A0, b, x2, l2)
                        do row = 1, 3
                            Jac(row, colum) = (l1(row) - l2(row))/(2.0d0*dc)
                        enddo
                    enddo

                    call Lyapunov_func(A0, b, x, l)  ! f(xi)を計算
                    ! write(*, '(6e12.4)') l(:)
                    
                    call gauss(Jac, l, r)  ! ガウスの消去法
                    x(:) = x(:) - r(:)
                    ! write(*, *) itr, sum(r**2), sum(x**2), sqrt(sum(r**2)/sum(x**2))

                    if (sqrt(sum(r**2)/sum(x**2)) < eps) exit
                enddo
                C(:, i, j) = x(:)
            enddo
        enddo
        call C_PBM(C)
    end subroutine Lyapunov


    subroutine Lyapunov_func(A0, b, x, l)  ! l(x)=Mx-bを計算
        real(8), intent(in) :: A0(2, 2), b(3), x(3)
        real(8), intent(out) :: l(3)
        real(8) A(2, 2)

        A(:, :) = A0(:, :)
        A(1, 1) = A(1, 1) + (1.0d0 + f(x)*dt/Wi)/2.0d0
        A(2, 2) = A(2, 2) + (1.0d0 + f(x)*dt/Wi)/2.0d0

        l(1) = (A(1, 1) + A(1, 1))*x(1) + A(1, 2)*x(2) + A(1, 2)*x(2) - b(1)
        l(2) = A(2, 1)*x(1) + (A(2, 2) + A(1, 1))*x(2) + A(1, 2)*x(3) - b(2)
        l(3) = A(2, 1)*x(2) + A(2, 1)*x(2) + (A(2, 2) + A(2, 2))*x(3) - b(3)

    end subroutine Lyapunov_func


    subroutine gauss(a, b, x)  ! ガウスの消去法（ピボット選択あり）
        integer, parameter :: n = 3
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


    subroutine polymer_stress(C, Tx, Ty)
        real(8), intent(in) :: C(3, 0:NX+1, 0:NY+1)
        real(8), intent(out) :: Tx(1:NX, 1:NY), Ty(1:NX, 1:NY)
        integer i, j
        real(8) C0(3, 0:NX+1, 0:NY+1), E(3)

        E(:) = (/1.0, 0.0, 1.0/)

        do j = 1, NY
            do i = 1, NX
                C0(:, i, j) = (f(C(:, i, j))*C(:, i, j) - E(:))
            enddo
        enddo
        call C_PBM(C0)

        do j = 1, NY
            do i = 1, NX
                Tx(i, j) = (-C0(1, i, j) + C0(1, i+1, j))/dX &
                          +(-C0(2, i, j-1) - C0(2, i+1, j-1) + C0(2, i, j+1) + C0(2, i+1, j+1))/(4*dY)
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                Ty(i, j) = (-C0(2, i-1, j) + C0(2, i+1, j) - C0(2, i-1, j+1) + C0(2, i+1, j+1))/(4*dX) &
                          +(-C0(3, i, j) + C0(3, i, j+1))/dY
            enddo
        enddo

        Tx(:, :) = (1.0d0-beta)/(Re*Wi)*Tx(:, :)
        Ty(:, :) = (1.0d0-beta)/(Re*Wi)*Ty(:, :)
    end subroutine polymer_stress


    subroutine convection(U, V, Ax, Ay)  ! 発散型
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
        integer i, j

        do j = 1, NY
            do i = 1, NX
                Ax(i, j) = (-((U(i-1, j)+U(i, j))/2)**2 + ((U(i, j)+U(i+1, j))/2)**2)/dX &
                          +(-(V(i, j-1)+V(i+1, j-1))/2*(U(i, j-1)+U(i, j))/2 + (V(i, j)+V(i+1, j))/2*(U(i, j)+U(i, j+1))/2)/dY
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                Ay(i, j) = (-(U(i-1, j)+U(i-1, j+1))/2*(V(i-1, j)+V(i, j))/2 + (U(i, j)+U(i, j+1))/2*(V(i, j)+V(i+1, j))/2)/dX &
                          +(-((V(i, j-1)+V(i, j))/2)**2 + ((V(i, j)+V(i, j+1))/2)**2)/dY
            enddo
        enddo

        Ax(:, :) = -Ax(:, :)
        Ay(:, :) = -Ay(:, :)
    end subroutine convection

    
    subroutine viscous(U, V, Bx, By)  ! ニュー一定
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        integer i, j

        do j = 1, NY
            do i = 1, NX
                Bx(i, j) = (U(i-1, j) - 2*U(i, j) + U(i+1, j)) / dX**2 &
                          +(U(i, j-1) - 2*U(i, j) + U(i, j+1)) / dY**2 
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                By(i, j) = (V(i-1, j) - 2*V(i, j) + V(i+1, j)) / dX**2 &
                          +(V(i, j-1) - 2*V(i, j) + V(i, j+1)) / dY**2 
            enddo
        enddo

        Bx(:, :) = beta/Re*Bx(:, :)
        By(:, :) = beta/Re*By(:, :)
    end subroutine viscous

    subroutine navier(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, &
                      Bx, By, Bx0, By0, Tx, Ty, Tx0, Ty0, Fx, Fy, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
        real(8), intent(in) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        real(8), intent(in) :: Tx(1:NX, 1:NY), Ty(1:NX, 1:NY)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
        real(8), intent(inout) :: Bx0(1:NX, 1:NY), By0(1:NX, 1:NY)
        real(8), intent(inout) :: Tx0(1:NX, 1:NY), Ty0(1:NX, 1:NY)
        real(8), intent(in) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        integer, intent(in) :: step
        integer i, j

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :) = Ax(:, :)
            Ay0(:, :) = Ay(:, :)
            Bx0(:, :) = Bx(:, :)
            By0(:, :) = By(:, :)
            Tx0(:, :) = Tx(:, :)
            Ty0(:, :) = Ty(:, :)
        endif

        ! NS方程式で速度の予測
        do j = 1, NY
            do i = 1, NX
                Up(i, j) = U(i, j) - dt*(-P(i, j)+P(i+1, j))/dX &
                                    +dt*(3.0d0*(Ax(i, j) + Bx(i, j)) - (Ax0(i, j) + Bx0(i, j)))/2.0d0 &
                                    +dt*(Tx(i, j) + Tx0(i, j))/2.0d0 &
                                    +dt*Fx(i, j)
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                Vp(i, j) = V(i, j) - dt*(-P(i, j)+P(i, j+1))/dY &
                                    +dt*(3.0d0*(Ay(i, j) + By(i, j)) - (Ay0(i, j) + By0(i, j)))/2.0d0 &
                                    +dt*(Ty(i, j) + Ty0(i, j))/2.0d0 &
                                    +dt*Fy(i, j)
            enddo
        enddo

        ! Phiを求めるときに0番目の値も必要
        call PBM(Up)
        call PBM(Vp)

        ! n-1ステップ目の保存
        Ax0(:, :) = Ax(:, :)
        Ay0(:, :) = Ay(:, :)
        Bx0(:, :) = Bx(:, :)
        By0(:, :) = By(:, :)
        Tx0(:, :) = Tx(:, :)
        Ty0(:, :) = Ty(:, :)
    end subroutine navier

    subroutine poisson(Up, Vp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Phi(0:NX+1, 0:NY+1)
        integer i, j, l, itr
        real(8) BXM, BXP, BYM, BYP, B0
        real(8) er, er0
        real(8) E(1:NX, 1:NY), Q(1:NX, 1:NY)
        real(8), parameter :: SOR = 2.0d0 / (1 + sin(PI*dX))
        real(8), parameter :: eps = 1.0d-5
        integer, parameter :: itrmax = 10000

        ! 右辺の計算
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = ((-Up(i-1, j) + Up(i, j)) / dX &
                          +(-Vp(i, j-1) + Vp(i, j)) / dY) / dt
            enddo
        enddo
        
        BXM = 1.0d0 / dX**2
        BXP = 1.0d0 / dX**2
        BYM = 1.0d0 / dY**2
        BYP = 1.0d0 / dY**2
        B0 = BYM + BXM + BXP + BYP
        Phi(:, :) = 0.0d0
        ! SOR法
        do itr = 1, itrmax
            do j = 1, NY
                do l = 1, 2  ! 奇数と偶数に分けて計算する
                    do i = l, NX, 2
                        E(i, j) = BYM*Phi(i, j-1) + BXM*Phi(i-1, j) - B0*Phi(i, j) &
                                 +BXP*Phi(i+1, j) + BYP*Phi(i, j+1) - Q(i, j)
                        Phi(i, j) = Phi(i, j) + SOR*E(i, j) / B0
                    enddo
                enddo
            enddo
            call PBM(Phi)

            er = sum(E**2)
            er0 = sum(Q**2) + 1.0d-16  ! Qが全て0のとき計算できない
            if (sqrt(er / er0) < eps) exit
        enddo
        ! write(*, *) 'Ph:', 'itr =', itr, 'er', er, 'er0', er0 

        ! Phiの平均を0にするため
        B0 = sum(Phi(1:NX, 1:NY))/(NX*NY)
        Phi(:, :) = Phi(:, :) - B0
    end subroutine poisson

    subroutine march(Up, Vp, U, V, Phi, P)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Phi(0:NX+1, 0:NY+1)
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1)
        integer i, j
        do j = 1, NY
            do i = 1, NX
                U(i, j) = Up(i, j) - dt*(-Phi(i, j)+Phi(i+1, j))/dX
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                V(i, j) = Vp(i, j) - dt*(-Phi(i, j)+Phi(i, j+1))/dY
            enddo
        enddo


        do j = 0, NY+1
            do i = 0, NX+1
                P(i, j) = P(i, j) + Phi(i, j)
            enddo
        enddo

        call PBM(U)
        call PBM(V)
        call PBM(P)
    end subroutine march

    subroutine taylor_debug(U, V, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        integer i, j
        real(8) err
        real(8) U0(0:NX+1, 0:NY+1), V0(0:NX+1, 0:NY+1)
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換


        do j = 1, NY
            do i = 1, NX
                U0(i, j) = -sin(i*dX) * cos((j-0.5d0)*dY)
                V0(i, j) = cos((i-0.5d0)*dX) * sin(j*dY)
            enddo
        enddo

        U0(:, :) = U0(:, :) * exp(-2.0d0/Re*step*dt)
        V0(:, :) = V0(:, :) * exp(-2.0d0/Re*step*dt)

        i = NX/4
        j = NY/4
        err = (sum((U0(1:NX, 1:NY)-U(1:NX, 1:NY))**2) &
              +sum((V0(1:NX, 1:NY)-V(1:NX, 1:NY))**2))/(NX*NY)
        open(10, file=trim(dir)//'taylor_debug.d', position='append')
        write(10, *) step, U0(NX/4, NY/4), U(NX/4, NY/4), U0(NX/4, NY/4) - U(NX/4, NY/4), sqrt(err)
        close(10)
    end subroutine taylor_debug


    subroutine get_data(U, V, C, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: C(3, 0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        integer i, j
        real(8) D(2, 2), Omega, S, trC, E(2)
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        open(10, file=trim(dir)//'z_'//str//'.d')

        do j = 1, NY
            do i = 1, NX
                D(1, 1) = (U(i, j) - U(i-1, j))/dX
                D(1, 2) = (U(i, j+1) - U(i, j-1) + U(i-1, j+1) - U(i-1, j-1))/(4*dY)
                D(2, 1) = (V(i+1, j) - V(i-1, j) + V(i+1, j-1) - V(i-1, j-1))/(4*dX)
                D(2, 2) = (V(i, j) - V(i, j-1))/dY
                D(:, :) = D(:, :)*U_C/L_C

                E(1) = (U(i-1, j)+U(i, j))/2*U_C
                E(2) = (V(i, j-1)+V(i, j))/2*U_C
                Omega = (D(2, 1) - D(1, 2))
                S = (D(2, 1) + D(1, 2))
                trC = C(1, i, j)+C(3, i, j)

                write(10, '(100e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, 0.0d0, E(1), E(2), 0.0d0, sum(E**2)/2, &
                                    0.0d0, 0.0d0, Omega, Omega**2/2, &
                                    trC, C(1, i, j), C(2, i, j), 0.0d0, C(3, i, j), 0.0d0, 0.0d0
            enddo
        enddo
        close(10)
    end subroutine get_data


    subroutine logging(U, V, C, step, t_start)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: C(3,0:NX+1, 0:NY+1 )
        integer, intent(in) :: step
        real(8), intent(in) :: t_start
        real(8) t_temp
        real(8) K_energy
        integer i, j, count
        real(8) trC(1:NX, 1:NY), re1, re2
        integer(8) t_pre, hour, min, sec
        real(8) mean, std

        if (mod(step, Dstep) == 0) then
            write(*, '(a, I6)', advance='no') 'step:', step

            ! 運動エネルギーとエネルギーカスケードの総和
            K_energy = sum(U(1:NX, 1:NY)**2) + sum(V(1:NX, 1:NY)**2)
            K_energy = K_energy*U_C**2/2
            K_energy = K_energy/(NX*NY)
            write(*, '(a, e12.4)', advance='no') '  | K_energy:', K_energy

            ! CFL条件
            write(*, '(a, F7.3)', advance='no') &
            '  | CFL:', max(maxval(U(:, :))*dt/dX, maxval(V(:, :))*dt/dY)*4.0d0

            ! write(*, '(a, F6.3)', advance='no') '  | index 0:', counter(0)*1.0d0/sum(counter)
            ! write(*, '(a, F6.3)', advance='no') '  1:', counter(1)*1.0d0/sum(counter)
            ! write(*, '(a, F6.3)', advance='no') '  2:', counter(2)*1.0d0/sum(counter)
            ! write(*, '(a, F6.3)', advance='no') '  3:', counter(3)*1.0d0/sum(counter)
            trC(:, :) = C(1, 1:NX, 1:NY) + C(3, 1:NX, 1:NY)
            mean = sum(trC) / (NX*NY)
            std = sqrt(sum((trC-mean)**2) / (NX*NY))
            write(*, '(a, F8.3, a, F7.3)', advance='no') '  | trC:', mean, ' +-', std
            ! write(*, '(a, F7.3, a, F8.3)', advance='no') '  |', minval(trC(:, :)), '< trC <', maxval(trC(:, :))

            ! count = 0
            ! do j = 1, NY
            !     do i = 1, NX
            !         call Eigenvalues(C(:, i, j), re1, re2)
            !         if (re1 > 0.0d0) count = count + 1
            !     enddo
            ! enddo
            ! write(*, '(a, F7.3)', advance='no') '  | SPD:', count*1.0d0/(NX*NY)

            call cpu_time(t_temp)
            t_pre = int((t_temp-t_start)*(Nstep-step)/step)
            hour = t_pre / 3600
            min = mod(t_pre, 3600) / 60
            sec = mod(t_pre, 60)
            write(*, '(a, I3, a, I2, a, I2)', advance='no') '  | time_left:', hour, ':', min, ':', sec

            write(*, *) ''
        endif
    end subroutine logging


    subroutine input(U, V, P)  ! initを実行した後に書く
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(out) :: P(0:NX+1, 0:NY+1)
        integer i, j
        character(8) str
        real(8) K_energy
        write(str, '(I8.8)') input_step  ! 数値を文字列に変換
        open(10, file=trim(dir)//str//'.d')

        do j = 1, NY
            do i = 1, NX
                read(10, '(4e12.4)') U(i, j), V(i, j), P(i, j)
            enddo
        enddo
        close(10)
        call PBM(U)
        call PBM(V)
        call PBM(P)

        K_energy = sum(U(1:NX, 1:NY)**2) + sum(V(1:NX, 1:NY)**2)
        K_energy = K_energy*U_C**2/2
        K_energy = K_energy/(NX*NY)
        write(*, '(a, a, e12.4)') 'input_file='//trim(dir)//str//'.d', '  | K_energy:', K_energy
    end subroutine input

    subroutine output(U, V, P, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        integer i, j
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換
        open(10, file=trim(dir)//str//'.d')

        do j = 1, NY
            do i = 1, NX
                write(10, '(4e12.4)') U(i, j), V(i, j), P(i, j)
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
    integer(8) plan1, plan2, plan3, plan4
    real(8) Re1(1:NX, 1:NY), Re2(1:NX)
    complex(8) Im1(1:NX/2+1, 1:NY), Im2(1:NX/2+1)
contains
    subroutine fft_init
        call dfftw_plan_dft_r2c_2d(plan1, NX, NY, Re1, Im1, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(plan2, NX, NY, Im1, Re1, FFTW_ESTIMATE)
        ! call dfftw_plan_dft_1d(plan3, NX, Re2, Im2, FFTW_FORWARD, FFTW_ESTIMATE)  ! なぜかエラー
        ! call dfftw_plan_dft_1d(plan4, NX, Im2, Re2, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_plan_dft_r2c_1d(plan3, NX, Re2, Im2, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_1d(plan4, NX, Im2, Re2, FFTW_ESTIMATE)
    end subroutine fft_init

    subroutine fftr2c_2d(input, output)
        real(8), intent(in) :: input(1:NX, 1:NY)
        complex(8), intent(out) :: output(1:NX/2+1, 1:NY)
        Re1(:, :) = input(:, :)
        call dfftw_execute(plan1, Re1, Im1)
        output(:, :) = Im1(:, :)
    end subroutine fftr2c_2d
    
    subroutine fftc2r_2d(input, output)
        complex(8), intent(in) :: input(1:NX/2+1, 1:NY)
        real(8), intent(out) :: output(1:NX, 1:NY)
        Im1(:, :) = input(:, :)
        call dfftw_execute(plan2, Im1, Re1)
        output(:, :) = Re1(:, :)
    end subroutine fftc2r_2d

    subroutine fftr2c_1d(input, output)
        real(8), intent(in) :: input(1:NX)
        complex(8), intent(out) :: output(1:NX/2+1)
        Re2(:) = input(:)
        call dfftw_execute(plan3, Re2, Im2)
        output(:) = Im2(:)
    end subroutine fftr2c_1d

    subroutine fftc2r_1d(input, output)
        complex(8), intent(in) :: input(1:NX/2+1)
        real(8), intent(out) :: output(1:NX)
        Im2(:) = input(:)
        call dfftw_execute(plan4, Im2, Re2)
        output(:) = Re2(:)
    end subroutine fftc2r_1d

    subroutine fft_finalize
        call dfftw_destroy_plan(plan1)
        call dfftw_destroy_plan(plan2)
        call dfftw_destroy_plan(plan3)
        call dfftw_destroy_plan(plan4)
    end subroutine fft_finalize

    
    subroutine fft_solve(Q, LHS, Phi)
        real(8), intent(in) :: Q(1:NX, 1:NY), LHS(1:NX, 1:NY)
        real(8), intent(out) :: Phi(1:NX, 1:NY)
        complex(8) Q_hat(1:NX/2+1, 1:NY), Phi_hat(1:NX/2+1, 1:NY)
        integer i, j
        integer(8) plan

        ! 右辺をフーリエ変換
        call fftr2c_2d(Q, Q_hat)
        
        Q_hat(:, :) = Q_hat(:, :)/(NX*NY)

        ! 定数で割ることで左辺を求める
        do j = 1, NY
            do i = 1, NX/2+1
                if (abs(LHS(i, j)) < 1.0d-16) then  ! マシン零以下
                    Phi_hat(i, j) = (0.0d0, 0.0d0)
                else
                    Phi_hat(i, j) = Q_hat(i, j) / LHS(i, j)
                endif
            enddo
        enddo

        ! 左辺を逆フーリエ変換
        call fftc2r_2d(Phi_hat, Phi)

    end subroutine fft_solve

    subroutine fft_navier(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, &
                          Bx, By, Tx, Ty, Tx0, Ty0, Fx, Fy, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
        real(8), intent(in) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        real(8), intent(in) :: Tx(1:NX, 1:NY), Ty(1:NX, 1:NY)
        real(8), intent(inout) :: Tx0(1:NX, 1:NY), Ty0(1:NX, 1:NY)
        real(8), intent(in) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        integer, intent(in) :: step
        integer i, j
        real(8) Q(1:NX, 1:NY), LHS(1:NX, 1:NY)

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :) = Ax(:, :)
            Ay0(:, :) = Ay(:, :)
            Tx0(:, :) = Tx(:, :)
            Ty0(:, :) = Ty(:, :)
        endif

        ! 左辺の定数部分の計算
        do j = 1, NY
            do i = 1, NX
                LHS(i, j) = 1.0d0 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                                     +2*(1-cos((j-1)*dY))/dY**2 )
            enddo
        enddo

        ! 速度場Upを予測
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = U(i, j) -dt*(-P(i, j)+P(i+1, j))/dX &
                                  +dt*(3.0d0*Ax(i, j) - Ax0(i, j))/2.0d0 &
                                  +dt*Bx(i, j)/2.0d0 &
                                  +dt*(Tx(i, j) + Tx0(i, j))/2.0d0 &
                                  +dt*Fx(i, j)
            enddo
        enddo

        call fft_solve(Q, LHS, Up(1:NX, 1:NY))
        
        ! 同様にVpについても求める
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = V(i, j) -dt*(-P(i, j)+P(i, j+1))/dY &
                                  +dt*(3.0d0*Ay(i, j) - Ay0(i, j))/2.0d0 &
                                  +dt*By(i, j)/2.0d0 &
                                  +dt*(Ty(i, j) + Ty0(i, j))/2.0d0 &
                                  +dt*Fy(i, j)
            enddo
        enddo

        call fft_solve(Q, LHS, Vp(1:NX, 1:NY))
        
        ! Phiを求めるときに0番目の値も必要
        call PBM(Up)
        call PBM(Vp)

        ! n-1ステップ目の保存
        Ax0(:, :) = Ax(:, :)
        Ay0(:, :) = Ay(:, :)
        Tx0(:, :) = Tx(:, :)
        Ty0(:, :) = Ty(:, :)
    end subroutine fft_navier

    subroutine fft_poisson(Up, Vp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: Phi(0:NX+1, 0:NY+1)
        integer i, j
        real(8) Q(1:NX, 1:NY), LHS(1:NX, 1:NY)

        do j = 1, NY
            do i = 1, NX
                LHS(i, j) = -(2*(1-cos((i-1)*dX))/dX**2 &
                             +2*(1-cos((j-1)*dY))/dY**2)
            enddo
        enddo
        
        ! 右辺の計算
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = ((-Up(i-1, j) + Up(i, j)) / dX &
                          +(-Vp(i, j-1) + Vp(i, j)) / dY) / dt
            enddo
        enddo

        call fft_solve(Q, LHS, Phi(1:NX, 1:NY))
        call PBM(Phi)

    end subroutine fft_poisson

    subroutine fft_march(Up, Vp, U, V, Phi, P)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Phi(0:NX+1, 0:NY+1)
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1)
        integer i, j

        do j = 1, NY
            do i = 1, NX
                U(i, j) = Up(i, j) - dt*(-Phi(i, j)+Phi(i+1, j))/dX
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                V(i, j) = Vp(i, j) - dt*(-Phi(i, j)+Phi(i, j+1))/dY
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                P(i, j) = P(i, j) + Phi(i, j) &
                        -dt*beta/Re/2.0d0*((Phi(i-1, j)-2*Phi(i, j)+Phi(i+1, j))/dX**2 &
                                          +(Phi(i, j-1)-2*Phi(i, j)+Phi(i, j+1))/dY**2)
            enddo
        enddo

        call PBM(U)
        call PBM(V)
        call PBM(P)
    end subroutine fft_march

    subroutine energy_single(U, V, step)  ! エネルギースペクトル
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        real(8) U_tmp(0:NX+1, 0:NY+1), V_tmp(0:NX+1, 0:NY+1)
        complex(8) U_hat(1:NX/2+1, 1:NY), V_hat(1:NX/2+1, 1:NY)
        real(8) E_tmp(1:NX/2+1, 1:NY)
        real(8) K_abs(1:NX/2+1, 1:NY)
        real(8) Energy(0:NX)  ! 格子数程度あれば十分
        real(8) K_energy
        integer i, j, index
        integer(8) plan
        character(8) str
        real(8) kx, ky
        real(8) Re(1:NX, 1:NY)
        complex(8) Im(1:NX/2+1, 1:NY)
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        ! そのまま用いると，Energy(0)=0ではないが，そのまま求めたK_energyと一致する．
        U_tmp(:, :) = U(:, :) * U_C
        V_tmp(:, :) = V(:, :) * U_C
        ! 変動成分にすると，Energy(0)=0になり，変動成分を使ったK_energyと一致する．
        ! U_tmp(:, :) = U(:, :) - sum(U(1:NX, 1:NY))/(NX*NY)
        ! V_tmp(:, :) = V(:, :) - sum(V(1:NX, 1:NY))/(NX*NY)
        ! U_tmp(:, :) = U_tmp(:, :) * U_C
        ! V_tmp(:, :) = V_tmp(:, :) * U_C

        ! 速度場をフーリエ変換
        call fftr2c_2d(U_tmp(1:NX, 1:NY), U_hat)
        call fftr2c_2d(V_tmp(1:NX, 1:NY), V_hat)

        ! 格子点数で割る
        U_hat(:, :) = U_hat(:, :)/dble(NX*NY)
        V_hat(:, :) = V_hat(:, :)/dble(NX*NY)

        ! 配列要素に対するエネルギー
        do j = 1, NY
            do i = 1, NX/2+1
                E_tmp(i, j) = abs(U_hat(i, j))**2.0d0 + abs(V_hat(i, j))**2.0d0
            enddo
        enddo

        ! 配列要素に対する波数
        do j = 1, NY
            do i = 1, NX/2+1
                ! 0, 1, 2, ..., 63, 64, -63, -62, ..., -3, -2, -1
                kx = sign(1.0d0, dble(NX + 3)/2.0d0 - dble(i)) * (- dble(abs(NX/2 + 1 - i)) + dble(NX)/2.0d0)  ! 増田さん参考
                ky = sign(1.0d0, dble(NY + 3)/2.0d0 - dble(j)) * (- dble(abs(NY/2 + 1 - j)) + dble(NY)/2.0d0)
                K_abs(i, j) = sqrt(kx**2.0d0 + ky**2.0d0)
            enddo
        enddo

        ! 波数を四捨五入し、対応する整数の波数にエネルギーを足し合わせる
        Energy(:) = 0.0d0
        ! kx = 2 ~ NX/2 は，絶対値が同じ波数が2つあるので2倍する．
        do j = 1, NY
            do i = 1, NX/2+1
                index = nint(K_abs(i, j))
                if (i==1 .or. i==NX/2+1) then
                    Energy(index) = Energy(index) + E_tmp(i, j)/2.d0
                else
                    Energy(index) = Energy(index) + E_tmp(i, j)/2.d0*2.0d0
                endif
            enddo
        enddo

        ! 出力
        ! open(10, file=trim(dir)//'fft_'//str//'.d')
        ! do i = 0, NX
        !     write(10, '(I4, e12.4)') i, Energy(i)
        ! enddo
        ! close(10)


        ! 運動エネルギーとエネルギーカスケードの総和
        K_energy = sum(U_tmp(1:NX, 1:NY)**2.0d0) + sum(V_tmp(1:NX, 1:NY)**2.0d0) 
        K_energy = K_energy/2.d0
        K_energy = K_energy/(NX*NY)
        ! write(*, '(I6, 5e12.4)') step, K_energy, sum(Energy(:)), K_energy-sum(Energy(:)), Energy(0), K_energy/sum(Energy(:))   ! 一致するはず
        ! open(30, file = 'fft/debug_energy.d', position='append')
        ! write(30, '(I6, 2e12.4)') step, K_energy, sum(Energy(:))
        ! close(30)

    end subroutine energy_single

end module fft

module kabe
    use smac
    use fft
    implicit none
    real(8) Uwall  ! 壁面速度
    integer, parameter :: NA = NY/4 ! 下壁
    integer, parameter :: NB = 3*NY/4 ! 上壁
contains
    subroutine kabe_init(U, V, P, C)
        real(8), intent(inout) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1), P(0:NX+1, 0:NY+1), C(3, 0:NX+1, 0:NY+1)
        if (flow_type == 1) Uwall = 1.0d0
        if (flow_type == 2) Uwall = 0.0d0
        call U_kabe(U)
        call V_kabe(V)
        call P_kabe(P)
        call C_kabe(C)
    end subroutine kabe_init

    subroutine U_kabe(U)
        real(8), intent(inout) :: U(0:NX+1, 0:NY+1)
        U(0, :) = U(NX, :)
        U(NX+1, :) = U(1, :)
        U(:, NA) = -U(:, NA+1)
        U(:, NB+1) = 2*Uwall - U(:, NB)
        ! その他0
        U(:, 0:NA-1) = 0.0d0
        U(:, NB+2:NY+1) = 0.0d0
    end subroutine U_kabe

    subroutine V_kabe(V)
        real(8), intent(inout) :: V(0:NX+1, 0:NY+1)
        V(0, :) = V(NX, :)
        V(NX+1, :) = V(1, :)
        V(:, NA) = 0.0d0
        V(:, NB) = 0.0d0
        ! その他0
        V(:, 0:NA-1) = 0.0d0
        V(:, NB+1:NY+1) = 0.0d0
    end subroutine V_kabe

    subroutine P_kabe(P)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1)
        P(0, :) = P(NX, :)
        P(NX+1, :) = P(1, :)
        P(:, NA) = P(:, NA+1)
        P(:, NB+1) = P(:, NB)
        ! その他
        P(:, 0:NA-1) = 0.0d0
        P(:, NB+2:NY+1) = 0.0d0
    end subroutine P_kabe

    subroutine C_kabe(C)
        real(8), intent(inout) :: C(3, 0:NX+1, 0:NY+1)
        C(:, 0, :) = C(:, NX, :)
        C(:, NX+1, :) = C(:, 1, :)
        C(:, :, NA) = 2*C(:, :, NA+1) - C(:, :, NA+2)
        C(:, :, NB+1) = 2*C(:, :, NB) - C(:, :, NB-1)
        ! その他0
        C(:, :, 0:NA-1) = 0.0d0
        C(:, :, NB+2:NY+1) = 0.0d0
        ! 下のつけると補間した値ばかり使われるようになるかもしれない。あえて0にしてSPDを満たさないようにしてる
        ! C(:, :, NA-1) = 2*C(:, :, NA) - C(:, :, NA+1)
        ! C(:, :, NB+2) = 2*C(:, :, NB+1) - C(:, :, NB)
    end subroutine C_kabe

    subroutine tdma(a, b, c, d, x, size)  ! 三重対角行列
        implicit none
        integer, intent(in) :: size
        complex(8), intent(in) :: a(size), b(size), c(size), d(size)
        complex(8), intent(out) :: x(size)
        complex(8) :: P(size), Q(size)
        integer :: i
      
        P(1) = -b(1)/a(1)
        Q(1) = d(1)/a(1)
      
        do i = 2, size
            if (abs(a(i)+c(i)*P(i-1))<1.0d-16) then  ! マシン零以下
                P(i) = 0.0d0
                Q(i) = 0.0d0
            else
                P(i) = -b(i)/(a(i)+c(i)*P(i-1))
                Q(i) = (d(i)-c(i)*Q(i-1))/(a(i)+c(i)*P(i-1))
            endif
            ! write(*, *) 'i', i, P(i), Q(i)
        enddo
        ! write(*, *) ''

        x(size) = Q(size)
      
        do i = size-1, 1, -1
            x(i) = P(i)*x(i+1) + Q(i)
            ! write(*, *) 'i', i, x(i)
        enddo
        ! write(*, *) x(:)
        ! write(*, *) ''
      
      end subroutine tdma

    subroutine fft_navier_kabe(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, &
                               Bx, By, Tx, Ty, Tx0, Ty0, Fx, Fy, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
        real(8), intent(in) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        real(8), intent(in) :: Tx(1:NX, 1:NY), Ty(1:NX, 1:NY)
        real(8), intent(inout) :: Tx0(1:NX, 1:NY), Ty0(1:NX, 1:NY)
        real(8), intent(in) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        integer, intent(in) :: step
        integer i, j
        real(8) Q(1:NX, NA+1:NB)
        complex(8) Q_hat(1:NX/2+1, NA+1:NB)
        complex(8) Up_hat(1:NX/2+1, NA:NB+1), Vp_hat(1:NX/2+1, NA:NB+1)
        integer, parameter :: size = NB-NA+2
        complex(8) a(size), b(size), c(size), d(size), x(size)

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :) = Ax(:, :)
            Ay0(:, :) = Ay(:, :)
            Tx0(:, :) = Tx(:, :)
            Ty0(:, :) = Ty(:, :)
        endif

        ! 速度場Upを予測
        do j = NA+1, NB
            do i = 1, NX
                Q(i, j) = U(i, j) -dt*(-P(i, j)+P(i+1, j))/dX &
                                  +dt*(3.0d0*Ax(i, j) - Ax0(i, j))/2.0d0 &
                                  +dt*Bx(i, j)/2.0d0 &
                                  +dt*(Tx(i, j) + Tx0(i, j))/2.0d0 &
                                  +dt*Fx(i, j)
            enddo
        enddo

        ! 右辺をフーリエ変換
        do j = NA+1, NB
            call fftr2c_1d(Q(:, j), Q_hat(:, j))
        enddo
        Q_hat(:, :) = Q_hat(:, :)/NX

        ! 三重対角行列を解く
        do i = 1, NX/2+1
            a(2:size-1) = 1.0d0 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 + 2.0d0/dY**2)
            b(2:size-1) = - dt*beta/Re/2.0d0 / dY**2
            c(2:size-1) = - dt*beta/Re/2.0d0 / dY**2
            d(2:size-1) = Q_hat(i, :)
            a(1) = 1.0d0
            b(1) = 1.0d0
            d(1) = 0.0d0
            a(size) = 1.0d0
            c(size) = 1.0d0
            ! if (i==1) d(size) = 2*Uwall*NX  ! δ(0)は格子点数になるはず。でも格子点を変えたらその分大きくなる。
            if (i==1) d(size) = 2*Uwall
            if (i/=1) d(size) = 0.0d0
            call tdma(a, b, c, d, x, size)
            Up_hat(i, :) = x(:)
        enddo
        ! 左辺を逆フーリエ変換
        do j = NA, NB+1
            call fftc2r_1d(Up_hat(:, j), Up(1:NX, j))
        enddo

        
        ! 同様にVpについても求める
        do j = NA+1, NB  ! NBは壁面
            do i = 1, NX
                Q(i, j) = V(i, j) -dt*(-P(i, j)+P(i, j+1))/dY &
                                  +dt*(3.0d0*Ay(i, j) - Ay0(i, j))/2.0d0 &
                                  +dt*By(i, j)/2.0d0 &
                                  +dt*(Ty(i, j) + Ty0(i, j))/2.0d0 &
                                  +dt*Fy(i, j)
            enddo
        enddo

        ! 右辺をフーリエ変換
        do j = NA+1, NB
            call fftr2c_1d(Q(:, j), Q_hat(:, j))
        enddo
        Q_hat(:, :) = Q_hat(:, :)/NX

        ! 三重対角行列を解く
        do i = 1, NX/2+1
            a(2:size-2) = 1.0d0 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 + 2.0d0/dY**2)
            b(2:size-2) = - dt*beta/Re/2.0d0 / dY**2
            c(2:size-2) = - dt*beta/Re/2.0d0 / dY**2
            d(2:size-1) = Q_hat(i, :)
            a(1) = 1.0d0
            b(1) = 0.0d0
            d(1) = 0.0d0
            a(size-1:size) = 1.0d0
            b(size-1:size) = 0.0d0
            c(size-1:size) = 0.0d0
            d(size-1:size) = 0.0d0
            call tdma(a, b, c, d, x, size)
            Vp_hat(i, :) = x(:)
        enddo
        ! 左辺を逆フーリエ変換
        do j = NA, NB+1
            call fftc2r_1d(Vp_hat(:, j), Vp(1:NX, j))
        enddo

        ! 予測速度に境界条件使っちゃう
        call U_kabe(Up)
        call V_kabe(Vp)

        ! n-1ステップ目の保存
        Ax0(:, :) = Ax(:, :)
        Ay0(:, :) = Ay(:, :)
        Tx0(:, :) = Tx(:, :)
        Ty0(:, :) = Ty(:, :)
    end subroutine fft_navier_kabe

    subroutine fft_poisson_kabe(Up, Vp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: Phi(0:NX+1, 0:NY+1)
        integer i, j
        real(8) Q(1:NX, NA+1:NB)
        complex(8) Q_hat(1:NX/2+1, NA+1:NB)
        complex(8) Phi_hat(1:NX/2+1, NA:NB+1)
        integer, parameter :: size = NB-NA+2
        complex(8) a(size), b(size), c(size), d(size), x(size)

        ! 右辺の計算
        do j = NA+1, NB
            do i = 1, NX
                Q(i, j) = ((-Up(i-1, j) + Up(i, j)) / dX &
                          +(-Vp(i, j-1) + Vp(i, j)) / dY) / dt
            enddo
        enddo

        ! 右辺をフーリエ変換
        do j = NA+1, NB
            call fftr2c_1d(Q(:, j), Q_hat(:, j))
        enddo
        Q_hat(:, :) = Q_hat(:, :)/NX

        ! 三重対角行列を解く
        do i = 1, NX/2+1
            a(2:size-1) = -(2.0d0*(1.0d0-cos((i-1)*dX))/dX**2 + 2.0d0/dY**2)
            b(2:size-1) = 1.0d0 / dY**2
            c(2:size-1) = 1.0d0 / dY**2
            d(2:size-1) = Q_hat(i, :)
            a(1) = 1.0d0
            b(1) = -1.0d0
            d(1) = 0.0d0
            a(size) = 1.0d0
            c(size) = -1.0d0
            d(size) = 0.0d0
            call tdma(a, b, c, d, x, size)
            Phi_hat(i, :) = x(:)
        enddo

        ! 左辺を逆フーリエ変換
        do j = NA, NB+1
            call fftc2r_1d(Phi_hat(:, j), Phi(1:NX, j))
        enddo

        call P_kabe(Phi)
        ! Phi(:, :) = Phi(:, :) - sum(Phi(1:NX, NA+1:NB))/(NX*(NB-NA))

    end subroutine fft_poisson_kabe

    subroutine fft_march_kabe(Up, Vp, U, V, Phi, P, C)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Phi(0:NX+1, 0:NY+1)
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: C(3, 0:NX+1, 0:NY+1)
        integer i, j

        do j = NA+1, NB
            do i = 1, NX
                U(i, j) = Up(i, j) - dt*(-Phi(i, j)+Phi(i+1, j))/dX
            enddo
        enddo

        do j = NA+1, NB
            do i = 1, NX
                V(i, j) = Vp(i, j) - dt*(-Phi(i, j)+Phi(i, j+1))/dY
            enddo
        enddo

        do j = NA+1, NB
            do i = 1, NX
                P(i, j) = P(i, j) + Phi(i, j) &
                        -dt*beta/Re/2.0d0*((Phi(i-1, j)-2*Phi(i, j)+Phi(i+1, j))/dX**2 &
                                          +(Phi(i, j-1)-2*Phi(i, j)+Phi(i, j+1))/dY**2)
            enddo
        enddo
        call U_kabe(U)
        call V_kabe(V)
        call P_kabe(P)
        call C_kabe(C)

        ! do j = NA, NB+1
        !     write(*, *) sum(U(:, j))
        ! enddo
        ! write(*, *) ''
        ! do j = NA, NB
        !     write(*, *) sum(V(:, j))
        ! enddo
    end subroutine fft_march_kabe

    subroutine koide_theory(U, V, C)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: C(3, 0:NX+1, 0:NY+1)
        real(8) Fp, Fn, Gp, Gn, B, A, tmp, Re_tmp, fx, y
        real(8) Fp_tmp, Fn_tmp, Gp_tmp, Gn_tmp
        real(8) u_theory(NA+1:NB)
        integer i, j

        ! はじめにy=-1での値が必要
        fx = 1.0d0
        y = -1.0d0
        A = Lp**2/(6.0d0*beta*Wi**2)
        Re_tmp = Re *(PI/2.0d0)
        B = -Lp**2*fx*Re_tmp/(4.0d0*beta*Wi**2)

        Fp_tmp = (B*y + sqrt(B**2*y**2+A**3))**(1.0d0/3.0d0)
        tmp = B*y - sqrt(B**2*y**2+A**3)
        Fn_tmp = sign(1.0d0, tmp)*abs(tmp)**(1.0d0/3.0d0)
        Gp_tmp = 3*B*y + sqrt(B**2*y**2+A**3)
        Gn_tmp = 3*B*y - sqrt(B**2*y**2+A**3)

        open(12, file=trim(dir)//'koide_theory.d')
        do j = NA+1, NB
            y = ((j-0.5d0)*dY - PI)/(PI/2.0d0)  ! (PI/2.0d0)はスケールを合わせるため
            A = Lp**2/(6.0d0*beta*Wi**2)
            Re_tmp = Re *(PI/2.0d0)
            B = -Lp**2*fx*Re_tmp/(4.0d0*beta*Wi**2)

            Fp = (B*y + sqrt(B**2*y**2+A**3))**(1.0d0/3.0d0)
            tmp = B*y - sqrt(B**2*y**2+A**3)
            Fn = sign(1.0d0, tmp)*abs(tmp)**(1.0d0/3.0d0)  ! 負の数の三乗根が計算できない
            Gp = 3*B*y + sqrt(B**2*y**2+A**3)
            Gn = 3*B*y - sqrt(B**2*y**2+A**3)

            u_theory(j) = Re_tmp/(2.0d0*beta)*fx*(1.0d0-y**2) - &
                          (1.0d0-beta)/beta*3.0d0/(8.0d0*B)*(Fp*Gn-Fp_tmp*Gn_tmp + Fn*Gp-Fn_tmp*Gp_tmp)
            u_theory(j) = u_theory(j) * (PI/2.0d0)

            write(12, '(12e12.4)') (j-0.5d0)*dY_C, u_theory(j)*U_C, sum(U(1:NX, j) - u_theory(j))/NX*U_C
        enddo
        close(12)

        tmp = 0.0d0
        do j = NA+1, NB
            do i = 1, NX
                tmp = tmp + (U(i, j) - u_theory(j))**2
            enddo
        enddo
        tmp = sqrt(tmp/(NX*(NB-NA))) * U_C
        write(*, '(a, E12.4, a, E20.12)') 'dX:', dX, '  err:', tmp


    end subroutine koide_theory
    
end module kabe

module ibm
    use smac
    use fft
    implicit none
    ! IBM用パラメータ
    integer NC, NL
    real(8) dV
    character(64) :: dir_ibm = './ibm/'
contains
    subroutine ibm_init(X, Y, Xc, Yc, Uc, Vc, Fxc, Fyc)
        real(8), intent(out) :: X(0:NX+1, 0:NY+1), Y(0:NX+1, 0:NY+1)
        real(8), allocatable :: Xc(:, :), Yc(:, :)
        real(8), allocatable :: Uc(:, :), Vc(:, :)
        real(8), allocatable :: Fxc(:, :), Fyc(:, :)
        integer i, j, index
        integer NL_lap  ! 円柱内部の殻数
        integer, allocatable :: NL_shell(:)  ! 殻に対する外力点の数
        character(2) str
        write(str, '(I2.2)') ibm_type
        dir_ibm = trim(dir_ibm)//'ibm'//str//'/'
        write(*, '(a, a)') 'dir_ibm:   ', dir_ibm
        call mk_dir(trim(dir_ibm))

        if (ibm_type == 1) NC = 2  ! 壁表面
        if (ibm_type == 2) NC = NY/2+1  ! 壁の中全て
        if (ibm_type == 3) NC = NY/4+1  ! 壁の中半分
        if (ibm_type < 10) then
            NL = NX
            dV = dX*dX
        endif

        if (ibm_type == 11 .or. ibm_type == 21) NC = 1  ! 円柱の数
        if (ibm_type == 12 .or. ibm_type == 22) NC = 4

        if (ibm_type == 11 .or. ibm_type == 12) then  ! 外力点が円柱表面のみの場合
            NL = nint(PI*DC/dX)
            dV = PI*DC/NL * dX
        endif
        if (ibm_type == 21 .or. ibm_type == 22) then  ! 外力点が円柱内部にもある場合
            NL_lap = int(DC/(2*dX))  ! 円柱内部の殻数
            allocate(NL_shell(0:NL_lap))
            do i = 0, NL_lap
                NL_shell(i) = nint(PI*(DC-2*i*dX)/dX)
            enddo
            if (NL_shell(NL_lap) == 0) NL_shell(NL_lap) = 1  ! 中心に外力点を1つ置く
            NL = sum(NL_shell)
            dV = PI*(DC+dX)**2/4 / NL
            write(*, '(a, 100I5)') 'NL_shell:  ', NL_shell(:)
        endif

        ! NLを定義したあとallocate
        allocate(Xc(1:NC, 1:NL), Yc(1:NC, 1:NL))
        allocate(Uc(1:NC, 1:NL), Vc(1:NC, 1:NL))
        allocate(Fxc(1:NC, 1:NL), Fyc(1:NC, 1:NL))

        ! コメントとファイル
        if (dX /= dY) then
            stop 'dX and dY are not equal'
        endif

        write(*, '(a, I5)') 'NL:  ', NL
        write(*, '(a, E12.4)') 'dX**2:', dX**2
        write(*, '(a, E12.4)') 'dV   :', dV


        ! 格子点中心座標
        do i = 0, NX+1
            X(i, :) = (i-0.5d0)*dX
        enddo
        do j = 0, NY+1
            Y(:, j) = (j-0.5d0)*dY
        enddo

        ! 外力点の定義
        Uc(:, :) = 0.0d0
        Vc(:, :) = 0.0d0
        if (ibm_type == 1) then
            do j = 1, NL
                Xc(:, j) = (j-0.25d0)*dX
            enddo
            Yc(1, :) = 2*PI/4
            Yc(2, :) = 2*PI*3/4
            if (flow_type == 1) then
                Uc(1, :) = 0.0d0
                Uc(2, :) = 1.0d0
            endif
        endif
        if (ibm_type == 2) then
            do j = 1, NL
                Xc(:, j) = (j-0.25d0)*dX
            enddo
            do i = 1, NY/4+1
                Yc(i, :) = (i-1)*dY
            enddo
            do i = NY/4+2, NY/2+1
                Yc(i, :) = (i-2)*dY + PI
            enddo
            if (flow_type == 1) then
                do i = 1, NY/4+1
                    Uc(i, :) = 0.5d0 - (i-0.5d0)*0.5d0/(NY/4)
                enddo
                do i = NY/4+2, NY/2+1
                    Uc(i, :) = 1.5d0 - (i-1.5d0)*0.5d0/(NY/4)
                enddo
            endif
        endif
        if (ibm_type == 3) then
            do j = 1, NL
                Xc(:, j) = (j-0.25d0)*dX
            enddo
            do i = 1, NY/8+1
                Yc(i, :) = (i-1)*dY + 2*PI/8
            enddo
            do i = NY/8+2, NY/4+1
                Yc(i, :) = (i-2)*dY + 2*PI*5/8
            enddo
            if (flow_type == 1) then
                do i = 1, NY/8+1
                    Uc(i, :) = 0.0d0
                enddo
                do i = NY/8+2, NY/4+1
                    Uc(i, :) = 1.0d0
                enddo
            endif
        endif
        if (ibm_type == 11) then
            do j = 1, NL
                Xc(1, j) = 2*PI/2 + DC/2 * cos(2*PI/NL*j)  ! 中心に配置
                Yc(1, j) = 2*PI/2 + DC/2 * sin(2*PI/NL*j)
            enddo
        endif
        if (ibm_type == 12) then
            do j = 1, NL
                Xc(1, j) = 2*PI/4 + DC/2 * cos(2*PI/NL*j)
                Yc(1, j) = 2*PI/4 + DC/2 * sin(2*PI/NL*j)
                Xc(2, j) = 2*PI*3/4 + DC/2 * cos(2*PI/NL*j)
                Yc(2, j) = 2*PI/4 + DC/2 * sin(2*PI/NL*j)
                Xc(3, j) = 2*PI/4 + DC/2 * cos(2*PI/NL*j)
                Yc(3, j) = 2*PI*3/4 + DC/2 * sin(2*PI/NL*j)
                Xc(4, j) = 2*PI*3/4 + DC/2 * cos(2*PI/NL*j)
                Yc(4, j) = 2*PI*3/4 + DC/2 * sin(2*PI/NL*j)
                Uc(1, j) = sin(2*PI/NL*j)
                Vc(1, j) = -cos(2*PI/NL*j)
                Uc(2, j) = -sin(2*PI/NL*j)
                Vc(2, j) = cos(2*PI/NL*j)
                Uc(3, j) = -sin(2*PI/NL*j)
                Vc(3, j) = cos(2*PI/NL*j)
                Uc(4, j) = sin(2*PI/NL*j)
                Vc(4, j) = -cos(2*PI/NL*j)
            enddo
        endif
        if (ibm_type == 21) then
            index = 0
            do i = 0, NL_lap
                do j = 1, NL_shell(i)
                    index = index + 1
                    Xc(1, index) = 2*PI/2   + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc(1, index) = 2*PI/2   + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Uc(1, index) =  (DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc(1, index) = -(DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                enddo
            enddo
        endif
        if (ibm_type == 22) then
            index = 0
            do i = 0, NL_lap
                do j = 1, NL_shell(i)
                    index = index + 1
                    Xc(1, index) = 2*PI/4   + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc(1, index) = 2*PI/4   + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Xc(2, index) = 2*PI*3/4 + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc(2, index) = 2*PI/4   + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Xc(3, index) = 2*PI/4   + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc(3, index) = 2*PI*3/4 + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Xc(4, index) = 2*PI*3/4 + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc(4, index) = 2*PI*3/4 + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Uc(1, index) =  (DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc(1, index) = -(DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                    Uc(2, index) = -(DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc(2, index) =  (DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                    Uc(3, index) = -(DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc(3, index) =  (DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                    Uc(4, index) =  (DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc(4, index) = -(DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                enddo
            enddo
        endif
        Xc(:, :) = Xc(:, :) + 10d-10*dX  ! 速度定義点から少しずらして外力点を定義
        Yc(:, :) = Yc(:, :) + 10d-10*dX

    end subroutine ibm_init

    subroutine ibm_get(Xc, Yc)
        real(8), intent(in) :: Xc(1:NC, 1:NL), Yc(1:NC, 1:NL)
        integer i, j
        character(8) str

        do i = 1, NC
            write(str, '(I8.8)') i
            open(10, file=trim(dir_ibm)//str//'.d')
            do j = 1, NL
                write(10, '(2e12.4)') Xc(i, j)*L_C, Yc(i, j)*L_C
            enddo
            write(10, '(2e12.4)') Xc(i, 1)*L_C, Yc(i, 1)*L_C  ! 1周させる
            close(10)
        enddo
    end subroutine ibm_get

    subroutine ibm_preliminary(U, V, P, Ax, Ay, Ax0, Ay0, Bx, By, Tx, Ty, Tx0, Ty0, Fx, Fy, &
                               Ua, Va, Fxc, Fyc, Up, Vp, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
        real(8), intent(in) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        real(8), intent(in) :: Tx(1:NX, 1:NY), Ty(1:NX, 1:NY)
        real(8), intent(in) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        real(8), intent(inout) :: Tx0(1:NX, 1:NY), Ty0(1:NX, 1:NY)
        real(8), intent(out) :: Ua(0:NX+1, 0:NY+1), Va(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: Fxc(1:NC, 1:NL), Fyc(1:NC, 1:NL)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        integer i, j
        
        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :) = Ax(:, :)
            Ay0(:, :) = Ay(:, :)
            Tx0(:, :) = Tx(:, :)
            Ty0(:, :) = Ty(:, :)
        endif

        do j = 1, NY
            do i = 1, NX
                Ua(i, j) = U(i, j) - dt*(-P(i, j)+P(i+1, j))/dX &
                                    +dt*(3.0d0*Ax(i, j) - Ax0(i, j))/2.0d0 &
                                    +dt*Bx(i, j) &
                                    +dt*(Tx(i, j) + Tx0(i, j))/2.0d0 &
                                    +dt*Fx(i, j)

                Va(i, j) = V(i, j) - dt*(-P(i, j)+P(i, j+1))/dY &
                                    +dt*(3.0d0*Ay(i, j) - Ay0(i, j))/2.0d0 &
                                    +dt*By(i, j) &
                                    +dt*(Ty(i, j) + Ty0(i, j))/2.0d0 &
                                    +dt*Fy(i, j)
            enddo
        enddo

        ! 周期境界
        call PBM(Ua)
        call PBM(Va)

        ! n-1ステップ目の保存
        Ax0(:, :) = Ax(:, :)
        Ay0(:, :) = Ay(:, :)
        Tx0(:, :) = Tx(:, :)
        Ty0(:, :) = Ty(:, :)

        ! 円柱上での外力初期化
        Fxc(:, :) = 0.0d0
        Fyc(:, :) = 0.0d0

        ! 反復回数によらずHelmholtzの式を使うため
        Up(:, :) = Ua(:, :)
        Vp(:, :) = Va(:, :)

    end subroutine ibm_preliminary

    function delta(x, y, xc, yc) result(out)
        real(8), intent(in) :: x, y, xc, yc
        real(8) out
        out = weight((x-xc)/dX) * weight((y-yc)/dY) /(dX*dY)
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

    subroutine ibm_Helmholtz(Up, Vp, X, Y, Xc, Yc, Uc, Vc, Fxc, Fyc, fxint, fyint, step)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: X(0:NX+1, 0:NY+1), Y(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Xc(1:NC, 1:NL), Yc(1:NC, 1:NL)
        real(8), intent(in) :: Uc(1:NC, 1:NL), Vc(1:NC, 1:NL)
        real(8), intent(inout) :: Fxc(1:NC, 1:NL), Fyc(1:NC, 1:NL)
        real(8), intent(out) :: fxint(1:NX, 1:NY), fyint(1:NX, 1:NY)
        integer, intent(in)  :: step
        integer i, j, l, m
        real(8) Ub(1:NC, 1:NL), Vb(1:NC, 1:NL)
        real(8) fxtmp(0:NX+1, 0:NY+1), fytmp(0:NX+1, 0:NY+1)
        ! real(8) er
        ! integer count2d(0:NX+1, 0:NY+1), countint(1:NX, 1:NY)
        ! count2d(:, :) = 0
        
        Ub(:, :) = 0.0d0
        Vb(:, :) = 0.0d0
        do m = 1, NL
            do l = 1, NC
                ! 外力点周囲の3*3点の和をとる
                do j = int(Yc(l, m)/dY), int(Yc(l, m)/dY) + 2
                    do i = int(Xc(l, m)/dX - 0.5d0), int(Xc(l, m)/dX - 0.5d0) + 2
                        Ub(l, m) = Ub(l, m) &
                                 + Up(i, j) * delta(X(i, j)+0.5d0*dX, Y(i, j), Xc(l, m), Yc(l, m)) * dX**2
                    enddo
                enddo
                do j = int(Yc(l, m)/dY - 0.5d0), int(Yc(l, m)/dY - 0.5d0) + 2
                    do i = int(Xc(l, m)/dX), int(Xc(l, m)/dX) + 2
                        Vb(l, m) = Vb(l, m) &
                                 + Vp(i, j) * delta(X(i, j), Y(i, j)+0.5d0*dY, Xc(l, m), Yc(l, m)) * dX**2
                    enddo
                enddo
            enddo
        enddo

        do m = 1, NL
            do l = 1, NC
                Fxc(l, m) = Fxc(l, m) + (Uc(l, m) - Ub(l, m))/dt
                Fyc(l, m) = Fyc(l, m) + (Vc(l, m) - Vb(l, m))/dt
            enddo
        enddo

        ! if(mod(step, Dstep) == 0) write(*, '(a, E12.4)') 'sum(Uc-U):', sum((Uc(:, :) - Ub(:, :)))

        fxtmp(:, :) = 0.0d0
        fytmp(:, :) = 0.0d0
        do m = 1, NL
            do l = 1, NC
                ! 格子点周辺の外力点の和をとる、つまり外力点周辺の格子点へ和をのこす
                do j = int(Yc(l, m)/dY), int(Yc(l, m)/dY) + 2
                    do i = int(Xc(l, m)/dX - 0.5d0), int(Xc(l, m)/dX - 0.5d0) + 2
                        fxtmp(i, j) = fxtmp(i, j) &
                                    + Fxc(l, m) * delta(X(i, j)+0.5d0*dX, Y(i, j), Xc(l, m), Yc(l, m)) * dX**2
                        ! count2d(i, j) = count2d(i, j) + 1
                    enddo
                enddo
                do j = int(Yc(l, m)/dY - 0.5d0), int(Yc(l, m)/dY - 0.5d0) + 2
                    do i = int(Xc(l, m)/dX), int(Xc(l, m)/dX) + 2
                        fytmp(i, j) = fytmp(i, j) &
                                    + Fyc(l, m) * delta(X(i, j), Y(i, j)+0.5d0*dY, Xc(l, m), Yc(l, m)) * dX**2
                    enddo
                enddo
            enddo
        enddo

        ! do j = 0, NY+1
        !     write(*, '(100I4)') count2d(:, j)
        ! enddo
        ! write(*, *) ''
        ! count2d(1, :) = count2d(1, :) + count2d(NX+1, :)
        ! count2d(NX, :) = count2d(NX, :) + count2d(0, :)
        ! count2d(:, 1) = count2d(:, 1) + count2d(:, NY+1)
        ! count2d(:, NY) = count2d(:, NY) + count2d(:, 0)
        ! do j = 1, NY
        !     do i = 1, NX
        !         write(*, '(100I4)', advance='no') count2d(i, j)
        !     enddo
        !     write(*, *) ''
        ! enddo

        fxtmp(1, :) = fxtmp(1, :) + fxtmp(NX+1, :)
        fxtmp(NX, :) = fxtmp(NX, :) + fxtmp(0, :)
        fxtmp(:, 1) = fxtmp(:, 1) + fxtmp(:, NY+1)
        fxtmp(:, NY) = fxtmp(:, NY) + fxtmp(:, 0)

        fytmp(1, :) = fytmp(1, :) + fytmp(NX+1, :)
        fytmp(NX, :) = fytmp(NX, :) + fytmp(0, :)
        fytmp(:, 1) = fytmp(:, 1) + fytmp(:, NY+1)
        fytmp(:, NY) = fytmp(:, NY) + fytmp(:, 0)

        fxint(:, :) = fxtmp(1:NX, 1:NY)
        fyint(:, :) = fytmp(1:NX, 1:NY)

    end subroutine ibm_Helmholtz

    subroutine ibm_predict(Ua, Va, fxint, fyint, Bx, By, Up, Vp)
        real(8), intent(in) :: Ua(0:NX+1, 0:NY+1), Va(0:NX+1, 0:NY+1)
        real(8), intent(in) :: fxint(1:NX, 1:NY), fyint(1:NX, 1:NY)
        real(8), intent(in) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        integer i, j
        real(8) Q(1:NX, 1:NY), LHS(1:NX, 1:NY)

        ! 左辺の定数部分の計算
        do j = 1, NY
            do i = 1, NX
                LHS(i, j) = 1.0d0 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                                    + 2*(1-cos((j-1)*dY))/dY**2)
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                Q(i, j) = Ua(i, j) + dt*fxint(i, j) - dt*Bx(i, j)/2.0d0
            enddo
        enddo

        call fft_solve(Q, LHS, Up(1:NX, 1:NY))

        do j = 1, NY
            do i = 1, NX
                Q(i, j) = Va(i, j) + dt*fyint(i, j) - dt*By(i, j)/2.0d0
            enddo
        enddo
        call fft_solve(Q, LHS, Vp(1:NX, 1:NY))

        ! do j = 1, NY
        !     do i = 1, NX
        !         Up(i, j) = Ua(i, j) + dt*fxint(i, j)
        !     enddo
        ! enddo

        ! do j = 1, NY
        !     do i = 1, NX
        !         Vp(i, j) = Va(i, j) + dt*fyint(i, j)
        !     enddo
        ! enddo

        ! Phiを求めるときに0番目の値も必要
        call PBM(Up)
        call PBM(Vp)
    end subroutine ibm_predict
end module ibm

program main
    use smac
    use fft
    use kabe
    use ibm
    implicit none
    real(8) U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
    real(8) C(3, 0:NX+1, 0:NY+1)
    real(8) Cx(3, NX, NY)
    real(8) Cpx(3, 0:NX, NY), Cnx(3, 0:NX, NY)
    real(8) Cpy(3, NX, 0:NY), Cny(3, NX, 0:NY)
    real(8) Tx(1:NX, 1:NY), Ty(1:NX, 1:NY)
    real(8) Tx0(1:NX, 1:NY), Ty0(1:NX, 1:NY)
    real(8) Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)  ! 対流項の計算
    real(8) Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
    real(8) Bx(1:NX, 1:NY), By(1:NX, 1:NY)  ! 粘性項の計算
    real(8) Bx0(1:NX, 1:NY), By0(1:NX, 1:NY)
    real(8) Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
    real(8) P(0:NX+1, 0:NY+1)
    real(8) Phi(0:NX+1, 0:NY+1)
    real(8) Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
    ! ibmのために追加した変数
    real(8) X(0:NX+1, 0:NY+1), Y(0:NX+1, 0:NY+1)
    real(8), allocatable :: Xc(:, :), Yc(:, :)
    real(8), allocatable :: Uc(:, :), Vc(:, :)
    real(8), allocatable :: Fxc(:, :), Fyc(:, :)
    real(8) Ua(0:NX+1, 0:NY+1), Va(0:NX+1, 0:NY+1)
    real(8) fxint(1:NX, 1:NY), fyint(1:NX, 1:NY)
    integer step, i, j
    real(8) t1, t2, t3, t4, t5, t12, t23, t34, t45
    real(8) total_time, t_start, t_end, others_time
    call cpu_time(t_start)
    t12 = 0.0d0
    
    call fft_init
    call init(U, V, P, Phi, C, Fx, Fy)  ! 初期条件
    if (input_step > 0) call input(U, V, P)
    ! call energy_single(U, V, W, 0)
    if (method == 2) call ibm_init(X, Y, Xc, Yc, Uc, Vc, Fxc, Fyc)
    if (method == 2) call ibm_get(Xc, Yc)
    if (method == 3) call kabe_init(U, V, P, C)

    do step = 1, Nstep
        call cpu_time(t1)
        if (beta == 1.0d0) then
            C(:, :, :) = 0.0d0
            Tx(:, :) = 0.0d0
            Ty(:, :) = 0.0d0
            call cpu_time(t2)
            call cpu_time(t3)
        else
            call CpxCnx(C, Cpx, Cnx)
            call CpyCny(C, Cpy, Cny)

            call cpu_time(t2)
            call Cstar(Cpx, Cnx, Cpy, Cny, U, V, Cx)
            call Lyapunov(Cx, U, V, C)
            call cpu_time(t3)

            call polymer_stress(C, Tx, Ty)
        endif
        call convection(U, V, Ax, Ay)
        call viscous(U, V, Bx, By)

        if (method == 0) then
            call navier(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, &
                        Bx, By, Bx0, By0, Tx, Ty, Tx0, Ty0, Fx, Fy, step)
            call cpu_time(t4)
            call poisson(Up, Vp, Phi)
            call march(Up, Vp, U, V, Phi, P)
        endif

        if (method == 1) then
            call fft_navier(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, &
                            Bx, By, Tx, Ty, Tx0, Ty0, Fx, Fy, step)
            call cpu_time(t4)
            call fft_poisson(Up, Vp, Phi)
            call fft_march(Up, Vp, U, V, Phi, P)
        endif

        if (method == 2) then
            call ibm_preliminary(U, V, P, Ax, Ay, Ax0, Ay0, Bx, By, Tx, Ty, Tx0, Ty0, Fx, Fy, &
                                 Ua, Va, Fxc, Fyc, Up, Vp, step)
            call ibm_Helmholtz(Up, Vp, X, Y, Xc, Yc, Uc, Vc, Fxc, Fyc, fxint, fyint, step)
            call ibm_predict(Ua, Va, fxint, fyint, Bx, By, Up, Vp)
            call cpu_time(t4)
            call fft_poisson(Up, Vp, Phi)
            call fft_march(Up, Vp, U, V, Phi, P)
        endif

        if (method == 3) then
            call fft_navier_kabe(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, &
                                 Bx, By, Tx, Ty, Tx0, Ty0, Fx, Fy, step)
            call cpu_time(t4)
            call fft_poisson_kabe(Up, Vp, Phi)
            call fft_march_kabe(Up, Vp, U, V, Phi, P, C)
        endif

        call cpu_time(t5)

        call logging(U, V, C, step, t_start)
        if (mod(step, Gstep)==0) call get_data(U, V, C, step)
        if (mod(step, Estep)==0) call energy_single(U, V, step)
        ! if (mod(step, output_step)==0) call output(U, V, P, step)
        if (sum(U**2)*0.0d0 /= 0.0d0) stop 'NaN value'  ! NaNの判定
        t12 = t12 + t2-t1
        t23 = t23 + t3-t2
        t34 = t34 + t4-t3
        t45 = t45 + t5-t4

        ! write(*, '(12e12.4)') U(1:12, 1)
        if (mod(step, Dstep)==0) then
            open(11, file=trim(dir)//'debug.d')
            if (method == 3) then
                do j = NA, NB+1
                    write(11, '(3e12.4)') (j-0.5d0)*dY_C, sum(U(1:NX, j))/NX*U_C, sum(C(1, 1:NX, j)+C(3, 1:NX, j))/NX
                enddo
            else
                do j = 1, NY
                    write(11, '(3e12.4)') (j-0.5d0)*dY_C, sum(U(1:NX, j))/NX*U_C, sum(C(1, 1:NX, j)+C(3, 1:NX, j))/NX
                enddo
            endif
            close(11)
        endif
        
        if (mod(step, Gstep)==0 .and. flow_type == 4) call taylor_debug(U, V, step)
        if (mod(step, Dstep)==0 .and. flow_type == 2) call koide_theory(U, V, C)
    enddo


    call fft_finalize
    call cpu_time(t_end)
    total_time = t_end-t_start
    others_time = total_time - (t12+t23+t34+t45)
    write(*, '(a9, F10.3, a3)') 'Total   :', total_time, '[s]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'C+-     :', t12, '[s]', t12/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Newton  :', t23, '[s]', t23/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Navier  :', t34, '[s]', t34/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Poisson :', t45, '[s]', t45/total_time*100, '[%]'
    write(*, '(a9, F10.3, a3, F8.3, a3)') 'Others  :', others_time, '[s]', others_time/total_time*100, '[%]'

end program main

