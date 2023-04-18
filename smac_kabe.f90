module smac
    implicit none
    ! gfortran smac.f90 -I/$HOME/local/include -L$HOME/local/lib -lfftw3 && ./a.out
    integer, parameter :: Nstep = 10000
    integer, parameter :: Gstep = 10000 ! データを取得する間隔
    integer, parameter :: method = 1  ! 0:陽解法、1:FFT、2:IBM，3:時間平均エネルギーカスケード
    integer, parameter :: pattern = 1  ! 0:周期境界条件、1:z壁面壁（速度0）
    integer, parameter :: NX = 16, NY = NX, NZ = NX  ! 70までifort ok
    real(8), parameter :: PI = acos(-1.0d0)
    real(8), parameter :: Xmax = 2*PI, Ymax = 2*PI, Zmax = 2*PI/NX*NZ
    real(8), parameter :: dX_C = Xmax/NX, dY_C = Ymax/NY, dZ_C = Zmax/NZ
    real(8), parameter :: dt_C = 0.01d0
    ! real(8), parameter :: nu = 1.77245d0  ! Re=1
    ! real(8), parameter :: U_C = 1.0d0/(2.0d0*nu)
    ! real(8), parameter :: Re = Xmax*U_C/nu
    real(8), parameter :: Re = 5.0d0  ! 無次元量
    real(8), parameter :: nu = sqrt(Xmax/(2.0d0*Re))
    ! real(8), parameter :: U_C = 0.5d0/nu  ! テイラーグリーン
    real(8), parameter :: Uwall =  0.0d0  ! 壁面の規格化速度
    ! real(8), parameter :: U_C = Uwall  ! クエット
    real(8), parameter :: U_C = Zmax**2*Re/8.0d0  ! ポアズイユ
    ! real(8), parameter :: A = 1.0d0
    real(8), parameter :: dX = dX_C/Xmax, dY = dY_C/Ymax, dZ = dZ_C/Zmax
    real(8), parameter :: dt = dt_C*U_C/Xmax
    real(8), parameter :: SOR = 2.0d0 / (1 + sin(PI*dX))
    real(8), parameter :: eps = 1.0d-5
    integer, parameter :: itrmax = 10000
    ! real(8), parameter :: f0 = 1.0d0
    real(8), parameter :: f0 = Xmax/(U_C)**2  ! 外力を1に固定し，無次元化するときの定数
    ! IBM用パラメータ
    real(8), parameter :: D_C = Xmax/4  ! 円柱直径
    integer, parameter :: NC = 4  ! 円柱の個数
    integer, parameter :: NS = 1000  ! 反復回数
    integer, parameter :: NL = nint(PI*D_C/dX_C)  ! 円周方向の分割数
    real(8), parameter :: dV = PI*D_C/Xmax/NL * dX * dX
    ! エネルギーカスケード
    integer, parameter :: NE = 1000  ! 時間平均反復回数
    ! 入出力ファイル
    ! character(*), parameter :: dir = 'N128_dt005_Re500_f01/'  ! なければ自動で作成される
    character(*), parameter :: dir = './data/' 
    integer, parameter :: input_step = 0  ! 0以外で初期条件をファイルから読み込む
    integer, parameter :: output_step = 100000000  ! 配列を保存する間隔

contains
    subroutine init(U, V, W, P, Phi, Fx, Fy, Fz)  ! はじめに実行
        real(8), intent(out) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: P(0:NX+1, 0:NY+1, 0:NZ+1), Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
        integer i, j, k
        real(8), parameter :: large_K = 2.0d0

        if (method == 0) then
            write(*, '(a, F8.3, F8.3, F8.3, a8)') 'Von Neumann:', 1.0d0/Re*dt/dX**2, 1.0d0/Re*dt/dY**2, 1.0d0/Re*dt/dZ**2, '< 0.167'
        endif
        write(*, '(a, F8.3, F8.3, F8.3, a8)') 'CFL:', 1.0d0*dt/dX, 1.0d0*dt/dY, 1.0d0*dt/dZ,'< 0.167'
        write(*, '(a, F8.3)') 'U_C =', U_C
        write(*, '(a, F8.3)') 'Re  =', Re
        write(*, '(a, E12.4)') 'dt  =', dt

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
        Fx(:, :, :) = 0.0d0
        Fy(:, :, :) = 0.0d0
        Fz(:, :, :) = 0.0d0
        P(:, :, :) = 0.0d0
        Phi(:, :, :) = 0.0d0
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    ! U(i, j, k) = -sin(i*dX_C) * cos((j-0.5d0)*dY_C)
                    ! U(i, j, k) = U(i, j, k) + 1.0d0
                    ! V(i, j, k) = cos((i-0.5d0)*dX_C) * sin(j*dY_C)
                    ! V(i, j, k) = V(i, j, k) + 1.0d0
                    ! W(i, j, k) = W(i, j, k) + 1.0d0
                    ! W(i, j, k) = W(i, j, k) + sin((i-0.5d0)*dX_C)

                    ! Fx(i, j, k) = -sin(i*dX_C) * cos((j-0.5d0)*dY_C) * cos((k-0.5d0)*dZ_C)  ! 荒木さん
                    ! Fy(i, j, k) = cos((i-0.5d0)*dX_C) * sin(j*dY_C) * cos((k-0.5d0)*dZ_C)
                    ! Fx(i, j, k) = -f0 * sin(i*dX_C) * cos((j-0.5d0)*dY_C)  ! 増田さん, 安房井さん
                    ! Fy(i, j, k) = f0 * cos((i-0.5d0)*dX_C) * sin(j*dY_C)
                    ! Fx(i, j, k) = -sin(large_K*(j-0.5d0)*dY_C)  ! 小井手さん
                    ! Fy(i, j, k) = sin(large_K*(i-0.5d0)*dX_C)
                    ! Fx(i, j, k) = f0 * A
                enddo
            enddo
        enddo
        if (pattern == 0) then
            call PBM(U)
            call PBM(V)
            call PBM(W)
            call PBM(P)
        else if (pattern == 1) then
            call U_kabe(U)
            call V_kabe(V)
            call W_kabe(W)
            call P_kabe(P)
        endif
    end subroutine init

    subroutine PBM(A)
        real(8), intent(inout) :: A(0:NX+1, 0:NY+1, 0:NZ+1)
        A(0, :, :) = A(NX, :, :)
        A(NX+1, :, :) = A(1, :, :)
        A(:, 0, :) = A(:, NY, :)
        A(:, NY+1, :) = A(:, 1, :)
        A(:, :, 0) = A(:, :, NZ)
        A(:, :, NZ+1) = A(:, :, 1)
    end subroutine PBM

    subroutine U_kabe(U)
        real(8), intent(inout) :: U(0:NX+1, 0:NY+1, 0:NZ+1)
        U(0, :, :) = U(NX, :, :)  ! x, y周期境界条件
        U(NX+1, :, :) = U(1, :, :)
        U(:, 0, :) = U(:, NY, :)
        U(:, NY+1, :) = U(:, 1, :)
        U(:, :, 0) = - U(:, :, 1)  ! z壁
        U(:, :, NZ+1) =  2*Uwall - U(:, :, NZ)
    end subroutine U_kabe

    subroutine V_kabe(V)
        real(8), intent(inout) :: V(0:NX+1, 0:NY+1, 0:NZ+1)
        V(0, :, :) = V(NX, :, :)  ! x, y周期境界条件
        V(NX+1, :, :) = V(1, :, :)
        V(:, 0, :) = V(:, NY, :)
        V(:, NY+1, :) = V(:, 1, :)
        V(:, :, 0) = - V(:, :, 1)  ! z壁
        V(:, :, NZ+1) = - V(:, :, NZ)
    end subroutine V_kabe

    subroutine W_kabe(W)
        real(8), intent(inout) :: W(0:NX+1, 0:NY+1, 0:NZ+1)
        W(0, :, :) = W(NX, :, :)  ! x, y周期境界条件
        W(NX+1, :, :) = W(1, :, :)
        W(:, 0, :) = W(:, NY, :)
        W(:, NY+1, :) = W(:, 1, :)
        W(:, :, 0) = 0.0d0  ! z壁
        W(:, :, NZ) = 0.0d0
        W(:, :, NZ+1) = 0.0d0
    end subroutine W_kabe

    subroutine P_kabe(P)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        P(0, :, :) = P(NX, :, :)  ! x, y周期境界条件
        P(NX+1, :, :) = P(1, :, :)
        P(:, 0, :) = P(:, NY, :)
        P(:, NY+1, :) = P(:, 1, :)
        P(:, :, 0) = P(:, :, 1)  ! z壁
        P(:, :, NZ+1) = P(:, :, NZ)
    end subroutine P_kabe

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

    subroutine viscous(U, V, W, Bx, By, Bz)  ! ニュー一定
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
        Bx(:, :, :) = 1.0d0/Re*Bx(:, :, :)
        By(:, :, :) = 1.0d0/Re*By(:, :, :)
        Bz(:, :, :) = 1.0d0/Re*Bz(:, :, :)
    end subroutine viscous

    subroutine navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Bx, By, Bz, Ax0, Ay0, Az0, Bx0, By0, Bz0, Fx, Fy, Fz, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Bx0(1:NX, 1:NY, 1:NZ), By0(1:NX, 1:NY, 1:NZ), Bz0(1:NX, 1:NY, 1:NZ)
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
        endif

        ! NS方程式で速度の予測
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Up(i, j, k) = U(i, j, k) - dt*(-P(i, j, k)+P(i+1, j, k))/dX &
                                 +dt*(3.0d0*(Ax(i, j, k) + Bx(i, j, k)) - (Ax0(i, j, k) + Bx0(i, j, k)))/2.0d0 &
                                 +dt*Fx(i, j, k)
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Vp(i, j, k) = V(i, j, k) - dt*(-P(i, j, k)+P(i, j+1, k))/dY &
                                 +dt*(3.0d0*(Ay(i, j, k) + By(i, j, k)) - (Ay0(i, j, k) + By0(i, j, k)))/2.0d0 &
                                 +dt*Fy(i, j, k)
                enddo
            enddo
        enddo
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Wp(i, j, k) = W(i, j, k) - dt*(-P(i, j, k)+P(i, j, k+1))/dZ &
                                 +dt*(3.0d0*(Az(i, j, k) + Bz(i, j, k)) - (Az0(i, j, k) + Bz0(i, j, k)))/2.0d0 &
                                 +dt*Fz(i, j, k)
                enddo
            enddo
        enddo
        ! Wp(:, :, :) = 0.0d0

        ! Phiを求めるときに0番目の値も必要
        Up(0, :, :) = Up(NX, :, :)
        Vp(:, 0, :) = Vp(:, NY, :)
        Wp(:, :, 0) = Wp(:, :, NZ)

        ! n-1ステップ目の保存
        Ax0(:, :, :) = Ax(:, :, :)
        Ay0(:, :, :) = Ay(:, :, :)
        Az0(:, :, :) = Az(:, :, :)
        Bx0(:, :, :) = Bx(:, :, :)
        By0(:, :, :) = By(:, :, :)
        Bz0(:, :, :) = Bz(:, :, :)
    end subroutine navier

    subroutine poisson(Up, Vp, Wp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k, l, itr
        real(8) BXM, BXP, BYM, BYP, BZM, BZP, B0
        real(8) er, er0
        real(8) E(1:NX, 1:NY, 1:NZ), Q(1:NX, 1:NY, 1:NZ)
        real(8) Phi0(1:NX, 1:NY, 1:NZ)

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

        ! ポアソン方程式で用いる定数の計算
        BXM = 1.0d0 / dX**2
        BXP = 1.0d0 / dX**2
        BYM = 1.0d0 / dY**2
        BYP = 1.0d0 / dY**2
        BZM = 1.0d0 / dZ**2
        BZP = 1.0d0 / dZ**2
        B0 = BXM + BXP + BYM + BYP + BZM + BZP
        Phi(:, :, :) = 0.0d0
        ! SOR法
        do itr = 1, itrmax
            do k = 1, NZ
                do j = 1, NY
                    do l = 1, 2
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
            er0 = sum(Q**2)
            if (er0 == 0.0d0) then  ! Qが全て0のとき計算できない
                er0 = sum(Q**2) + 1.0d-16
            endif
            if (sqrt(er / er0) < eps) then
                ! write(*, *) 'step =', step, 'itr =', itr, 'err =', sqrt(er / er0)
                exit
            endif
        enddo
        ! Phiの平均を0にするため
        B0 = sum(Phi(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
        Phi(:, :, :) = Phi(:, :, :) - B0

        ! デバッグ用
        er = 0.0d0
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX/2+1
                    ! write(*, '(3I3, 4e12.4)') i, j, k, Phi(i, j, k), Phi0(i, j, k), Phi(i, j, k)-Phi0(i, j, k)
                    er = er + (Phi(i, j, k)-Phi0(i, j, k))**2
                enddo
            enddo
        enddo

        ! write(*, '(12e12.4)') sqrt(er/(NX*NY*NZ))
        l = 1  ! 使ってません対策
        ! write(*, *) 'step =', step, 'itr =', itr, 'er', er, 'er0', er0 
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
        ! W(:, :, :) = 0.0d0
        do k = 0, NZ+1
            do j = 0, NY+1
                do i = 0, NX+1
                    P(i, j, k) = P(i, j, k) + Phi(i, j, k)
                enddo
            enddo
        enddo

        call PBM(U)
        call PBM(V)
        call PBM(W)
        call PBM(P)
    end subroutine march


    ! subroutine get_data(U, V, W, step)
    !     real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
    !     integer, intent(in) :: step
    !     integer i, j, k, l
    !     real(8) D(3, 3), Omega(3), S(3), Qti
    !     character(8) str
    !     write(str, '(I8.8)') step  ! 数値を文字列に変換

    !     ! open(10, file='data/y_'//str//'.d')
    !     open(10, file=dir//str//'.d')

    !     do k = 1, NZ
    !         do j = 1, NY
    !             do i = 1, NX
    !                 D(1, 1) = (U(i, j, k) - U(i-1, j, k))/dX
    !                 D(1, 2) = (U(i, j+1, k) - U(i, j-1, k) + U(i-1, j+1, k) - U(i-1, j-1, k))/(4*dY)
    !                 D(1, 3) = (U(i, j, k+1) - U(i, j, k-1) + U(i-1, j, k+1) - U(i-1, j, k-1))/(4*dZ)
    !                 D(2, 1) = (V(i+1, j, k) - V(i-1, j, k) + V(i+1, j-1, k) - V(i-1, j-1, k))/(4*dX)
    !                 D(2, 2) = (V(i, j, k) - V(i, j-1, k))/dY
    !                 D(2, 3) = (V(i, j, k+1) - V(i, j, k-1) + V(i, j-1, k+1) - V(i, j-1, k-1))/(4*dZ)
    !                 D(3, 1) = (W(i+1, j, k) - W(i-1, j, k) + W(i+1, j, k-1) - W(i-1, j, k-1))/(4*dX)
    !                 D(3, 2) = (W(i, j+1, k) - W(i, j-1, k) + W(i, j+1, k-1) - W(i, j-1, k-1))/(4*dY)
    !                 D(3, 1) = (W(i, j, k) - W(i, j, k-1))/dZ
    !                 D(:, :) = D(:, :)*U_C

    !                 Omega(1) = (D(3, 2) - D(2, 3))/2
    !                 Omega(2) = (D(1, 3) - D(3, 1))/2
    !                 Omega(3) = (D(2, 1) - D(1, 2))/2
    !                 S(1) = (D(3, 2) + D(2, 3))/2
    !                 S(2) = (D(1, 3) + D(3, 1))/2
    !                 S(3) = (D(2, 1) + D(1, 2))/2
    !                 Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
    !                 write(10, '(11e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
    !                                     (U(i-1, j, k)+U(i, j, k))/2*U_C, &
    !                                     (V(i, j-1, k)+V(i, j, k))/2*U_C, &
    !                                     (W(i, j, k-1)+W(i, j, k))/2*U_C, &
    !                                     Omega(3), sum(Omega**2)/2, S(3), sum(S**2)/2, Qti
    !             enddo
    !         enddo
    !     enddo
    !     close(10)
    !     l = 1  ! 使ってません対策
    ! end subroutine get_data

    subroutine get_data(U, V, W, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k, l
        real(8) D(3, 3), Omega(3), S(3), Qti
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換
        ! open(10, file='data/y_'//str//'.d')
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

                Omega(1) = (D(3, 2) - D(2, 3))/2
                Omega(2) = (D(1, 3) - D(3, 1))/2
                Omega(3) = (D(2, 1) - D(1, 2))/2
                write(10, '(11e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
                                        (U(i-1, j, k)+U(i, j, k))/2*U_C, &
                                        (V(i, j-1, k)+V(i, j, k))/2*U_C, &
                                        (W(i, j, k-1)+W(i, j, k))/2*U_C, &
                                        Omega(2)
            enddo
        enddo
        close(10)
        l = 1  ! 使ってません対策
    end subroutine get_data

    
    subroutine input(U, V, W, P)  ! initを実行した後に書く
        real(8), intent(out) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        character(8) str
        write(str, '(I8.8)') input_step  ! 数値を文字列に変換
        open(10, file=dir//str//'.d')

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    read(10, '(10e12.4)') U(i, j, k), V(i, j, k), W(i, j, k), P(i, j, k)
                enddo
            enddo
        enddo
        close(10)
        call PBM(U)
        call PBM(V)
        call PBM(W)
        call PBM(P)
    end subroutine input

    subroutine output(U, V, W, P, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換
        open(10, file=dir//str//'.d')

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    write(10, '(10e12.4)') U(i, j, k), V(i, j, k), W(i, j, k), P(i, j, k)
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
    integer(8) plan1, plan2, plan3, plan4
    real(8) Re1(1:NX, 1:NY, 1:NZ), Re2(1:NX, 1:NY)
    complex(8) Im1(1:NX/2+1, 1:NY, 1:NZ), Im2(1:NX/2+1, 1:NY)
contains
    subroutine fft_init
        call dfftw_plan_dft_r2c_3d(plan1, NX, NY, NZ, Re1, Im1, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_3d(plan2, NX, NY, NZ, Im1, Re1, FFTW_ESTIMATE)
        call dfftw_plan_dft_r2c_2d(plan3, NX, NY, Re2, Im2, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(plan4, NX, NY, Im2, Re2, FFTW_ESTIMATE)
    end subroutine fft_init

    subroutine fftr2c_3d(input, output)
        real(8), intent(in) :: input(1:NX, 1:NY, 1:NZ)
        complex(8), intent(out) :: output(1:NX/2+1, 1:NY, 1:NZ)
        Re1(:, :, :) = input(:, :, :)
        call dfftw_execute(plan1, Re1, Im1)
        output(:, :, :) = Im1(:, :, :)
    end subroutine fftr2c_3d

    subroutine fftc2r_3d(input, output)
        complex(8), intent(in) :: input(1:NX/2+1, 1:NY, 1:NZ)
        real(8), intent(out) :: output(1:NX, 1:NY, 1:NZ)
        Im1(:, :, :) = input(:, :, :)
        call dfftw_execute(plan2, Im1, Re1)
        output(:, :, :) = Re1(:, :, :)
    end subroutine fftc2r_3d

    subroutine fftr2c_2d(input, output)
        real(8), intent(in) :: input(1:NX, 1:NY)
        complex(8), intent(out) :: output(1:NX/2+1, 1:NY)
        Re2(:, :) = input(:, :)
        call dfftw_execute(plan3, Re2, Im2)
        output(:, :) = Im2(:, :)
    end subroutine fftr2c_2d

    subroutine fftc2r_2d(input, output)
        complex(8), intent(in) :: input(1:NX/2+1, 1:NY)
        real(8), intent(out) :: output(1:NX, 1:NY)
        Im2(:, :) = input(:, :)
        call dfftw_execute(plan4, Im2, Re2)
        output(:, :) = Re2(:, :)
    end subroutine fftc2r_2d

    subroutine fft_finalize
        call dfftw_destroy_plan(plan1)
        call dfftw_destroy_plan(plan2)
        call dfftw_destroy_plan(plan3)
        call dfftw_destroy_plan(plan4)
    end subroutine fft_finalize

    
    subroutine fft_solve(Q, LHS, Phi)
        real(8), intent(in) :: Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)
        real(8), intent(out) :: Phi(1:NX, 1:NY, 1:NZ)
        complex(8) Q_hat(1:NX/2+1, 1:NY, 1:NZ), Phi_hat(1:NX/2+1, 1:NY, 1:NZ)
        integer i, j, k
        integer(8) plan

        ! 右辺をフーリエ変換
        call fftr2c_3d(Q, Q_hat)
        
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

        ! 左辺を逆フーリエ変換
        call fftc2r_3d(Phi_hat, Phi)

    end subroutine fft_solve

    subroutine fft_navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Fx, Fy, Fz, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
        integer, intent(in) :: step
        integer i, j, k
        real(8) Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :, :) = Ax(:, :, :)
            Ay0(:, :, :) = Ay(:, :, :)
            Az0(:, :, :) = Az(:, :, :)
        endif

        ! 左辺の定数部分の計算
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    LHS(i, j, k) = 1 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX*2*PI))/dX**2 &
                                                  + 2*(1-cos((j-1)*dY*2*PI))/dY**2 &
                                                  + 2*(1-cos((k-1)*dZ*2*PI))/dZ**2)
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
    end subroutine fft_navier

    subroutine fft_poisson(Up, Vp, Wp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        real(8) Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    LHS(i, j, k) = -(2*(1-cos((i-1)*dX*2*PI))/dX**2 + &
                                     2*(1-cos((j-1)*dY*2*PI))/dY**2 + &
                                     2*(1-cos((k-1)*dZ*2*PI))/dZ**2)
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

    subroutine fft_navier_kabe(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Fx, Fy, Fz, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
        integer, intent(in) :: step
        integer i, j, k
        real(8) Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)
        complex(8) Q_hat(1:NX/2+1, 1:NY, 1:NZ), E_hat(1:NX/2+1, 1:NY, 1:NZ)
        complex(8) Up_hat(1:NX/2+1, 1:NY, 0:NZ+1), Vp_hat(1:NX/2+1, 1:NY, 0:NZ+1), Wp_hat(1:NX/2+1, 1:NY, 0:NZ+1)
        real(8) BZM, BZP, B0, er, er0
        integer itr

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :, :) = Ax(:, :, :)
            Ay0(:, :, :) = Ay(:, :, :)
            Az0(:, :, :) = Az(:, :, :)
        endif


        ! 速度場Upを予測
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = U(i, j, k) &
                                -dt*(-P(i, j, k)+P(i+1, j, k))/dX &
                                +dt*(3.0d0*Ax(i, j, k) - Ax0(i, j, k))/2.0d0 &
                                +dt*Bx(i, j, k)/2.0d0 &
                                +dt*Fx(i, j, k)
                enddo
            enddo
        enddo

        ! 右辺をフーリエ変換
        do k = 1, NZ
            call fftr2c_2d(Q(:, :, k), Q_hat(:, :, k))
        enddo
        Q_hat(:, :, :) = Q_hat(:, :, :)/(NX*NY)

        ! SOR法で3重対角行列を解く
        Up_hat(:, :, :) = 0.0d0
        do j = 1, NY
            do i = 1, NX/2+1
                BZM = - dt/Re/2.0d0 / dZ**2
                BZP = - dt/Re/2.0d0 / dZ**2
                B0 = -(1 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX*2*PI))/dX**2 + 2*(1-cos((j-1)*dY*2*PI))/dY**2 + 2.0d0/dZ**2))
                do itr = 1, itrmax
                    do k = 1, NZ
                        E_hat(i, j, k) = BZM*Up_hat(i, j, k-1) - B0*Up_hat(i, j, k) + BZP*Up_hat(i, j, k+1) - Q_hat(i, j, k)
                        Up_hat(i, j, k) = Up_hat(i, j, k) + SOR*E_hat(i, j, k) / B0
                        ! Up_hat(i, j, k) = Up_hat(i, j, k) + E_hat(i, j, k) / B0
                    enddo
                    ! 周期境界条件
                    if (pattern == 0) Up_hat(i, j, 0) = Up_hat(i, j, NZ)
                    if (pattern == 0) Up_hat(i, j, NZ+1) = Up_hat(i, j, 1)
                    ! z壁面条件
                    if (pattern == 1) Up_hat(i, j, 0) = -Up_hat(i, j, 1)
                    if (pattern == 1) then
                        if (i == 1 .and. j == 1) then
                            Up_hat(i, j, NZ+1) = 2*Uwall - Up_hat(i, j, NZ)  ! デルタ関数
                            ! Up_hat(i, j, NZ+1) = - Up_hat(i, j, NZ)
                        else 
                            Up_hat(i, j, NZ+1) = - Up_hat(i, j, NZ)
                        endif
                    endif

                    er = sum(abs(E_hat(i, j, :))**2)
                    er0 = sum(abs(Q_hat(i, j, :))**2) + 1.0d-16
                    if (sqrt(er / er0) < eps) then
                        ! write(*, '(2I4, a, I6, a, E12.4)') i, j, '  itr =', itr, '  err =', sqrt(er / er0)
                        exit
                    endif
                enddo
            enddo
        enddo
        ! 左辺を逆フーリエ変換
        do k = 0, NZ+1  ! 範囲に注意
            call fftc2r_2d(Up_hat(:, :, k), Up(1:NX, 1:NY, k))
        enddo

        
        ! 同様にVpについても求める
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = V(i, j, k) &
                                -dt*(-P(i, j, k)+P(i, j+1, k))/dY &
                                +dt*(3.0d0*Ay(i, j, k) - Ay0(i, j, k))/2.0d0 &
                                +dt*By(i, j, k)/2.0d0 &
                                +dt*Fy(i, j, k)
                enddo
            enddo
        enddo
        ! 右辺をフーリエ変換
        do k = 1, NZ
            call fftr2c_2d(Q(:, :, k), Q_hat(:, :, k))
        enddo
        Q_hat(:, :, :) = Q_hat(:, :, :)/(NX*NY)
        ! SOR法で3重対角行列を解く
        Vp_hat(:, :, :) = 0.0d0
        do j = 1, NY
            do i = 1, NX/2+1
                BZM = - dt/Re/2.0d0 / dZ**2
                BZP = - dt/Re/2.0d0 / dZ**2
                B0 = -(1 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX*2*PI))/dX**2 + 2*(1-cos((j-1)*dY*2*PI))/dY**2 + 2.0d0/dZ**2))
                do itr = 1, itrmax
                    do k = 1, NZ
                        E_hat(i, j, k) = BZM*Vp_hat(i, j, k-1) - B0*Vp_hat(i, j, k) + BZP*Vp_hat(i, j, k+1) - Q_hat(i, j, k)
                        Vp_hat(i, j, k) = Vp_hat(i, j, k) + SOR*E_hat(i, j, k) / B0
                        ! Vp_hat(i, j, k) = Vp_hat(i, j, k) + E_hat(i, j, k) / B0
                    enddo
                    ! 周期境界条件
                    if (pattern == 0) Vp_hat(i, j, 0) = Vp_hat(i, j, NZ)
                    if (pattern == 0) Vp_hat(i, j, NZ+1) = Vp_hat(i, j, 1)
                    ! z壁面条件
                    if (pattern == 1) Vp_hat(i, j, 0) = -Vp_hat(i, j, 1)
                    if (pattern == 1) Vp_hat(i, j, NZ+1) = -Vp_hat(i, j, NZ)

                    er = sum(abs(E_hat(i, j, :))**2)
                    er0 = sum(abs(Q_hat(i, j, :))**2) + 1.0d-16
                    if (sqrt(er / er0) < eps) then
                        ! write(*, '(2I4, a, I6, a, E12.4)') i, j, '  itr =', itr, '  err =', sqrt(er / er0)
                        exit
                    endif
                enddo
            enddo
        enddo
        ! 左辺を逆フーリエ変換
        do k = 0, NZ+1  ! 範囲に注意
            call fftc2r_2d(Vp_hat(:, :, k), Vp(1:NX, 1:NY, k))
        enddo


        ! 同様にWpについても求める
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = W(i, j, k) &
                                -dt*(-P(i, j, k)+P(i, j, k+1))/dZ &
                                +dt*(3.0d0*Az(i, j, k) - Az0(i, j, k))/2.0d0 &
                                +dt*Bz(i, j, k)/2.0d0 &
                                +dt*Fz(i, j, k)
                enddo
            enddo
        enddo
        ! 右辺をフーリエ変換
        do k = 1, NZ
            call fftr2c_2d(Q(:, :, k), Q_hat(:, :, k))
        enddo
        Q_hat(:, :, :) = Q_hat(:, :, :)/(NX*NY)
        ! SOR法で3重対角行列を解く
        Wp_hat(:, :, :) = 0.0d0
        do j = 1, NY
            do i = 1, NX/2+1
                BZM = - dt/Re/2.0d0 / dZ**2
                BZP = - dt/Re/2.0d0 / dZ**2
                B0 = -(1 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX*2*PI))/dX**2 + 2*(1-cos((j-1)*dY*2*PI))/dY**2 + 2.0d0/dZ**2))
                do itr = 1, itrmax
                    ! 周期境界条件
                    if (pattern == 0) then
                        do k = 1, NZ
                            E_hat(i, j, k) = BZM*Wp_hat(i, j, k-1) - B0*Wp_hat(i, j, k) + BZP*Wp_hat(i, j, k+1) - Q_hat(i, j, k)
                            Wp_hat(i, j, k) = Wp_hat(i, j, k) + SOR*E_hat(i, j, k) / B0
                            ! Wp_hat(i, j, k) = Wp_hat(i, j, k) + E_hat(i, j, k) / B0
                        enddo
                        Wp_hat(i, j, 0) = Wp_hat(i, j, NZ)
                        Wp_hat(i, j, NZ+1) = Wp_hat(i, j, 1)

                    ! z壁面条件
                    else if (pattern == 1) then
                        Wp_hat(i, j, 0) = 0.0d0
                        Wp_hat(i, j, NZ) = 0.0d0
                        Wp_hat(i, j, NZ+1) = 0.0d0
                        E_hat(i, j, NZ) = 0.0d0
                        Q_hat(i, j, NZ) = 0.0d0
                        do k = 1, NZ-1
                            E_hat(i, j, k) = BZM*Wp_hat(i, j, k-1) - B0*Wp_hat(i, j, k) + BZP*Wp_hat(i, j, k+1) - Q_hat(i, j, k)
                            Wp_hat(i, j, k) = Wp_hat(i, j, k) + SOR*E_hat(i, j, k) / B0
                            ! Wp_hat(i, j, k) = Wp_hat(i, j, k) + E_hat(i, j, k) / B0
                        enddo
                    endif

                    er = sum(abs(E_hat(i, j, :))**2)
                    er0 = sum(abs(Q_hat(i, j, :))**2) + 1.0d-16
                    if (sqrt(er / er0) < eps) then
                        ! write(*, '(2I4, a, I6, a, E12.4)') i, j, '  itr =', itr, '  err =', sqrt(er / er0)
                        exit
                    endif
                enddo
            enddo
        enddo
        ! 左辺を逆フーリエ変換
        do k = 0, NZ+1  ! 範囲に注意
            call fftc2r_2d(Wp_hat(:, :, k), Wp(1:NX, 1:NY, k))
        enddo


        ! 周期境界条件
        if (pattern == 0) then
            call PBM(Up)
            call PBM(Vp)
            call PBM(Wp)
        ! z壁面条件
        else if (pattern == 1) then 
            call U_kabe(Up)
            call V_kabe(Vp)
            call W_kabe(Wp)
        endif

        ! n-1ステップ目の保存
        Ax0(:, :, :) = Ax(:, :, :)
        Ay0(:, :, :) = Ay(:, :, :)
        Az0(:, :, :) = Az(:, :, :)
    end subroutine fft_navier_kabe

    subroutine fft_poisson_kabe(Up, Vp, Wp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(out) :: Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        real(8) Q(1:NX, 1:NY, 1:NZ), LHS(1:NX, 1:NY, 1:NZ)
        real(8) Re(1:NX, 1:NY)
        complex(8) Im(1:NX/2+1, 1:NY)
        complex(8) Q_hat(1:NX/2+1, 1:NY, 1:NZ), Phi_hat(1:NX/2+1, 1:NY, 0:NZ+1)
        real(8) BZM, BZP, B0, er, er0
        integer itr
        complex(8) E_hat(1:NX/2+1, 1:NY, 1:NZ)


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

        ! 右辺をフーリエ変換
        do k = 1, NZ
            call fftr2c_2d(Q(:, :, k), Q_hat(:, :, k))
        enddo
        Q_hat(:, :, :) = Q_hat(:, :, :)/(NX*NY)
        ! write(*, *) sum(abs(Q_hat(:, :, :))**2)

        ! do k = 1, NZ
        !     call fftc2r_2d(Q_hat(:, :, k), LHS(:, :, k))
        ! enddo
        ! write(*, *) sum((Q(:, :, :)-LHS(:, :, :))**2)

        ! SOR法で3重対角行列を解く
        Phi_hat(:, :, :) = 0.0d0
        do j = 1, NY
            do i = 1, NX/2+1
                BZM = 1.0d0 / dZ**2
                BZP = 1.0d0 / dZ**2
                B0 = 2*(1-cos((i-1)*dX*2*PI))/dX**2 + 2*(1-cos((j-1)*dY*2*PI))/dY**2 + 2.0d0/dZ**2
                do itr = 1, itrmax
                    do k = 1, NZ
                        E_hat(i, j, k) = BZM*Phi_hat(i, j, k-1) - B0*Phi_hat(i, j, k) + BZP*Phi_hat(i, j, k+1) - Q_hat(i, j, k)
                        Phi_hat(i, j, k) = Phi_hat(i, j, k) + SOR*E_hat(i, j, k) / B0
                        ! Phi_hat(i, j, k) = Phi_hat(i, j, k) + E_hat(i, j, k) / B0
                    enddo
                    ! 周期境界条件
                    if (pattern == 0) Phi_hat(i, j, 0) = Phi_hat(i, j, NZ)
                    if (pattern == 0) Phi_hat(i, j, NZ+1) = Phi_hat(i, j, 1)
                    ! z壁面条件
                    if (pattern == 1) Phi_hat(i, j, 0) = Phi_hat(i, j, 1)
                    if (pattern == 1) Phi_hat(i, j, NZ+1) = Phi_hat(i, j, NZ)

                    er = sum(abs(E_hat(i, j, :))**2)
                    er0 = sum(abs(Q_hat(i, j, :))**2) + 1.0d-16
                    if (sqrt(er / er0) < eps) then
                        ! write(*, '(2I4, a, I6, a, E12.4)') i, j, '  itr =', itr, '  err =', sqrt(er / er0)
                        exit
                    endif
                enddo
            enddo
        enddo

        ! 左辺を逆フーリエ変換
        do k = 0, NZ+1  ! 範囲に注意
            call fftc2r_2d(Phi_hat(:, :, k), Phi(1:NX, 1:NY, k))
        enddo

        ! Phiの平均を0にするため
        B0 = sum(Phi(1:NX, 1:NY, 1:NZ))/(NX*NY*NZ)
        Phi(:, :, :) = Phi(:, :, :) - B0

        if (pattern == 0) call PBM(Phi)
        if (pattern == 1) call P_kabe(Phi)

    end subroutine fft_poisson_kabe

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
        ! W(:, :, :) = 0.0d0
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    P(i, j, k) = P(i, j, k) + Phi(i, j, k) &
                                -dt/Re/2.0d0*((Phi(i-1, j, k)-2*Phi(i, j, k)+Phi(i+1, j, k))/dX**2 &
                                             +(Phi(i, j-1, k)-2*Phi(i, j, k)+Phi(i, j+1, k))/dY**2 &
                                             +(Phi(i, j, k-1)-2*Phi(i, j, k)+Phi(i, j, k+1))/dZ**2)
                enddo
            enddo
        enddo
        if (pattern == 0) then
            call PBM(U)
            call PBM(V)
            call PBM(W)
            call PBM(P)
        else if (pattern == 1) then
            call U_kabe(U)
            call V_kabe(V)
            call W_kabe(W)
            call P_kabe(P)
        endif
    end subroutine fft_march


    subroutine fft_forward(Q, Q_hat)  ! FFT用関数
        real(8), intent(in) :: Q(1:NX, 1:NY, 1:NZ)
        complex(8), intent(out) :: Q_hat(1:NX/2+1, 1:NY, 1:NZ)
        integer(8) plan

        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, Q, Q_hat, FFTW_ESTIMATE)
        call dfftw_execute(plan, Q, Q_hat)
        call dfftw_destroy_plan(plan)
        ! Q_hat(:, :, :) = Q_hat(:, :, :)/(NX*NY*NZ)
        
    end subroutine fft_forward

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
        real(8) Re(1:NX, 1:NY, 1:NZ)
        complex(8) Im(1:NX/2+1, 1:NY, 1:NZ)
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
        ! call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, U_tmp(1:NX, 1:NY, 1:NZ), U_hat, FFTW_ESTIMATE)  ! 分けて書かないとバグった
        ! call dfftw_execute(plan, U_tmp(1:NX, 1:NY, 1:NZ), U_hat)
        ! call dfftw_destroy_plan(plan)

        ! call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, V_tmp(1:NX, 1:NY, 1:NZ), V_hat, FFTW_ESTIMATE)
        ! call dfftw_execute(plan, V_tmp(1:NX, 1:NY, 1:NZ), V_hat)
        ! call dfftw_destroy_plan(plan)

        ! call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, W_tmp(1:NX, 1:NY, 1:NZ), W_hat, FFTW_ESTIMATE)
        ! call dfftw_execute(plan, W_tmp(1:NX, 1:NY, 1:NZ), W_hat)
        ! call dfftw_destroy_plan(plan)

        ! call fft_forward(U_tmp(1:NX, 1:NY, 1:NZ), U_hat)
        ! call fft_forward(V_tmp(1:NX, 1:NY, 1:NZ), V_hat)
        ! call fft_forward(W_tmp(1:NX, 1:NY, 1:NZ), W_hat)

        call fftr2c_3d(U_tmp(1:NX, 1:NY, 1:NZ), U_hat)
        call fftr2c_3d(V_tmp(1:NX, 1:NY, 1:NZ), V_hat)
        call fftr2c_3d(W_tmp(1:NX, 1:NY, 1:NZ), W_hat)

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
        ! open(10, file='fft/'//str//'.d')
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

    subroutine energy_sum(U, V, W, U_ave, V_ave, W_ave, Energy)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: U_ave(0:NX+1, 0:NY+1, 0:NZ+1), V_ave(0:NX+1, 0:NY+1, 0:NZ+1), W_ave(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: Energy(0:NX)
        complex(8) U_hat(1:NX/2+1, 1:NY, 1:NZ), V_hat(1:NX/2+1, 1:NY, 1:NZ), W_hat(1:NX/2+1, 1:NY, 1:NZ)
        real(8) U_tmp(0:NX+1, 0:NY+1, 0:NZ+1), V_tmp(0:NX+1, 0:NY+1, 0:NZ+1), W_tmp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) E_tmp(1:NX/2+1, 1:NY, 1:NZ)
        real(8) K_abs(1:NX/2+1, 1:NY, 1:NZ)
        integer i, j, k, index
        integer(8) plan

        ! 平均を引き変動成分を求める
        U_tmp(:, :, :) = U_tmp(:, :, :) - U_ave(:, :, :)
        V_tmp(:, :, :) = V_tmp(:, :, :) - V_ave(:, :, :)
        W_tmp(:, :, :) = W_tmp(:, :, :) - W_ave(:, :, :)

        ! 速度場をフーリエ変換
        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, U(1:NX, 1:NY, 1:NZ), U_hat, FFTW_ESTIMATE)
        call dfftw_execute(plan, U(1:NX, 1:NY, 1:NZ), U_hat)
        call dfftw_destroy_plan(plan)

        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, V(1:NX, 1:NY, 1:NZ), V_hat, FFTW_ESTIMATE)
        call dfftw_execute(plan, V(1:NX, 1:NY, 1:NZ), V_hat)
        call dfftw_destroy_plan(plan)

        call dfftw_plan_dft_r2c_3d(plan, NX, NY, NZ, W(1:NX, 1:NY, 1:NZ), W_hat, FFTW_ESTIMATE)
        call dfftw_execute(plan, W(1:NX, 1:NY, 1:NZ), W_hat)
        call dfftw_destroy_plan(plan)
        
        U_hat(:, :, :) = U_hat(:, :, :)/dble(NX*NY*NZ/2)  ! 格子点数を2で割ると速度の2乗の平均に一致する
        V_hat(:, :, :) = V_hat(:, :, :)/dble(NX*NY*NZ/2)
        W_hat(:, :, :) = W_hat(:, :, :)/dble(NX*NY*NZ/2)


        ! 変動成分の二乗を足す
        E_tmp(:, :, :) = 0.0d0
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX/2+1
                    E_tmp(i, j, k) = E_tmp(i, j, k) + real(U_hat(i, j, k))**2 + imag(U_hat(i, j, k))**2
                    E_tmp(i, j, k) = E_tmp(i, j, k) + real(V_hat(i, j, k))**2 + imag(V_hat(i, j, k))**2
                    E_tmp(i, j, k) = E_tmp(i, j, k) + real(W_hat(i, j, k))**2 + imag(W_hat(i, j, k))**2
                enddo
            enddo
        enddo

        ! 配列要素に対する波数
        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX/2+1
                    if (i==1 .and. j==1 .and. k==1) then
                        K_abs(i, j, k) = 0.0d0
                    else
                        K_abs(i, j, k) = sqrt((i-1)**2 + (j-1)**2 + (k-1)**2.0d0)
                    endif
                enddo
            enddo
        enddo

        ! 波数を四捨五入し、対応する整数の波数にエネルギーを足し合わせる
        do k = 1, NZ/2+1  ! 波数は半分しか計算できない
            do j = 1, NY/2+1
                do i = 1, NX/2+1
                    index = nint(K_abs(i, j, k))
                    Energy(index) = Energy(index) + E_tmp(i, j, k)/2
                enddo
            enddo
        enddo


    end subroutine energy_sum

end module fft


module ibm
    use smac
    use fft
    implicit none
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
            X(i, :, :) = (i-0.5d0)*dX_C
        enddo
        do j = 0, NY+1
            Y(:, j, :) = (j-0.5d0)*dY_C
        enddo
        do k = 0, NZ+1
            Z(:, :, k) = (k-0.5d0)*dZ_C
        enddo

        do j = 1, NL
            Xc(1, j, :) = Xmax/4 + D_C/2 * cos(2*PI/NL*j)
            Yc(1, j, :) = Ymax/2 + D_C/2 * sin(2*Pi/NL*j)
        enddo

        ! 円柱上の座標
        ! if (NC == 4) then
        !     do j = 1, NL
        !         Xc(1, j, :) = Xmax/4 + D_C/2 * cos(2*PI/NL*j)
        !         Yc(1, j, :) = Ymax/4 + D_C/2 * sin(2*Pi/NL*j)
        !         Xc(2, j, :) = Xmax*3/4 + D_C/2 * cos(2*PI/NL*j)
        !         Yc(2, j, :) = Ymax/4 + D_C/2 * sin(2*Pi/NL*j)
        !         Xc(3, j, :) = Xmax/4 + D_C/2 * cos(2*PI/NL*j)
        !         Yc(3, j, :) = Ymax*3/4 + D_C/2 * sin(2*Pi/NL*j)
        !         Xc(4, j, :) = Xmax*3/4 + D_C/2 * cos(2*PI/NL*j)
        !         Yc(4, j, :) = Ymax*3/4 + D_C/2 * sin(2*Pi/NL*j)
        !     enddo
        ! else if (NC == 1) then
        !     do j = 1, NL
        !         Xc(1, j, :) = Xmax/4 + D_C/2 * cos(2*PI/NL*j)
        !         Yc(1, j, :) = Ymax/2 + D_C/2 * sin(2*Pi/NL*j)
        !     enddo
        ! endif

        ! do k = 1, NZ
        !     Zc(:, :, k) = (k-0.5d0)*dZ_C
        ! enddo

        ! ! 座標を0から1に規格化
        ! X(:, :, :) = X(:, :, :) / Xmax
        ! Y(:, :, :) = Y(:, :, :) / Ymax
        ! Z(:, :, :) = Z(:, :, :) / Zmax
        ! Xc(:, :, :) = Xc(:, :, :) / Xmax
        ! Yc(:, :, :) = Yc(:, :, :) / Ymax
        ! Zc(:, :, :) = Zc(:, :, :) / Zmax

        ! ! 円柱上の座標での速度
        ! Uc(:, :, :) = 0.0d0
        ! Vc(:, :, :) = 0.0d0
        ! Wc(:, :, :) = 0.0d0
        ! if (NC == 4) then
        !     do j = 1, NL
        !         Uc(1, j, :) = -U_C * sin(2*Pi/NL*j)
        !         Vc(1, j, :) = U_C * cos(2*Pi/NL*j)
        !         Uc(2, j, :) = U_C * sin(2*Pi/NL*j)
        !         Vc(2, j, :) = -U_C * cos(2*Pi/NL*j)
        !         Uc(3, j, :) = U_C * sin(2*Pi/NL*j)
        !         Vc(3, j, :) = -U_C * cos(2*Pi/NL*j)
        !         Uc(4, j, :) = -U_C * sin(2*Pi/NL*j)
        !         Vc(4, j, :) = U_C * cos(2*Pi/NL*j)
        !     enddo
        ! ! else if (NC == 1) then
        ! !     do j = 1, NL
        ! !         Uc(1, j, :) = -U_C * sin(2*Pi/NL*j)
        ! !         Vc(1, j, :) = U_C * cos(2*Pi/NL*j)
        ! !     enddo
        ! endif
        ! ! 速度を0から1に規格化
        ! Uc(:, :, :) = Uc(:, :, :) / U_C
        ! Vc(:, :, :) = Vc(:, :, :) / U_C
        ! Wc(:, :, :) = Wc(:, :, :) / U_C
    end subroutine ibm_init

    subroutine ibm_get(Xc, Yc)
        real(8), intent(in) :: Xc(1:NC, 1:NL, 1:NZ), Yc(1:NC, 1:NL, 1:NZ)
        integer i, j, k
        character(*), parameter :: str1='ibm_dat/'
        character(*), parameter :: str2='.d'
        character(8) str
        character(64) filename
        
        k = NZ/2
        do i = 1, NC
            write(str, '(I8.8)') i
            filename = str1//str//str2
            open(10, file=filename)
            do j = 1, NL
                write(10, '(2e12.4)') Xc(i, j, k)*Xmax, Yc(i, j, k)*Ymax
            enddo
            write(10, '(2e12.4)') Xc(i, 1, k)*Xmax, Yc(i, 1, k)*Ymax  ! 1周させる
            close(10)
        enddo
    end subroutine ibm_get

    subroutine ibm_preliminary(U, V, W, P, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Ua, Va, Wa, Fxc, Fyc, Fzc, Up, Vp, Wp, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
        real(8), intent(in) :: Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)
        real(8), intent(out) :: Ua(0:NX+1, 0:NY+1, 0:NZ+1), Va(0:NX+1, 0:NY+1, 0:NZ+1), Wa(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8), intent(inout) :: Fxc(1:NC, 1:NL, 1:NZ), Fyc(1:NC, 1:NL, 1:NZ), Fzc(1:NC, 1:NL, 1:NZ)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1, 0:NZ+1), Vp(0:NX+1, 0:NY+1, 0:NZ+1), Wp(0:NX+1, 0:NY+1, 0:NZ+1)
        integer, intent(in) :: step
        integer i, j, k
        
        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :, :) = Ax(:, :, :)
            Ay0(:, :, :) = Ay(:, :, :)
            Az0(:, :, :) = Az(:, :, :)
        endif

        do k = 1, NZ
            do j = 1, NY
                do i = 1, NX
                    Ua(i, j, k) = U(i, j, k) - dt*(-P(i, j, k)+P(i+1, j, k))/dX &
                                 +dt*(3.0d0*Ax(i, j, k) - Ax0(i, j, k))/2.0d0 &
                                 +dt*Bx(i, j, k)

                    Va(i, j, k) = V(i, j, k) - dt*(-P(i, j, k)+P(i, j+1, k))/dY &
                                 +dt*(3.0d0*Ay(i, j, k) - Ay0(i, j, k))/2.0d0 &
                                 +dt*By(i, j, k)

                    Wa(i, j, k) = W(i, j, k) - dt*(-P(i, j, k)+P(i, j, k+1))/dZ &
                                 +dt*(3.0d0*Az(i, j, k) - Az0(i, j, k))/2.0d0 &
                                 +dt*Bz(i, j, k)
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

    subroutine ibm_Helmholtz_kabe(Up, Vp, Wp, X, Y, Z, Xc, Yc, Zc, Uc, Vc, Wc, Fxc, Fyc, Fzc, fxint, fyint, fzint)
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
        
        Ub(:, :, :) = 0.0d0
        Vb(:, :, :) = 0.0d0
        Wb(:, :, :) = 0.0d0

        ! do n = 1, NZ
        !     do m = 1, NL
        !         do l = 1, NC
        !             ! 外力点周囲の3*3*3点の和をとる
        !             do k = int(Zc(l, m, n)/dZ), int(Zc(l, m, n)/dZ) + 2
        !                 do j = int(Yc(l, m, n)/dY), int(Yc(l, m, n)/dY) + 2
        !                     do i = int(Xc(l, m, n)/dX - 0.5d0), int(Xc(l, m, n)/dX - 0.5d0) + 2
        !                         Ub(l, m, n) = Ub(l, m, n) &
        !                                     + Up(i, j, k) * delta(X(i, j, k)+0.5d0*dX, Y(i, j, k), Z(i, j, k), &
        !                                                           Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
        !                     enddo
        !                 enddo
        !             enddo
        !             do k = int(Zc(l, m, n)/dZ), int(Zc(l, m, n)/dZ) + 2
        !                 do j = int(Yc(l, m, n)/dY - 0.5d0), int(Yc(l, m, n)/dY - 0.5d0) + 2
        !                     do i = int(Xc(l, m, n)/dX), int(Xc(l, m, n)/dX) + 2
        !                         Vb(l, m, n) = Vb(l, m, n) &
        !                                     + Vp(i, j, k) * delta(X(i, j, k), Y(i, j, k)+0.5d0*dY, Z(i, j, k), &
        !                                                           Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
        !                     enddo
        !                 enddo
        !             enddo
        !             do k = int(Zc(l, m, n)/dZ - 0.5d0), int(Zc(l, m, n)/dZ - 0.5d0) + 2
        !                 do j = int(Yc(l, m, n)/dY) , int(Yc(l, m, n)/dY) + 2
        !                     do i = int(Xc(l, m, n)/dX), int(Xc(l, m, n)/dX) + 2
        !                         Wb(l, m, n) = Wb(l, m, n) &
        !                                     + Wp(i, j, k) * delta(X(i, j, k), Y(i, j, k), Z(i, j, k)+0.5d0*dZ, &
        !                                                           Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
        !                     enddo
        !                 enddo
        !             enddo
        !         enddo
        !     enddo
        ! enddo

        do n = 1, NZ
            do m = 1, NL
                do l = 1, NC
                    ! 外力点周囲の3*3*3点の和をとる
                    do k = n-1, n+1
                        do j = 1, NY
                            do i = 1, NX
                                Ub(l, m, n) = Ub(l, m, n) &
                                            + Up(i, j, k) * delta(X(i, j, k)+0.5d0*dX, Y(i, j, k), Z(i, j, k), &
                                                                  Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
                                Vb(l, m, n) = Vb(l, m, n) &
                                            + Vp(i, j, k) * delta(X(i, j, k), Y(i, j, k)+0.5d0*dY, Z(i, j, k), &
                                                                  Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
                                Wb(l, m, n) = Wb(l, m, n) &
                                            + Wp(i, j, k) * delta(X(i, j, k), Y(i, j, k), Z(i, j, k)+0.5d0*dZ, &
                                                                  Xc(l, m, n), Yc(l, m, n), Zc(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! write(*, *) sum(Ub)/(NZ*NL*NC)
        ! write(*, *) sum(Vb)/(NZ*NL*NC)
        ! write(*, *) sum(Wb)/(NZ*NL*NC)

        do n = 1, NZ
            do m = 1, NL
                do l = 1, NC
                    Fxc(l, m, n) = Fxc(l, m, n) + (Uc(l, m, n) - Ub(l, m, n))/dt
                    Fyc(l, m, n) = Fyc(l, m, n) + (Vc(l, m, n) - Vb(l, m, n))/dt
                    Fzc(l, m, n) = Fzc(l, m, n) + (Wc(l, m, n) - Wb(l, m, n))/dt
                enddo
            enddo
        enddo

        ! ! デバッグ用
        ! er = 0.0d0
        ! do n = 1, NZ
        !     do m = 1, NL
        !         do l = 1, NC
        !             er = er + (Uc(l, m, n) - Ub(l, m, n))**2
        !             er = er + (Vc(l, m, n) - Vb(l, m, n))**2
        !             er = er + (Wc(l, m, n) - Wb(l, m, n))**2
        !         enddo
        !     enddo
        ! enddo
        ! write(20, '(1e12.4)') sqrt(er/(3*NC*NL*NC))



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


    end subroutine ibm_Helmholtz_kabe

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
                    LHS(i, j, k) = 1 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX*2*PI))/dX**2 &
                                                  + 2*(1-cos((j-1)*dY*2*PI))/dY**2 &
                                                  + 2*(1-cos((k-1)*dZ*2*PI))/dZ**2)
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
    real(8) Ax(1:NX, 1:NY, 1:NZ), Ay(1:NX, 1:NY, 1:NZ), Az(1:NX, 1:NY, 1:NZ)  ! 対流項の計算
    real(8) Ax0(1:NX, 1:NY, 1:NZ), Ay0(1:NX, 1:NY, 1:NZ), Az0(1:NX, 1:NY, 1:NZ)
    real(8) Bx(1:NX, 1:NY, 1:NZ), By(1:NX, 1:NY, 1:NZ), Bz(1:NX, 1:NY, 1:NZ)  ! 粘性項の計算
    real(8) Bx0(1:NX, 1:NY, 1:NZ), By0(1:NX, 1:NY, 1:NZ), Bz0(1:NX, 1:NY, 1:NZ)
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
    ! エネルギーカスケードのために保存
    real(8) U_tmp(0:NX+1, 0:NY+1, 0:NZ+1), V_tmp(0:NX+1, 0:NY+1, 0:NZ+1), W_tmp(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) P_tmp(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) U_ave(0:NX+1, 0:NY+1, 0:NZ+1), V_ave(0:NX+1, 0:NY+1, 0:NZ+1), W_ave(0:NX+1, 0:NY+1, 0:NZ+1)
    real(8) Energy(0:NX)
    integer step, i, j, k, s
    real(8) t_start, t_end, total_time, t1, t2, t12
    call mk_dir(dir)
    call cpu_time(t_start)
    t12 = 0.0d0

    call fft_init
    call init(U, V, W, P, Phi, Fx, Fy, Fz)  ! 初期条件
    if (input_step > 0) call input(U, V, W, P)
    if (method == 2) call ibm_init(X, Y, Z, Xc, Yc, Zc, Uc, Vc, Wc)
    if (method == 2) call ibm_get(Xc, Yc)
    ! call energy_single(U, V, W, 0)

    do step = 1, Nstep
        ! if (mod(step, Gstep)==0) write(*, *) 'step =', step
        call convection(U, V, W, Ax, Ay, Az)
        call viscous(U, V, W, Bx, By, Bz)

        if (method == 0) then
            call navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Bx, By, Bz, Ax0, Ay0, Az0, Bx0, By0, Bz0, Fx, Fy, Fz, step)
            call poisson(Up, Vp, Wp, Phi)
            call march(Up, Vp, Wp, U, V, W, Phi, P)
        endif

        if (method == 1 .or. method == 3) then
            ! call fft_navier(U, V, W, P, U_tmp, V_tmp, W_tmp, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Fx, Fy, Fz, step)
            call fft_navier_kabe(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Fx, Fy, Fz, step)
            ! write(*, *) 'Up_er =', sqrt(sum((U_tmp(1:NX, 1:NY, 1:NZ)-Up(1:NX, 1:NY, 1:NZ))**2)/(NX*NY*NZ))
            ! write(*, *) 'Vp_er =', sqrt(sum((V_tmp(1:NX, 1:NY, 1:NZ)-Vp(1:NX, 1:NY, 1:NZ))**2)/(NX*NY*NZ))
            ! write(*, *) 'Wp_er =', sqrt(sum((W_tmp(1:NX, 1:NY, 1:NZ)-Wp(1:NX, 1:NY, 1:NZ))**2)/(NX*NY*NZ))

            ! call fft_poisson(Up, Vp, Wp, P)
            call fft_poisson_kabe(Up, Vp, Wp, Phi)

            ! write(*, *) 'Phier =', sqrt(sum((P(1:NX, 1:NY, 1:NZ)-Phi(1:NX, 1:NY, 1:NZ))**2)/(NX*NY*NZ))
            ! do k = 1, NZ
            !     do j = 1, NY
            !         do i = 1, NX
            !             write(*, '(3I4, 3E12.4)') i, j, k, P(i, j, k), Phi(i, j, k), P(i, j, k)-Phi(i, j, k)
            !         enddo
            !     enddo
            ! enddo

            call fft_march(Up, Vp, Wp, U, V, W, Phi, P)
            if (mod(step, 50)==0) then
                write(*, '(I7, 10F8.3)') step, sum(Up(1:NX, 1:NY, 0))/(NX*NY)*U_C, sum(Up(1:NX, 1:NY, 1))/(NX*NY)*U_C, &
                                               sum(Up(1:NX, 1:NY, NZ/2))/(NX*NY)*U_C, sum(Up(1:NX, 1:NY, NZ/2+1))/(NX*NY)*U_C, &
                                               sum(Up(1:NX, 1:NY, NZ))/(NX*NY)*U_C, sum(Up(1:NX, 1:NY, NZ+1))/(NX*NY)*U_C
            endif
        endif

        if (method == 2) then
            call ibm_preliminary(U, V, W, P, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Ua, Va, Wa, Fxc, Fyc, Fzc, Up, Vp, Wp, step)
            call ibm_Helmholtz(Up, Vp, Wp, X, Y, Z, Xc, Yc, Zc, Uc, Vc, Wc, Fxc, Fyc, Fzc, fxint, fyint, fzint)
            call ibm_predict(Ua, Va, Wa, fxint, fyint, fzint, Bx, By, Bz, Up, Vp, Wp)
            call fft_poisson(Up, Vp, Wp, Phi)
            call fft_march(Up, Vp, Wp, U, V, W, Phi, P)
        endif

        if (mod(step, Gstep)==0) call get_data(U, V, W, step)
        ! if (mod(step, Gstep)==0) call energy_single(U, V, W, step)
        if (mod(step, output_step)==0) call output(U, V, W, P, step)
        t12 = t12 + t2-t1

        ! if (mod(step, 100)==0) then
        !     write(*, '(I5, 6e12.4)') step, maxval(abs(U(:, :, :)))*U_C, maxval(abs(V(:, :, :)))*U_C, maxval(abs(W(:, :, :)))*U_C
        ! endif
    enddo

    ! call get_data(U, V, W, 0)

    if (method == 3) then  ! 定常状態に落ち着いていると仮定し、エネルギーカスケードを求める
        ! 平均値を得るため初期条件を保存して数ステップ回す
        U_tmp(:, :, :) = U(:, :, :)
        V_tmp(:, :, :) = V(:, :, :)
        W_tmp(:, :, :) = W(:, :, :)
        P_tmp(:, :, :) = P(:, :, :)
        U_ave(:, :, :) = 0.0d0
        V_ave(:, :, :) = 0.0d0
        W_ave(:, :, :) = 0.0d0
        ! do step = 1, NE
        !     call convection(U, V, W, Ax, Ay, Az)
        !     call viscous(U, V, W, Bx, By, Bz)
        !     call fft_navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Fx, Fy, Fz, step)
        !     call fft_poisson(Up, Vp, Wp, Phi)
        !     call fft_march(Up, Vp, Wp, U, V, W, Phi, P)
        !     ! 平均を計算
        !     U_ave(:, :, :) = U_ave(:, :, :) + U(:, :, :)
        !     V_ave(:, :, :) = V_ave(:, :, :) + V(:, :, :)
        !     W_ave(:, :, :) = W_ave(:, :, :) + W(:, :, :)
        ! enddo
        ! U_ave(:, :, :) = U_ave(:, :, :)/NE
        ! V_ave(:, :, :) = V_ave(:, :, :)/NE
        ! W_ave(:, :, :) = W_ave(:, :, :)/NE


        ! 平均を引いた変動成分からエネルギーを求める
        U(:, :, :) = U_tmp(:, :, :)
        V(:, :, :) = V_tmp(:, :, :)
        W(:, :, :) = W_tmp(:, :, :)
        P(:, :, :) = P_tmp(:, :, :)
        Energy(:) = 0.0d0
        do step = 1, NE
            call convection(U, V, W, Ax, Ay, Az)
            call viscous(U, V, W, Bx, By, Bz)
            call fft_navier(U, V, W, P, Up, Vp, Wp, Ax, Ay, Az, Ax0, Ay0, Az0, Bx, By, Bz, Fx, Fy, Fz, step)
            call fft_poisson(Up, Vp, Wp, Phi)
            call fft_march(Up, Vp, Wp, U, V, W, Phi, P)
            call energy_sum(U, V, W, U_ave, V_ave, W_ave, Energy)  ! エネルギーの足し算
        enddo
        Energy(:) = Energy(:)/NE

        open(10, file='fft/energy_v1.d')
        do i = 0, NX
            write(10, '(I4, e12.4)') i, Energy(i)
        enddo
        close(10)

    endif

    call fft_finalize

    call cpu_time(t_end)
    total_time = t_end - t_start
    write(*, '(a9, F10.3, a3)') 'Total   :', total_time, '[s]'
    if (method == 2) write(*, '(a9, F10.3, a3, F8.3, a3)') 'Helmholtz:', t12, '[s]', t12/total_time*100, '[%]'
    
end program main