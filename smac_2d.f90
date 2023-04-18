module smac
    implicit none
    ! gfortran smac_2d.f90 -I$HOME/local/include -L$HOME/local/lib -lfftw3 && ./a.out
    ! ステップ数
    integer, parameter :: Nstep = 10
    integer, parameter :: Gstep = 1  ! データを取得する間隔
    integer, parameter :: Estep = 1  ! エネルギースペクトルを取得する間隔
    integer, parameter :: Dstep = 1  ! デバッグする間隔
    character(*), parameter :: dir = './2d/'
    integer, parameter :: input_step = 0  ! 0以外で初期条件をファイルから読み込む
    integer, parameter :: output_step = 1000  ! 配列を保存する間隔
    ! 手法
    integer, parameter :: method = 1  ! 0:陽解法、1:FFT、2:IBM
    real(8), parameter :: PI = acos(-1.0d0)
    ! パラメータ
    integer, parameter :: NX = 128, NY = NX
    real(8), parameter :: dX = 2*PI/NX, dY = 2*PI/NY  ! 規格化長さは2*PI
    real(8), parameter :: dt = 0.01d0
    ! method = 1
    real(8), parameter :: Re = 1.0d0
    real(8), parameter :: L_C = 1.0d0  ! 長さが100なら100/2*PI
    real(8), parameter :: U_C = 1.0d0  ! 本来は乱流テイラーグリーン渦の平均流の速さ  ! ibmでの円柱回転速度
    real(8), parameter :: nu = L_C*U_C/Re
    ! method = 2
    ! real(8), parameter :: L_C = 8.0d0 * 0.03d0 / (2.0d0*PI)
    ! real(8), parameter :: D_C = 0.03d0
    ! real(8), parameter :: U_C = 0.06d0 * PI
    ! real(8), parameter :: nu = 1.0d-6
    ! real(8), parameter :: Re = U_C*L_C/nu

    real(8), parameter :: f0 = 1.0d0  ! 無次元化したときに1となるように
    real(8), parameter :: dX_C = dX*L_C, dY_C = dY*L_C
    real(8), parameter :: dt_C = dt*L_C/U_C
    ! real(8), parameter :: f0 = L_C/(U_C)**2 ! 外力を1に固定し，無次元化するときの定数

    ! エネルギーカスケード
    ! integer, parameter :: NE = 1000  ! 時間平均反復回数

contains
    subroutine init(U, V, P, Phi, Fx, Fy)  ! はじめに実行
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(out) :: P(0:NX+1, 0:NY+1), Phi(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        integer i, j
        real(8), parameter :: large_K = 2.0d0

        if (method == 0) then
            write(*, '(a, F8.3, F8.3, a8)') 'Von Neumann:', 1.0d0/Re*dt/dX**2, 1.0d0/Re*dt/dY**2, '< 0.167'
        endif
        write(*, '(a, F8.3, F8.3, a8)') 'CFL:', 1.0d0*dt/dX, 1.0d0*dt/dY,'< 0.167'
        ! write(*, '(a, F8.3)') 'L_C =', L_C
        ! write(*, '(a, F8.3)') 'D_C =', D_C
        ! write(*, '(a, F8.3)') 'U_C =', U_C
        ! write(*, '(a, E12.4)') 'nu =', nu
        ! write(*, '(a, F8.3)') 'Re  =', Re
        ! write(*, '(a, E12.4)') 'dt_C =', dt_C

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
        do j = 1, NY
            do i = 1, NX
                U(i, j) = -sin(i*dX) * cos((j-0.5d0)*dY)
                V(i, j) = cos((i-0.5d0)*dX) * sin(j*dY)
                ! Fx(i, j) = -sin(i*dX) * cos((j-0.5d0)*dY)  ! 増田さん, 安房井さん
                ! Fy(i, j) = cos((i-0.5d0)*dX) * sin(j*dY)
                ! Fx(i, j) = -sin(large_K*(j-0.5d0)*dY)  ! 小井手さん
                ! Fy(i, j) = sin(large_K*(i-0.5d0)*dX)
            enddo
        enddo
        Fx(:, :) = f0 * Fx(:, :)
        Fy(:, :) = f0 * Fy(:, :)
        call PBM(U)
        call PBM(V)
        call PBM(P)
     end subroutine init

    subroutine PBM(A)
        real(8), intent(inout) :: A(0:NX+1, 0:NY+1)
        A(0, :) = A(NX, :)
        A(NX+1, :) = A(1, :)
        A(:, 0) = A(:, NY)
        A(:, NY+1) = A(:, 1)
    end subroutine PBM

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

        Bx(:, :) = 1.0d0/Re*Bx(:, :)
        By(:, :) = 1.0d0/Re*By(:, :)
    end subroutine viscous

    subroutine navier(U, V, P, Up, Vp, Ax, Ay, Bx, By, Ax0, Ay0, Bx0, By0, Fx, Fy, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
        real(8), intent(in) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
        real(8), intent(inout) :: Bx0(1:NX, 1:NY), By0(1:NX, 1:NY)
        real(8), intent(in) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        integer, intent(in) :: step
        integer i, j

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :) = Ax(:, :)
            Ay0(:, :) = Ay(:, :)
            Bx0(:, :) = Bx(:, :)
            By0(:, :) = By(:, :)
        endif

        ! NS方程式で速度の予測
        do j = 1, NY
            do i = 1, NX
                Up(i, j) = U(i, j) - dt*(-P(i, j)+P(i+1, j))/dX &
                                    +dt*(3.0d0*(Ax(i, j) + Bx(i, j)) - (Ax0(i, j) + Bx0(i, j)))/2.0d0 &
                                    +dt*Fx(i, j)
            enddo
        enddo

        do j = 1, NY
            do i = 1, NX
                Vp(i, j) = V(i, j) - dt*(-P(i, j)+P(i, j+1))/dY &
                                    +dt*(3.0d0*(Ay(i, j) + By(i, j)) - (Ay0(i, j) + By0(i, j)))/2.0d0 &
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

    subroutine taylor_debag(U, V, step)
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
        open(10, file=dir//'debag2.d', position='append')
        write(10, *) step, U0(NX/4, NY/4), U(NX/4, NY/4), U0(NX/4, NY/4) - U(NX/4, NY/4), sqrt(err)
        close(10)
    end subroutine taylor_debag


    subroutine get_data(U, V, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        integer i, j
        real(8) D(2, 2), Omega, S
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        open(10, file=dir//'z_'//str//'.d')

        do j = 1, NY
            do i = 1, NX
                D(1, 1) = (U(i, j) - U(i-1, j))/dX
                D(1, 2) = (U(i, j+1) - U(i, j-1) + U(i-1, j+1) - U(i-1, j-1))/(4*dY)
                D(2, 1) = (V(i+1, j) - V(i-1, j) + V(i+1, j-1) - V(i-1, j-1))/(4*dX)
                D(2, 2) = (V(i, j) - V(i, j-1))/dY
                D(:, :) = D(:, :)*U_C
                Omega = (D(2, 1) - D(1, 2))
                S = (D(2, 1) + D(1, 2))
                write(10, '(100e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, 0.0d0, &
                                    (U(i-1, j)+U(i, j))/2*U_C, (V(i, j-1)+V(i, j))/2*U_C, 0.0d0, &
                                    Omega, 0.0d0, S, 0.0d0, 0.0d0
            enddo
        enddo
        close(10)
    end subroutine get_data


    subroutine logging(U, V, step, t_start)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        real(8), intent(in) :: t_start
        real(8) t_temp
        real(8) K_energy
        integer(8) t_pre, hour, min, sec

        if (mod(step, Dstep) == 0) then
            write(*, '(a, I6)', advance='no') 'step:', step

            ! 運動エネルギーとエネルギーカスケードの総和
            K_energy = sum(U(1:NX, 1:NY)**2) + sum(V(1:NX, 1:NY)**2)
            K_energy = K_energy*U_C**2/2
            K_energy = K_energy/(NX*NY)
            write(*, '(a, e12.4)', advance='no') '  | K_energy:', K_energy

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
        open(10, file=dir//str//'.d')

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
        write(*, '(a, a, e12.4)') 'input_file='//dir//str//'.d', '  | K_energy:', K_energy
    end subroutine input

    subroutine output(U, V, P, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        integer, intent(in) :: step
        integer i, j
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換
        open(10, file=dir//str//'.d')

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
    integer(8) plan1, plan2
    real(8) Re1(1:NX, 1:NY)
    complex(8) Im1(1:NX/2+1, 1:NY)
contains
    subroutine fft_init
        call dfftw_plan_dft_r2c_2d(plan1, NX, NY, Re1, Im1, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(plan2, NX, NY, Im1, Re1, FFTW_ESTIMATE)
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

    subroutine fft_finalize
        call dfftw_destroy_plan(plan1)
        call dfftw_destroy_plan(plan2)
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

    subroutine fft_navier(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, Bx, By, Fx, Fy, step)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(in) :: P(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
        real(8), intent(inout) :: Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
        real(8), intent(in) :: Bx(1:NX, 1:NY), By(1:NX, 1:NY)
        real(8), intent(in) :: Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
        integer, intent(in) :: step
        integer i, j
        real(8) Q(1:NX, 1:NY), LHS(1:NX, 1:NY)

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :) = Ax(:, :)
            Ay0(:, :) = Ay(:, :)
        endif

        ! 左辺の定数部分の計算
        do j = 1, NY
            do i = 1, NX
                LHS(i, j) = 1 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                            +2*(1-cos((j-1)*dY))/dY**2 )
            enddo
        enddo

        ! 速度場Upを予測
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = U(i, j) &
                            -dt*(-P(i, j)+P(i+1, j))/dX &
                            +dt*(3.0d0*Ax(i, j) - Ax0(i, j))/2.0d0 &
                            +dt*Bx(i, j)/2.0d0 &
                            +dt*Fx(i, j)
            enddo
        enddo

        call fft_solve(Q, LHS, Up(1:NX, 1:NY))
        
        ! 同様にVpについても求める
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = V(i, j) &
                            -dt*(-P(i, j)+P(i, j+1))/dY &
                            +dt*(3.0d0*Ay(i, j) - Ay0(i, j))/2.0d0 &
                            +dt*By(i, j)/2.0d0 &
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
                        -dt/Re/2.0d0*((Phi(i-1, j)-2*Phi(i, j)+Phi(i+1, j))/dX**2 &
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
        ! open(10, file=dir//'fft_'//str//'.d')
        ! do i = 0, NX
        !     write(10, '(I4, e12.4)') i, Energy(i)
        ! enddo
        ! close(10)


        ! 運動エネルギーとエネルギーカスケードの総和
        K_energy = sum(U_tmp(1:NX, 1:NY)**2.0d0) + sum(V_tmp(1:NX, 1:NY)**2.0d0) 
        K_energy = K_energy/2.d0
        K_energy = K_energy/(NX*NY)
        write(*, '(I6, 5e12.4)') step, K_energy, sum(Energy(:)), K_energy-sum(Energy(:)), Energy(0), K_energy/sum(Energy(:))   ! 一致するはず
        ! open(30, file = 'fft/debag_energy.d', position='append')
        ! write(30, '(I6, 2e12.4)') step, K_energy, sum(Energy(:))
        ! close(30)

    end subroutine energy_single

end module fft


program main
    use smac
    use fft
    implicit none
    real(8) U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
    real(8) Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)  ! 対流項の計算
    real(8) Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
    real(8) Bx(1:NX, 1:NY), By(1:NX, 1:NY)  ! 粘性項の計算
    real(8) Bx0(1:NX, 1:NY), By0(1:NX, 1:NY)
    real(8) Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
    real(8) P(0:NX+1, 0:NY+1)
    real(8) Phi(0:NX+1, 0:NY+1)
    real(8) Fx(1:NX, 1:NY), Fy(1:NX, 1:NY)
    integer step
    real(8) t_start, t_end, total_time, t1, t2, t12
    call mk_dir(dir)
    call cpu_time(t_start)
    t12 = 0.0d0

    call fft_init
    call init(U, V, P, Phi, Fx, Fy)  ! 初期条件
    if (input_step > 0) call input(U, V, P)
    ! call energy_single(U, V, W, 0)

    do step = 1, Nstep
        call convection(U, V, Ax, Ay)
        call viscous(U, V, Bx, By)

        if (method == 0) then
            call navier(U, V, P, Up, Vp, Ax, Ay, Bx, By, Ax0, Ay0, Bx0, By0, Fx, Fy, step)
            call poisson(Up, Vp, Phi)
            call march(Up, Vp, U, V, Phi, P)
        endif

        if (method == 1) then
            call fft_navier(U, V, P, Up, Vp, Ax, Ay, Ax0, Ay0, Bx, By, Fx, Fy, step)
            call fft_poisson(Up, Vp, Phi)
            call fft_march(Up, Vp, U, V, Phi, P)
        endif

        ! if (mod(step, Gstep)==0) call get_data(U, V, step)
        call taylor_debag(U, V, step)
        call logging(U, V, step, t_start)
        if (mod(step, Estep)==0) call energy_single(U, V, step)
        if (mod(step, output_step)==0) call output(U, V, P, step)
        ! if (sum(U**2)*0.0d0 /= 0.0d0) stop 'NaN value'  ! NaNの判定
        t12 = t12 + t2-t1

        ! if (mod(step, 100)==0) then
        !     write(*, '(I5, 6e12.4)') step, maxval(abs(U(:, :)))*U_C, maxval(abs(V(:, :)))*U_C, maxval(abs(W(:, :)))*U_C
        ! endif
    enddo


    call fft_finalize
    call cpu_time(t_end)
    total_time = t_end - t_start
    write(*, '(a9, F10.3, a3)') 'Total   :', total_time, '[s]'
    if (method == 2) write(*, '(a9, F10.3, a3, F8.3, a3)') 'Helmholtz:', t12, '[s]', t12/total_time*100, '[%]'
    
end program main

