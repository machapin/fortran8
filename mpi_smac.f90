module smac
    implicit none
    ! module load fftw
    ! mpifrtpx mpi_smac.f90 -lfftw3 -Kfast
    ! jxsub run.sh
    include 'mpif.h'
    ! ステップ数
    integer, parameter :: Nstep =   100
    integer, parameter :: Gstep =   40000  ! データを取得する間隔
    integer, parameter :: Estep =   40000  ! エネルギースペクトルを取得する間隔
    integer, parameter :: Dstep =    10  ! デバッグする間隔
    character(*), parameter :: dir = './'
    integer, parameter :: input_step = 0  ! 0以外で初期条件をファイルから読み込む
    integer, parameter :: output_step = 40000  ! 配列を保存する間隔
    ! 手法
    integer, parameter :: method = 1  ! 0:陽解法、1:FFT，2:IBM, 3:実験と比較(inputを0以外に)
    real(8), parameter :: PI = acos(-1.0d0)
    ! パラメータ
    integer, parameter :: NX = 128, NY = NX, NZ = NX
    real(8), parameter :: dX = 2*PI/NX, dY = 2*PI/NY, dZ = 2*PI/NZ  ! 規格化長さは2*PI
    real(8), parameter :: dt = 0.0004d0
    ! 無次元パラメータ
    real(8), parameter :: Re = 14000.0d0
    ! 有次元パラメータ
    real(8), parameter :: L_C = 1.0d0  ! 長さが100なら100/(2*PI)
    real(8), parameter :: U_C = 1.0d0  ! 本来は乱流テイラーグリーン渦の平均流の速さ
    ! real(8), parameter :: L_C = 8.0d0 * 0.03d0 / (2.0d0*PI)
    ! real(8), parameter :: U_C = Re/0.03d0*1.0d-6  ! L_Cは直径ではなく系の長さ．つまり代表長さではない
    real(8), parameter :: f0 = 1.0d0  ! 無次元化したときに1となるように
    real(8), parameter :: dX_C = dX*L_C, dY_C = dY*L_C, dZ_C = dZ*L_C
    real(8), parameter :: dt_C = dt*L_C/U_C
    real(8), parameter :: nu = L_C*U_C/Re

    ! MPI用変数
    integer N_procs
    integer ierr, procs, myrank
    integer next_rank, former_rank
    integer req1s, req1r, req2s, req2r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s, sta1r, sta2s, sta2r

contains
    subroutine init(U_procs, V_procs, W_procs, P_procs, Phi_procs, &
                    Ax_procs, Ay_procs, Az_procs, Bx_procs, By_procs, Bz_procs, &
                    Ax0_procs, Ay0_procs, Az0_procs, Bx0_procs, By0_procs, Bz0_procs, &
                    Up_procs, Vp_procs, Wp_procs, Fx_procs, Fy_procs, Fz_procs)  ! はじめに実行
        real(8), allocatable :: U_procs(:, :, :), V_procs(:, :, :), W_procs(:, :, :)
        real(8), allocatable :: P_procs(:, :, :), Phi_procs(:, :, :)
        real(8), allocatable :: Ax_procs(:, :, :), Ay_procs(:, :, :), Az_procs(:, :, :)
        real(8), allocatable :: Bx_procs(:, :, :), By_procs(:, :, :), Bz_procs(:, :, :)
        real(8), allocatable :: Ax0_procs(:, :, :), Ay0_procs(:, :, :), Az0_procs(:, :, :)
        real(8), allocatable :: Bx0_procs(:, :, :), By0_procs(:, :, :), Bz0_procs(:, :, :)
        real(8), allocatable :: Up_procs(:, :, :), Vp_procs(:, :, :), Wp_procs(:, :, :)
        real(8), allocatable :: Fx_procs(:, :, :), Fy_procs(:, :, :), Fz_procs(:, :, :)
        real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) P(0:NX+1, 0:NY+1, 0:NZ+1), Phi(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) Fx(1:NX, 1:NY, 1:NZ), Fy(1:NX, 1:NY, 1:NZ), Fz(1:NX, 1:NY, 1:NZ)
        integer i, j, k
        real(8), parameter :: large_K = 2.0d0

        ! コメント
        N_procs = NZ / procs
        if (mod(NZ, procs)/=0 .or. mod(NY, procs)/=0) stop 'NZ or NY is not divisible by procs'
        ! if (1.0d0*dt/dX>1.0d0/6 .or. 1.0d0*dt/dY>1.0d0/6 .or. 1.0d0*dt/dZ>1.0d0/6) stop 'CFL condition is not met.'
        if (myrank == 0) then
            write(*, '(a, F8.3, F8.3, F8.3)') 'CFL:', 1.0d0*dt/dX, 1.0d0*dt/dY, 1.0d0*dt/dZ
            ! write(*, '(a, F8.3)') 'U_C =', U_C
            ! write(*, '(a, F8.3)') 'Re  =', Re
            ! write(*, '(a, E12.4)') 'dt  =', dt
            ! write(*, '(a, E12.4)') 'nu  =', nu
            write(*, *) 'procs:', procs
        endif

        ! 初期条件
        if (myrank == 0) then
            U(:, :, :) = 0.0d0
            V(:, :, :) = 0.0d0
            W(:, :, :) = 0.0d0
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
                        ! U(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)
                        ! V(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)
                        ! U(i, j, k) = -sin(i*dX) * cos((k-0.5d0)*dZ)
                        ! W(i, j, k) = cos((i-0.5d0)*dX) * sin(k*dZ)

                        ! Fx(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY) * cos((k-0.5d0)*dZ)  ! 荒木さん
                        ! Fy(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY) * cos((k-0.5d0)*dZ)
                        Fx(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)  ! 増田さん, 安房井さん
                        Fy(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)
                        ! Fx(i, j, k) = -sin(large_K*(j-0.5d0)*dY)  ! 小井手さん
                        ! Fy(i, j, k) = sin(large_K*(i-0.5d0)*dX)
                    enddo
                enddo
            enddo
            Fx(:, :, :) = f0 * Fx(:, :, :)
            Fy(:, :, :) = f0 * Fy(:, :, :)
            Fz(:, :, :) = f0 * Fz(:, :, :)
        endif

        ! 配列のallocate
        N_procs = NZ / procs
        allocate(U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(P_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Phi_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Bx0_procs(1:NX, 1:NY, 1:N_procs), By0_procs(1:NX, 1:NY, 1:N_procs), Bz0_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Fx_procs(1:NX, 1:NY, 1:N_procs), Fy_procs(1:NX, 1:NY, 1:N_procs), Fz_procs(1:NX, 1:NY, 1:N_procs))
        
        ! 初期条件を全プロセスに通信
        call MPI_Scatter(U(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, U_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(V(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, V_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(W(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, W_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(P(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, P_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(Phi(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, Phi_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(Fx(1, 1, 1), NX*NY*N_procs, MPI_REAL8, Fx_procs(1, 1, 1), NX*NY*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(Fy(1, 1, 1), NX*NY*N_procs, MPI_REAL8, Fy_procs(1, 1, 1), NX*NY*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(Fz(1, 1, 1), NX*NY*N_procs, MPI_REAL8, Fz_procs(1, 1, 1), NX*NY*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        ! のりしろ境界の通信先
        next_rank = myrank + 1  ! 送るランク
        former_rank = myrank - 1  ! 受け取るランク
        if (myrank == procs - 1) then  ! 周期境界条件が一つでもあるならここは必要
            next_rank = 0
        else if (myrank == 0) then
            former_rank = procs - 1
        endif

        ! MPI通信
        call MPI_Boundary(U_procs)
        call MPI_Boundary(V_procs)
        call MPI_Boundary(W_procs)
        call MPI_Boundary(P_procs)
        call MPI_Boundary(Phi_procs)

        ! 境界条件
        call PBM(U_procs)
        call PBM(V_procs)
        call PBM(W_procs)
        call PBM(P_procs)
        call PBM(Phi_procs)
    end subroutine init


    subroutine MPI_Boundary(A_procs)  ! 通信したい行列を代入
        real(8), intent(inout) :: A_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        ! 次のランクへ送信
        call MPI_Isend(A_procs(0, 0, N_procs), (NX+2)*(NY+2), MPI_REAL8, next_rank, 1, MPI_COMM_WORLD, req1s, ierr)
        call MPI_Irecv(A_procs(0, 0, 0), (NX+2)*(NY+2), MPI_REAL8, former_rank, 1, MPI_COMM_WORLD, req1r, ierr)
  
        ! 前のランクへ送信
        call MPI_Isend(A_procs(0, 0, 1), (NX+2)*(NY+2), MPI_REAL8, former_rank, 2, MPI_COMM_WORLD, req2s, ierr)
        call MPI_Irecv(A_procs(0, 0, N_procs+1), (NX+2)*(NY+2), MPI_REAL8, next_rank, 2, MPI_COMM_WORLD, req2r, ierr)
  
        ! 待機
        call MPI_Wait(req1s, sta1s, ierr)
        call MPI_Wait(req1r, sta1r, ierr)
        call MPI_Wait(req2s, sta2s, ierr)
        call MPI_Wait(req2r, sta2r, ierr)
    end subroutine MPI_Boundary

    subroutine PBM(A_procs)  ! 一般的な周期境界条件
        real(8), intent(inout) :: A_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        A_procs(0, :, :) = A_procs(NX, :, :)  ! 周期境界条件
        A_procs(NX+1, :, :) = A_procs(1, :, :)
        A_procs(:, 0, :) = A_procs(:, NY, :)
        A_procs(:, NY+1, :) = A_procs(:, 1, :)
    end subroutine PBM


    subroutine convection(U_procs, V_procs, W_procs, Ax_procs, Ay_procs, Az_procs)   ! 発散型
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs)
        integer i, j, k
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Ax_procs(i, j, k) = (-((U_procs(i-1, j, k)+U_procs(i, j, k))/2)**2 + ((U_procs(i, j, k)+U_procs(i+1, j, k))/2)**2)/dX &
                                       +(-(V_procs(i, j-1, k)+V_procs(i+1, j-1, k))/2*(U_procs(i, j-1, k)+U_procs(i, j, k))/2 &
                                       +(V_procs(i, j, k)+V_procs(i+1, j, k))/2*(U_procs(i, j, k)+U_procs(i, j+1, k))/2)/dY &
                                       +(-(W_procs(i, j, k-1)+W_procs(i+1, j, k-1))/2*(U_procs(i, j, k-1)+U_procs(i, j, k))/2 &
                                       +(W_procs(i, j, k)+W_procs(i+1, j, k))/2*(U_procs(i, j, k)+U_procs(i, j, k+1))/2)/dZ
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Ay_procs(i, j, k) = (-(U_procs(i-1, j, k)+U_procs(i-1, j+1, k))/2*(V_procs(i-1, j, k)+V_procs(i, j, k))/2 &
                                       +(U_procs(i, j, k)+U_procs(i, j+1, k))/2*(V_procs(i, j, k)+V_procs(i+1, j, k))/2)/dX &
                                       +(-((V_procs(i, j-1, k)+V_procs(i, j, k))/2)**2 + ((V_procs(i, j, k)+V_procs(i, j+1, k))/2)**2)/dY &
                                       +(-(W_procs(i, j, k-1)+W_procs(i, j+1, k-1))/2*(V_procs(i, j, k-1)+V_procs(i, j, k))/2 &
                                       +(W_procs(i, j, k)+W_procs(i, j+1, k))/2*(V_procs(i, j, k)+V_procs(i, j, k+1))/2)/dZ
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Az_procs(i, j, k) = (-(U_procs(i-1, j, k)+U_procs(i-1, j, k+1))/2*(W_procs(i-1, j, k)+W_procs(i, j, k))/2 &
                                       +(U_procs(i, j, k)+U_procs(i, j, k+1))/2*(W_procs(i, j, k)+W_procs(i+1, j, k))/2)/dX &
                                       +(-(V_procs(i, j-1, k)+V_procs(i, j-1, k+1))/2*(W_procs(i, j-1, k)+W_procs(i, j, k))/2 &
                                       +(V_procs(i, j, k)+V_procs(i, j, k+1))/2*(W_procs(i, j, k)+W_procs(i, j+1, k))/2)/dY &
                                       +(-((W_procs(i, j, k-1)+W_procs(i, j, k))/2)**2 + ((W_procs(i, j, k)+W_procs(i, j, k+1))/2)**2)/dZ
                enddo
            enddo
        enddo

        Ax_procs(:, :, :) = -Ax_procs(:, :, :)
        Ay_procs(:, :, :) = -Ay_procs(:, :, :)
        Az_procs(:, :, :) = -Az_procs(:, :, :)

    end subroutine convection

    subroutine viscous(U_procs, V_procs, W_procs, Bx_procs, By_procs, Bz_procs)  ! ニュー一定
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        integer i, j, k
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Bx_procs(i, j, k) = (U_procs(i-1, j, k) - 2*U_procs(i, j, k) + U_procs(i+1, j, k)) / dX**2 &
                                       +(U_procs(i, j-1, k) - 2*U_procs(i, j, k) + U_procs(i, j+1, k)) / dY**2 &
                                       +(U_procs(i, j, k-1) - 2*U_procs(i, j, k) + U_procs(i, j, k+1)) / dZ**2
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    By_procs(i, j, k) = (V_procs(i-1, j, k) - 2*V_procs(i, j, k) + V_procs(i+1, j, k)) / dX**2 &
                                       +(V_procs(i, j-1, k) - 2*V_procs(i, j, k) + V_procs(i, j+1, k)) / dY**2 &
                                       +(V_procs(i, j, k-1) - 2*V_procs(i, j, k) + V_procs(i, j, k+1)) / dZ**2
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Bz_procs(i, j, k) = (W_procs(i-1, j, k) - 2*W_procs(i, j, k) + W_procs(i+1, j, k)) / dX**2 &
                                       +(W_procs(i, j-1, k) - 2*W_procs(i, j, k) + W_procs(i, j+1, k)) / dY**2 &
                                       +(W_procs(i, j, k-1) - 2*W_procs(i, j, k) + W_procs(i, j, k+1)) / dZ**2
                enddo
            enddo
        enddo
        Bx_procs(:, :, :) = 1.0d0/Re*Bx_procs(:, :, :)
        By_procs(:, :, :) = 1.0d0/Re*By_procs(:, :, :)
        Bz_procs(:, :, :) = 1.0d0/Re*Bz_procs(:, :, :)

    end subroutine viscous

    subroutine navier(U_procs, V_procs, W_procs, P_procs, Up_procs, Vp_procs, Wp_procs, Ax_procs, Ay_procs, Az_procs, &
                      Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, Bx0_procs, By0_procs, Bz0_procs, &
                      Fx_procs, Fy_procs, Fz_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Bx0_procs(1:NX, 1:NY, 1:N_procs), By0_procs(1:NX, 1:NY, 1:N_procs), Bz0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Fx_procs(1:NX, 1:NY, 1:N_procs), Fy_procs(1:NX, 1:NY, 1:N_procs), Fz_procs(1:NX, 1:NY, 1:N_procs)
        integer, intent(in) :: step
        integer i, j, k

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0_procs(:, :, :) = Ax_procs(:, :, :)
            Ay0_procs(:, :, :) = Ay_procs(:, :, :)
            Az0_procs(:, :, :) = Az_procs(:, :, :)
            Bx0_procs(:, :, :) = Bx_procs(:, :, :)
            By0_procs(:, :, :) = By_procs(:, :, :)
            Bz0_procs(:, :, :) = Bz_procs(:, :, :)
        endif

        ! NS方程式で速度の予測
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Up_procs(i, j, k) = U_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i+1, j, k))/dX &
                                       +dt*(3.0d0*(Ax_procs(i, j, k) + Bx_procs(i, j, k)) - (Ax0_procs(i, j, k) + Bx0_procs(i, j, k)))/2.0d0 &
                                       +dt*Fx_procs(i, j, k)
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Vp_procs(i, j, k) = V_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j+1, k))/dY &
                                       +dt*(3.0d0*(Ay_procs(i, j, k) + By_procs(i, j, k)) - (Ay0_procs(i, j, k) + By0_procs(i, j, k)))/2.0d0 &
                                       +dt*Fy_procs(i, j, k)
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Wp_procs(i, j, k) = W_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j, k+1))/dZ &
                                        +dt*(3.0d0*(Az_procs(i, j, k) + Bz_procs(i, j, k)) - (Az0_procs(i, j, k) + Bz0_procs(i, j, k)))/2.0d0 &
                                        +dt*Fz_procs(i, j, k)
                enddo
            enddo
        enddo

        call MPI_Boundary(Up_procs)
        call MPI_Boundary(Vp_procs)
        call MPI_Boundary(Wp_procs)

        ! Phiを求めるときに0番目の値も必要
        call PBM(Up_procs)
        call PBM(Vp_procs)
        call PBM(Wp_procs)

        ! n-1ステップ目の保存
        Ax0_procs(:, :, :) = Ax_procs(:, :, :)
        Ay0_procs(:, :, :) = Ay_procs(:, :, :)
        Az0_procs(:, :, :) = Az_procs(:, :, :)
        Bx0_procs(:, :, :) = Bx_procs(:, :, :)
        By0_procs(:, :, :) = By_procs(:, :, :)
        Bz0_procs(:, :, :) = Bz_procs(:, :, :)
    end subroutine navier


    subroutine poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs, step)
        real(8), intent(in) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(inout) :: Phi_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        real(8), parameter :: SOR = 1.7d0
        real(8), parameter :: eps = 1.0d-5
        integer, parameter :: itrmax = 10000
        integer i, j, k, l, itr
        real(8) BXM, BXP, BYM, BYP, BZM, BZP, B0
        real(8) er, er0, er_sum, er0_sum
        real(8) E(1:NX, 1:NY, 1:N_procs), Q(1:NX, 1:NY, 1:N_procs)
        
        ! ポアソン方程式で用いる定数の計算(今回は等間隔な格子)
        BXM = 1.0d0 / dX**2
        BXP = 1.0d0 / dX**2
        BYM = 1.0d0 / dY**2
        BYP = 1.0d0 / dY**2
        BZM = 1.0d0 / dZ**2
        BZP = 1.0d0 / dZ**2
        B0 = BXM + BXP + BYM + BYP + BZM + BZP
        ! 右辺の計算
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q(i, j, k) = ((-Up_procs(i-1, j, k) + Up_procs(i, j, k)) / dX &
                                 +(-Vp_procs(i, j-1, k) + Vp_procs(i, j, k)) / dY &
                                 +(-Wp_procs(i, j, k-1) + Wp_procs(i, j, k)) / dZ) / dt
                enddo
            enddo
        enddo
        er0 = sum(Q**2)
        call MPI_Allreduce(er0, er0_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (abs(er0_sum) < 1.0d-16) then  ! マシン零以下
            er0_sum = er0_sum + 1.0d-16
        endif

        ! SOR法
        ! Phi_procs(:, :, :) = 0.0d0
        do itr = 1, itrmax
            do k = 1, N_procs
                do j = 1, NY
                    do l = 1, 2  ! 奇数と偶数に分けて計算する
                        do i = l, NX, 2
                            E(i, j, k) = BZM*Phi_procs(i, j, k-1) + BYM*Phi_procs(i, j-1, k) + BXM*Phi_procs(i-1, j, k) - B0*Phi_procs(i, j, k) &
                                        +BXP*Phi_procs(i+1, j, k) + BYP*Phi_procs(i, j+1, k) + BZP*Phi_procs(i, j, k+1) - Q(i, j, k)
                            Phi_procs(i, j, k) = Phi_procs(i, j, k) + SOR * E(i, j, k) / B0
                        enddo
                    enddo
                enddo
            enddo

            call MPI_Boundary(Phi_procs)
            call PBM(Phi_procs)
            
            er = sum(E**2)
            call MPI_Allreduce(er, er_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

            if (sqrt(er_sum / er0_sum) < eps) then
                exit
            endif
        enddo

        er0 = sum(Phi_procs(1:NX, 1:NY, 1:N_procs))/(NX*NY*NZ)
        call MPI_Allreduce(er0, er, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        Phi_procs(:, :, :) = Phi_procs(:, :, :) - er

        er0 = sum(Phi_procs(1:NX, 1:NY, 1:N_procs))/(NX*NY*NZ)
        call MPI_Allreduce(er0, er, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        
        if (myrank == 0) then
            open(20, file = './mpi_dat/debag_smac.d', position='append')
            write(20, '(a, I6, a, I6, a, E12.4, a, E12.4, a, E12.4)') 'step:', step, '  itr:', itr, '  er_sum:', er_sum, '  er0_sum:', er0_sum, '  Phi_sum:', er
            close(20)
        endif
    end subroutine poisson

    subroutine march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
        real(8), intent(in) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Phi_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(inout) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    U_procs(i, j, k) = Up_procs(i, j, k) - dt*(-Phi_procs(i, j, k)+Phi_procs(i+1, j, k))/dX
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    V_procs(i, j, k) = Vp_procs(i, j, k) - dt*(-Phi_procs(i, j, k)+Phi_procs(i, j+1, k))/dY
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    W_procs(i, j, k) = Wp_procs(i, j, k) - dt*(-Phi_procs(i, j, k)+Phi_procs(i, j, k+1))/dZ
                enddo
            enddo
        enddo

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    P_procs(i, j, k) = P_procs(i, j, k) + Phi_procs(i, j, k)
                enddo
            enddo
        enddo
        call MPI_Boundary(U_procs)
        call MPI_Boundary(V_procs)
        call MPI_Boundary(W_procs)
        call MPI_Boundary(P_procs)

        call PBM(U_procs)
        call PBM(V_procs)
        call PBM(W_procs)
        call PBM(P_procs)
    end subroutine march

    subroutine taylor_debag(U_procs, V_procs, W_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        integer i, j, k
        real(8) err
        real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) U0(0:NX+1, 0:NY+1, 0:NZ+1), V0(0:NX+1, 0:NY+1, 0:NZ+1), W0(0:NX+1, 0:NY+1, 0:NZ+1)
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        call MPI_Gather(U_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, U(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(V_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, V(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(W_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, W(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        
        if (myrank == 0) then
            U0(:, :, :) = 0.0d0
            V0(:, :, :) = 0.0d0
            W0(:, :, :) = 0.0d0
            do k = 1, NZ
                do j = 1, NY
                    do i = 1, NX
                        ! U0(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)
                        ! V0(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)
                        U0(i, j, k) = -sin(i*dX) * cos((k-0.5d0)*dZ)
                        W0(i, j, k) = cos((i-0.5d0)*dX) * sin(k*dZ)
                    enddo
                enddo
            enddo
            U0(:, :, :) = U0(:, :, :) * exp(-2.0d0/Re*step*dt)  ! 解析解
            V0(:, :, :) = V0(:, :, :) * exp(-2.0d0/Re*step*dt)
            W0(:, :, :) = W0(:, :, :) * exp(-2.0d0/Re*step*dt)

            ! i = NX/4
            ! j = NY/4
            ! k = NZ/4
            err = (sum((U0(1:NX, 1:NY, 1:NZ)-U(1:NX, 1:NY, 1:NZ))**2) &
                  +sum((V0(1:NX, 1:NY, 1:NZ)-V(1:NX, 1:NY, 1:NZ))**2) &
                  +sum((W0(1:NX, 1:NY, 1:NZ)-W(1:NX, 1:NY, 1:NZ))**2))/(NX*NY*NZ)
            open(10, file=dir//'taylor_debag.d', position='append')
            ! write(10, *) step, U0(i, j, k), U(i, j, k), U0(i, j, k) - U(i, j, k), sqrt(err)
            write(10, *) step, sqrt(err)
            close(10)
        endif
    end subroutine taylor_debag
    
    
    subroutine get_data(U_procs, V_procs, W_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        integer i, j, k
        real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) D(3, 3), Omega(3), S(3), Qti, E(3)
        character(8) chstep
        
        call MPI_Gather(U_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, U(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(V_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, V(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(W_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, W(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        U(:, :, 0) = U(:, :, NZ)
        U(:, :, NZ+1) = U(:, :, 1)
        V(:, :, 0) = V(:, :, NZ)
        V(:, :, NZ+1) = V(:, :, 1)
        W(:, :, 0) = W(:, :, NZ)
        W(:, :, NZ+1) = W(:, :, 1)

        if (myrank == 0) then
            write(chstep, '(I8.8)') step
            open(30, file = dir//'all_'//chstep//'.d')
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
                        E(1) = (U(i-1, j, k)+U(i, j, k))/2*U_C
                        E(2) = (V(i, j-1, k)+V(i, j, k))/2*U_C
                        E(3) = (W(i, j, k-1)+W(i, j, k))/2*U_C
                        Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                        write(30, '(11e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
                                               E(1), E(2), E(3), &
                                               Omega(3), sum(Omega**2)/2, sum(S**2)/2, sum(E**2)/2, Qti
                    enddo
                enddo
            enddo
            close(30)


            open(10, file = dir//'z_'//chstep//'.d')
            k = NZ/2  ! x-y平面
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
                    E(1) = (U(i-1, j, k)+U(i, j, k))/2*U_C
                    E(2) = (V(i, j-1, k)+V(i, j, k))/2*U_C
                    E(3) = (W(i, j, k-1)+W(i, j, k))/2*U_C
                    Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                    write(10, '(11e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
                                           E(1), E(2), E(3), &
                                           Omega(3), sum(Omega**2)/2, sum(S**2)/2, sum(E**2)/2, Qti
                enddo
            enddo
            close(10)


            open(20, file = dir//'y_'//chstep//'.d')
            j = NY/2  ! x-z平面
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
                    S(1) = (D(3, 2) + D(2, 3))
                    S(2) = (D(1, 3) + D(3, 1))
                    S(3) = (D(2, 1) + D(1, 2))
                    E(1) = (U(i-1, j, k)+U(i, j, k))/2*U_C
                    E(2) = (V(i, j-1, k)+V(i, j, k))/2*U_C
                    E(3) = (W(i, j, k-1)+W(i, j, k))/2*U_C
                    Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                    write(20, '(11e12.4)') (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, &
                                           E(1), E(2), E(3), &
                                           Omega(2), sum(Omega**2)/2, sum(S**2)/2, sum(E**2)/2, Qti
                enddo
            enddo
            close(20)
        endif
    end subroutine get_data

    subroutine logging(U_procs, V_procs, W_procs, step, time1)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        real(8), intent(in) :: time1
        real(8) K_energy, K_energy_sum
        real(8) time0
        integer time2, hour, min, sec
        integer i, j, k 
        real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) U_tmp(0:NX+1, 0:NY+1, 0:NZ+1), V_tmp(0:NX+1, 0:NY+1, 0:NZ+1), W_tmp(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) U_grad(NX, NY, NZ), V_grad(NX, NY, NZ), W_grad(NX, NY, NZ)
        real(8) U_rms, lamda, Re_lamda, tmp, tmp_sum, D(3, 3), S(3, 3), epsilon, eta

        if (mod(step, 10) == 0) then
            call MPI_Gather(U_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, U(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_Gather(V_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, V(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            call MPI_Gather(W_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, W(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
            if (myrank == 0) then
                U(:, :, 0) = U(:, :, NZ)
                U(:, :, NZ+1) = U(:, :, 1)
                V(:, :, 0) = V(:, :, NZ)
                V(:, :, NZ+1) = V(:, :, 1)
                W(:, :, 0) = W(:, :, NZ)
                W(:, :, NZ+1) = W(:, :, 1)

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
            endif

            ! コルモゴロフ長
            tmp = 0.0d0
            do k = 1, N_procs
                do j = 1, NY
                    do i = 1, NX
                        D(1, 1) = (U_procs(i, j, k) - U_procs(i-1, j, k))/dX
                        D(1, 2) = (U_procs(i, j+1, k) - U_procs(i, j-1, k) + U_procs(i-1, j+1, k) - U_procs(i-1, j-1, k))/(4*dY)
                        D(1, 3) = (U_procs(i, j, k+1) - U_procs(i, j, k-1) + U_procs(i-1, j, k+1) - U_procs(i-1, j, k-1))/(4*dZ)
                        D(2, 1) = (V_procs(i+1, j, k) - V_procs(i-1, j, k) + V_procs(i+1, j-1, k) - V_procs(i-1, j-1, k))/(4*dX)
                        D(2, 2) = (V_procs(i, j, k) - V_procs(i, j-1, k))/dY
                        D(2, 3) = (V_procs(i, j, k+1) - V_procs(i, j, k-1) + V_procs(i, j-1, k+1) - V_procs(i, j-1, k-1))/(4*dZ)
                        D(3, 1) = (W_procs(i+1, j, k) - W_procs(i-1, j, k) + W_procs(i+1, j, k-1) - W_procs(i-1, j, k-1))/(4*dX)
                        D(3, 2) = (W_procs(i, j+1, k) - W_procs(i, j-1, k) + W_procs(i, j+1, k-1) - W_procs(i, j-1, k-1))/(4*dY)
                        D(3, 3) = (W_procs(i, j, k) - W_procs(i, j, k-1))/dZ
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
            call MPI_Allreduce(tmp, tmp_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (myrank == 0) then
                epsilon = 0.5d0*nu*tmp_sum/(NX*NY*NZ)
                eta = (nu**3/epsilon)**0.25d0
                ! write(*, '(12e12.4)') Re_lamda, epsilon, eta, dX
            endif
        endif
        
        if (mod(step, 10) == 0) then
            K_energy = sum(U_procs(1:NX, 1:NY, 1:N_procs)**2) + sum(V_procs(1:NX, 1:NY, 1:N_procs)**2) + sum(W_procs(1:NX, 1:NY, 1:N_procs)**2)
            K_energy = K_energy*U_C**2/2
            K_energy = K_energy/(NX*NY*NZ)
            call MPI_Allreduce(K_energy, K_energy_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

            if (myrank == 0 .and. K_energy_sum*0.0d0 /= 0.0d0) stop 'NaN value'  ! NaNの判定

            if (myrank == 0) then
                open(30, file = dir//'debag_energy2.d', position='append')
                write(30, *) step*dt_C, K_energy_sum, Re_lamda, epsilon, eta
                close(30)
            endif
        endif

        if (mod(step, Dstep) == 0) then
            if (myrank == 0) then
                write(*, '(a, I7)', advance='no') 'step:', step
                write(*, '(a, e12.4)', advance='no') '  | K_energy:', K_energy_sum

                time0 = MPI_Wtime()
                time2 = int((time0-time1)*(Nstep-step)/step)
                hour = time2 / 3600
                min = mod(time2, 3600) / 60
                sec = mod(time2, 60)
                write(*, '(a, I3, a, I2, a, I2)', advance='no') '  | time_left:', hour, ':', min, ':', sec

                write(*, *) ''
            endif
        endif
    end subroutine logging

    subroutine input(U_procs, V_procs, W_procs, P_procs)  ! initを実行した後に書く
        real(8), intent(out) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1), P(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        character(8) str
        write(str, '(I8.8)') input_step  ! 数値を文字列に変換

        if (myrank == 0) then
            open(10, file=dir//str//'.d')
            do k = 1, NZ
                do j = 1, NY
                    do i = 1, NX
                        read(10, *) U(i, j, k), V(i, j, k), W(i, j, k), P(i, j, k)
                    enddo
                enddo
            enddo
            close(10)
        endif

        call MPI_Scatter(U(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, U_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(V(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, V_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(W(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, W_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Scatter(P(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, P_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        call MPI_Boundary(U_procs)
        call MPI_Boundary(V_procs)
        call MPI_Boundary(W_procs)
        call MPI_Boundary(P_procs)
        call PBM(U_procs)
        call PBM(V_procs)
        call PBM(W_procs)
        call PBM(P_procs)
        if (myrank == 0) write(*, *) 'input_file='//dir//str//'.d'
    end subroutine input

    subroutine output(U_procs, V_procs, W_procs, P_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1), P(0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        character(8) str
        write(str, '(I8.8)') step  ! 数値を文字列に変換

        call MPI_Gather(U_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, U(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(V_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, V(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(W_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, W(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(P_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, P(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        if (myrank == 0) then
            open(10, file=dir//str//'.d')
            do k = 1, NZ
                do j = 1, NY
                    do i = 1, NX
                        write(10, *) U(i, j, k), V(i, j, k), W(i, j, k), P(i, j, k)
                    enddo
                enddo
            enddo
            close(10)
        endif
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
    real(8) Re1(1:NX, 1:NY)
    complex(8) Im1(1:NX/2+1, 1:NY), Im2(1:NZ), Im3(1:NZ)
contains
    subroutine fft_init
        call dfftw_plan_dft_r2c_2d(plan1, NX, NY, Re1, Im1, FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(plan2, NZ, Im2, Im3, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(plan3, NZ, Im3, Im2, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(plan4, NX, NY, Im1, Re1, FFTW_ESTIMATE)
    end subroutine fft_init

    subroutine fftr2c_2d(input, output)
        real(8), intent(in) :: input(1:NX, 1:NY)
        complex(8), intent(out) :: output(1:NX/2+1, 1:NY)
        Re1(:, :) = input(:, :)
        call dfftw_execute(plan1, Re1, Im1)
        output(:, :) = Im1(:, :)
    end subroutine fftr2c_2d

    subroutine fftc2c_1d_forward(input, output)
        complex(8), intent(in) :: input(1:NZ)
        complex(8), intent(out) :: output(1:NZ)
        Im2(:) = input(:)
        call dfftw_execute(plan2, Im2, Im3)
        output(:) = Im3(:)
    end subroutine fftc2c_1d_forward

    subroutine fftc2c_1d_backward(input, output)
        complex(8), intent(in) :: input(1:NZ)
        complex(8), intent(out) :: output(1:NZ)
        Im3(:) = input(:)
        call dfftw_execute(plan3, Im3, Im2)
        output(:) = Im2(:)
    end subroutine fftc2c_1d_backward
    
    subroutine fftc2r_2d(input, output)
        complex(8), intent(in) :: input(1:NX/2+1, 1:NY)
        real(8), intent(out) :: output(1:NX, 1:NY)
        Im1(:, :) = input(:, :)
        call dfftw_execute(plan4, Im1, Re1)
        output(:, :) = Re1(:, :)
    end subroutine fftc2r_2d

    subroutine fft_finalize
        call dfftw_destroy_plan(plan1)
        call dfftw_destroy_plan(plan2)
        call dfftw_destroy_plan(plan3)
        call dfftw_destroy_plan(plan4)
    end subroutine fft_finalize

    subroutine fft_solve(Q_procs, LHS_procs, Phi_procs)
        real(8), intent(in) :: Q_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: LHS_procs(1:NX/2+1, 1:NY/procs, 1:NZ)
        real(8), intent(out) :: Phi_procs(1:NX, 1:NY, 1:N_procs)
        complex(8), allocatable :: Q_hat_procs(:, :, :), Phi_hat_procs(:, :, :)
        complex(8), allocatable :: Q_hat_hat_procs(:, :, :), Phi_hat_hat_procs(:, :, :)
        complex(8), allocatable :: Z_procs(:, :, :)
        complex(8), allocatable :: T1_procs(:, :, :, :), T2_procs(:, :, :, :)
        integer i, j, k
        integer NY_procs
        call fft_init

        NY_procs = NY / procs
        allocate(Q_hat_procs(1:NX/2+1, 1:NY, 1:N_procs), Phi_hat_procs(1:NX/2+1, 1:NY, 1:N_procs))
        allocate(Q_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ), Phi_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ))
        allocate(Z_procs(1:NX/2+1, 1:NY_procs, 1:NZ))
        allocate(T1_procs(1:NX/2+1, 1:NY_procs, 1:procs, 1:N_procs), T2_procs(1:NX/2+1, 1:NY_procs, 1:N_procs, 1:procs))


        ! 右辺をx,y方向にfft
        do k = 1, N_procs
            call fftr2c_2d(Q_procs(:, :, k), Q_hat_procs(:, :, k))
        enddo

        ! 軸変換
        T1_procs(:, :, :, :) = reshape(Q_hat_procs, shape(T1_procs))
        T2_procs(:, :, :, :) = reshape(T1_procs, shape(T2_procs), order=[1, 2, 4, 3])
        call MPI_Alltoall(T2_procs(1, 1, 1, 1), (NX/2+1)*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, &
                          Z_procs(1, 1, 1), (NX/2+1)*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        ! z方向にfft
        do j = 1, NY_procs
            do i = 1, NX/2+1
                call fftc2c_1d_forward(Z_procs(i, j, :), Q_hat_hat_procs(i, j, :))
            enddo
        enddo
        Q_hat_hat_procs(:, :, :) = Q_hat_hat_procs(:, :, :)/(NX*NY*NZ)

        ! 定数で割ることで左辺を求める
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    if (abs(LHS_procs(i, j, k)) < 1.0d-16) then  ! マシン零以下
                        Phi_hat_hat_procs(i, j, k) = (0.0d0, 0.0d0)
                    else
                        Phi_hat_hat_procs(i, j, k) = Q_hat_hat_procs(i, j, k) / LHS_procs(i, j, k)
                    endif
                enddo
            enddo
        enddo

        ! 左辺をz方向に逆fft
        do j = 1, NY_procs
            do i = 1, NX/2+1
                call fftc2c_1d_backward(Phi_hat_hat_procs(i, j, :), Z_procs(i, j, :))
            enddo
        enddo

        ! 軸変換
        call MPI_Alltoall(Z_procs(1, 1, 1), (NX/2+1)*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, &
                          T2_procs(1, 1, 1, 1), (NX/2+1)*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
        T1_procs(:, :, :, :) = reshape(T2_procs, shape(T1_procs), order=[1, 2, 4, 3])
        Phi_hat_procs(:, :, :) = reshape(T1_procs, shape(Phi_hat_procs))

        ! x,y方向に逆fft
        do k = 1, N_procs
            call fftc2r_2d(Phi_hat_procs(:, :, k),  Phi_procs(:, :, k))
        enddo

        call fft_finalize

    end subroutine fft_solve

    subroutine fft_navier(U_procs, V_procs, W_procs, P_procs, Up_procs, Vp_procs, Wp_procs, Ax_procs, Ay_procs, Az_procs, &
                          Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, Fx_procs, Fy_procs, Fz_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Fx_procs(1:NX, 1:NY, 1:N_procs), Fy_procs(1:NX, 1:NY, 1:N_procs), Fz_procs(1:NX, 1:NY, 1:N_procs)
        integer, intent(in) :: step
        real(8) Q_procs(1:NX, 1:NY, 1:N_procs)
        real(8) LHS_procs(1:NX/2+1, 1:NY/procs, 1:NZ)
        integer i, j, k
        integer NY_procs
        NY_procs = NY / procs

        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0_procs(:, :, :) = Ax_procs(:, :, :)
            Ay0_procs(:, :, :) = Ay_procs(:, :, :)
            Az0_procs(:, :, :) = Az_procs(:, :, :)
        endif

        ! 左辺の定数部分の計算
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    ! LHS_procs(i, j, k) = 1.0d0 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX*2*PI))/dX**2 &
                    !                                     + 2*(1-cos((myrank*NY_procs + j-1)*dY*2*PI))/dY**2 &
                    !                                     + 2*(1-cos((k-1)*dZ*2*PI))/dZ**2)
                    LHS_procs(i, j, k) = 1.0d0 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                                        + 2*(1-cos((myrank*NY_procs + j-1)*dY))/dY**2 &
                                                        + 2*(1-cos((k-1)*dZ))/dZ**2)
                enddo
            enddo
        enddo

        ! NS方程式で速度の予測
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = U_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i+1, j, k))/dX &
                                      +dt*(3.0d0*Ax_procs(i, j, k) - Ax0_procs(i, j, k))/2.0d0 &
                                      +dt*Bx_procs(i, j, k)/2.0d0 &
                                      +dt*Fx_procs(i, j, k)
                enddo
            enddo
        enddo

        call fft_solve(Q_procs, LHS_procs, Up_procs(1:NX, 1:NY, 1:N_procs))

        
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = V_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j+1, k))/dY &
                                      +dt*(3.0d0*Ay_procs(i, j, k) - Ay0_procs(i, j, k))/2.0d0 &
                                      +dt*By_procs(i, j, k)/2.0d0 &
                                      +dt*Fy_procs(i, j, k)
                enddo
            enddo
        enddo

        call fft_solve(Q_procs, LHS_procs, Vp_procs(1:NX, 1:NY, 1:N_procs))

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = W_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j, k+1))/dZ &
                                       +dt*(3.0d0*Az_procs(i, j, k) - Az0_procs(i, j, k))/2.0d0 &
                                       +dt*Bz_procs(i, j, k)/2.0d0 &
                                       +dt*Fz_procs(i, j, k)
                enddo
            enddo
        enddo

        call fft_solve(Q_procs, LHS_procs, Wp_procs(1:NX, 1:NY, 1:N_procs))

        call MPI_Boundary(Up_procs)
        call MPI_Boundary(Vp_procs)
        call MPI_Boundary(Wp_procs)

        ! Phiを求めるときに0番目の値も必要
        call PBM(Up_procs)
        call PBM(Vp_procs)
        call PBM(Wp_procs)

        ! n-1ステップ目の保存
        Ax0_procs(:, :, :) = Ax_procs(:, :, :)
        Ay0_procs(:, :, :) = Ay_procs(:, :, :)
        Az0_procs(:, :, :) = Az_procs(:, :, :)
    end subroutine fft_navier

    subroutine fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
        real(8), intent(in) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Phi_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8) Q_procs(1:NX, 1:NY, 1:N_procs)
        real(8) LHS_procs(1:NX/2+1, 1:NY/procs, 1:NZ)
        integer i, j, k
        integer NY_procs

        NY_procs = NY / procs

        ! 右辺の計算
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = ((-Up_procs(i-1, j, k) + Up_procs(i, j, k)) / dX &
                                       +(-Vp_procs(i, j-1, k) + Vp_procs(i, j, k)) / dY &
                                       +(-Wp_procs(i, j, k-1) + Wp_procs(i, j, k)) / dZ) / dt
                enddo
            enddo
        enddo

        ! 左辺の定数の計算
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    LHS_procs(i, j, k) = -(2*(1-cos((i-1)*dX))/dX**2 &
                                         + 2*(1-cos((myrank*NY_procs + j-1)*dY))/dY**2 &
                                         + 2*(1-cos((k-1)*dZ))/dZ**2)
                enddo
            enddo
        enddo

        ! FFT
        call fft_solve(Q_procs, LHS_procs, Phi_procs(1:NX, 1:NY, 1:N_procs))

        call MPI_Boundary(Phi_procs)
        call PBM(Phi_procs)

    end subroutine fft_poisson

    subroutine fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
        real(8), intent(in) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Phi_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(inout) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    U_procs(i, j, k) = Up_procs(i, j, k) - dt*(-Phi_procs(i, j, k)+Phi_procs(i+1, j, k))/dX
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    V_procs(i, j, k) = Vp_procs(i, j, k) - dt*(-Phi_procs(i, j, k)+Phi_procs(i, j+1, k))/dY
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    W_procs(i, j, k) = Wp_procs(i, j, k) - dt*(-Phi_procs(i, j, k)+Phi_procs(i, j, k+1))/dZ
                enddo
            enddo
        enddo

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    P_procs(i, j, k) = P_procs(i, j, k) + Phi_procs(i, j, k) &
                                      -dt/Re/2.0d0*((Phi_procs(i-1, j, k)-2*Phi_procs(i, j, k)+Phi_procs(i+1, j, k))/dX**2 &
                                                   +(Phi_procs(i, j-1, k)-2*Phi_procs(i, j, k)+Phi_procs(i, j+1, k))/dY**2 &
                                                   +(Phi_procs(i, j, k-1)-2*Phi_procs(i, j, k)+Phi_procs(i, j, k+1))/dZ**2)
                enddo
            enddo
        enddo
        
        call MPI_Boundary(U_procs)
        call MPI_Boundary(V_procs)
        call MPI_Boundary(W_procs)
        call MPI_Boundary(P_procs)

        call PBM(U_procs)
        call PBM(V_procs)
        call PBM(W_procs)
        call PBM(P_procs)
    end subroutine fft_march


    subroutine fft_forward(Q_procs, Q_hat_hat_procs)
        real(8), intent(in) :: Q_procs(1:NX, 1:NY, 1:N_procs)
        complex(8), allocatable :: Q_hat_procs(:, :, :)
        complex(8), allocatable :: Q_hat_hat_procs(:, :, :)
        complex(8), allocatable :: Z_procs(:, :, :)
        complex(8), allocatable :: T1_procs(:, :, :, :), T2_procs(:, :, :, :)
        integer i, j, k
        integer NY_procs
        call fft_init

        NY_procs = NY / procs
        allocate(Q_hat_procs(1:NX/2+1, 1:NY, 1:N_procs))
        allocate(Q_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ))
        allocate(Z_procs(1:NX/2+1, 1:NY_procs, 1:NZ))
        allocate(T1_procs(1:NX/2+1, 1:NY_procs, 1:procs, 1:N_procs), T2_procs(1:NX/2+1, 1:NY_procs, 1:N_procs, 1:procs))


        ! 右辺をx,y方向にfft
        do k = 1, N_procs
            call fftr2c_2d(Q_procs(:, :, k), Q_hat_procs(:, :, k))
        enddo

        ! 軸変換
        T1_procs(:, :, :, :) = reshape(Q_hat_procs, shape(T1_procs))
        T2_procs(:, :, :, :) = reshape(T1_procs, shape(T2_procs), order=[1, 2, 4, 3])
        call MPI_Alltoall(T2_procs(1, 1, 1, 1), (NX/2+1)*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, &
                          Z_procs(1, 1, 1), (NX/2+1)*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        ! z方向にfft
        do j = 1, NY_procs
            do i = 1, NX/2+1
                call fftc2c_1d_forward(Z_procs(i, j, :), Q_hat_hat_procs(i, j, :))
            enddo
        enddo
        Q_hat_hat_procs(:, :, :) = Q_hat_hat_procs(:, :, :)/(NX*NY*NZ)

        deallocate(Q_hat_procs, Z_procs, T1_procs, T2_procs)  ! 214000目で発生するエラー対策
        ! The program was terminated abnormally with signal number SIGBUS.
        ! signal identifier = BUS_ADRERR, non-existent physical address
        ! もしかしたら，W=0となっていたことが原因かも．

        call fft_finalize

    end subroutine fft_forward


    subroutine energy_sum(U_procs, V_procs, W_procs, Energy)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(inout) :: Energy(0:NX)
        complex(8), allocatable :: U_hat_hat_procs(:, :, :), V_hat_hat_procs(:, :, :), W_hat_hat_procs(:, :, :)
        real(8), allocatable :: E_tmp_procs(:, :, :)
        real(8), allocatable :: K_abs_procs(:, :, :)
        real(8) Energy_procs(0:NX), Energy_procs_sum(0:NX)
        integer i, j, k, index
        real(8) kx, ky, kz
        integer NY_procs
        NY_procs = NY / procs

        allocate(E_tmp_procs(1:NX/2+1, 1:NY_procs, 1:NZ), K_abs_procs(1:NX/2+1, 1:NY_procs, 1:NZ))

        ! 速度場をフーリエ変換
        call fft_forward(U_procs(1:NX, 1:NY, 1:N_procs)*U_C, U_hat_hat_procs)
        call fft_forward(V_procs(1:NX, 1:NY, 1:N_procs)*U_C, V_hat_hat_procs)
        call fft_forward(W_procs(1:NX, 1:NY, 1:N_procs)*U_C, W_hat_hat_procs)

        ! 配列要素に対するエネルギー
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    E_tmp_procs(i, j, k) = abs(U_hat_hat_procs(i, j, k))**2.0d0 &
                                         + abs(V_hat_hat_procs(i, j, k))**2.0d0 &
                                         + abs(W_hat_hat_procs(i, j, k))**2.0d0
                enddo
            enddo
        enddo

        ! 配列要素に対する波数
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    kx = sign(1.0d0, dble(NX + 3)/2.0d0 - dble(i)) * (- dble(abs(NX/2 + 1 - i)) + dble(NX)/2.0d0)  ! 増田さん参考
                    ky = sign(1.0d0, dble(NY + 3)/2.0d0 - dble(myrank*NY_procs + j)) * (- dble(abs(NY/2 + 1 - myrank*NY_procs - j)) + dble(NY)/2.0d0)
                    kz = sign(1.0d0, dble(NZ + 3)/2.0d0 - dble(k)) * (- dble(abs(NZ/2 + 1 - k)) + dble(NZ)/2.0d0)
                    K_abs_procs(i, j, k) = sqrt(kx**2.0d0 + ky**2.0d0 + kz**2.0d0)
                enddo
            enddo
        enddo

        ! 波数を四捨五入し、対応する整数の波数にエネルギーを足し合わせる
        Energy_procs(:) = 0.0d0
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    index = nint(K_abs_procs(i, j, k))
                    if (i==1 .or. i==NX/2+1) then
                        Energy_procs(index) = Energy_procs(index) + E_tmp_procs(i, j, k)/2.0d0
                    else
                        Energy_procs(index) = Energy_procs(index) + E_tmp_procs(i, j, k)/2.0d0*2.0d0
                    endif
                enddo
            enddo
        enddo

        call MPI_Allreduce(Energy_procs(0), Energy_procs_sum(0), (NX+1), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        Energy(:) = Energy(:) + Energy_procs_sum(:)

        deallocate(E_tmp_procs, K_abs_procs)
    end subroutine energy_sum

    subroutine energy_reset(U_procs, V_procs, W_procs, step, Energy)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        real(8), intent(inout) :: Energy(0:NX)
        real(8) K_energy, K_energy_sum
        integer i
        character(8) str
        write(str, '(I8.8)') step

        Energy(:) = Energy(:)/Estep

        if (myrank == 0) then
            open(30, file = dir//'energy_'//str//'.d')
            do i = 0, NX
                write(30, '(I4, e12.4)') i, Energy(i)
            enddo
            close(30)
        endif

        K_energy = sum(U_procs(1:NX, 1:NY, 1:N_procs)**2) + sum(V_procs(1:NX, 1:NY, 1:N_procs)**2) + sum(W_procs(1:NX, 1:NY, 1:N_procs)**2)
        K_energy = K_energy*U_C**2/2
        K_energy = K_energy/(NX*NY*NZ)
        call MPI_Allreduce(K_energy, K_energy_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (myrank == 0) then
            open(30, file = dir//'energy_check.d', position='append')  ! スペクトルは平均なので一致しない
            write(30, '(I8, 4e12.4)') step, K_energy_sum, sum(Energy(:)), K_energy_sum-sum(Energy(:)), K_energy_sum/sum(Energy(:))
            close(30)
        endif

        Energy(:) = 0.0d0

    end subroutine energy_reset
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
    subroutine ibm_init(X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, Uc_procs, Vc_procs, Wc_procs, &
                        Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
        real(8), allocatable :: X_procs(:, :, :), Y_procs(:, :, :), Z_procs(:, :, :)
        real(8), allocatable :: Xc_procs(:, :, :), Yc_procs(:, :, :), Zc_procs(:, :, :)
        real(8), allocatable :: Uc_procs(:, :, :), Vc_procs(:, :, :), Wc_procs(:, :, :)
        real(8), allocatable :: Ua_procs(:, :, :), Va_procs(:, :, :), Wa_procs(:, :, :)
        real(8), allocatable :: Fxc_procs(:, :, :), Fyc_procs(:, :, :), Fzc_procs(:, :, :)
        real(8), allocatable :: fxint_procs(:, :, :), fyint_procs(:, :, :), fzint_procs(:, :, :)
        integer i, j, k

        ! コメント
        if (dX /= dY .or. dX /= dZ) stop 'dX, dY and dZ are not equal'
        if (myrank == 0) then
            write(*, '(a, I5)') 'NL   :', NL
            write(*, '(a, E12.4)') 'dX**3:', dX**3
            write(*, '(a, E12.4)') 'dV   :', dV
        endif

        ! 配列のallocate
        N_procs = NZ / procs
        allocate(X_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Y_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Z_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Xc_procs(1:NC, 1:NL, 1:N_procs), Yc_procs(1:NC, 1:NL, 1:N_procs), Zc_procs(1:NC, 1:NL, 1:N_procs))
        allocate(Uc_procs(1:NC, 1:NL, 1:N_procs), Vc_procs(1:NC, 1:NL, 1:N_procs), Wc_procs(1:NC, 1:NL, 1:N_procs))
        allocate(Ua_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Va_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wa_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Fxc_procs(1:NX, 1:NY, 1:N_procs), Fyc_procs(1:NX, 1:NY, 1:N_procs), Fzc_procs(1:NX, 1:NY, 1:N_procs))
        allocate(fxint_procs(1:NX, 1:NY, 1:N_procs), fyint_procs(1:NX, 1:NY, 1:N_procs), fzint_procs(1:NX, 1:NY, 1:N_procs))

        ! 格子点中心座標
        do i = 0, NX+1
            X_procs(i, :, :) = (i-0.5d0)*dX
        enddo
        do j = 0, NY+1
            Y_procs(:, j, :) = (j-0.5d0)*dY
        enddo
        do k = 0, N_procs+1
            Z_procs(:, :, k) = (myrank*N_procs + k-0.5d0)*dZ
        enddo

        ! 円柱上の座標
        if (NC == 4) then
            do j = 1, NL
                Xc_procs(1, j, :) = 2*PI* 5/16 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(1, j, :) = 2*PI* 5/16 + DC/2 * sin(2*PI/NL*j)
                Xc_procs(2, j, :) = 2*PI*11/16 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(2, j, :) = 2*PI* 5/16 + DC/2 * sin(2*PI/NL*j)
                Xc_procs(3, j, :) = 2*PI* 5/16 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(3, j, :) = 2*PI*11/16 + DC/2 * sin(2*PI/NL*j)
                Xc_procs(4, j, :) = 2*PI*11/16 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(4, j, :) = 2*PI*11/16 + DC/2 * sin(2*PI/NL*j)
            enddo
        else if (NC == 1) then
            do j = 1, NL
                Xc_procs(1, j, :) = 2*PI/4 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(1, j, :) = 2*PI/2 + DC/2 * sin(2*PI/NL*j)
            enddo
        endif
        do k = 1, N_procs
            Zc_procs(:, :, k) = (myrank*N_procs + k-0.5d0)*dZ
        enddo

        ! 円柱上の座標での速度
        Uc_procs(:, :, :) = 0.0d0
        Vc_procs(:, :, :) = 0.0d0
        Wc_procs(:, :, :) = 0.0d0
        if (NC == 4) then
            do j = 1, NL
                Uc_procs(1, j, :) = -sin(2*PI/NL*j)
                Vc_procs(1, j, :) = cos(2*PI/NL*j)
                Uc_procs(2, j, :) = sin(2*PI/NL*j)
                Vc_procs(2, j, :) = -cos(2*PI/NL*j)
                Uc_procs(3, j, :) = sin(2*PI/NL*j)
                Vc_procs(3, j, :) = -cos(2*PI/NL*j)
                Uc_procs(4, j, :) = -sin(2*PI/NL*j)
                Vc_procs(4, j, :) = cos(2*PI/NL*j)
            enddo
        endif

    end subroutine ibm_init

    subroutine ibm_preliminary(U_procs, V_procs, W_procs, P_procs, Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                               Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, Up_procs, Vp_procs, Wp_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(out) :: Ua_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Va_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wa_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(inout) :: Fxc_procs(1:NX, 1:NY, 1:N_procs), Fyc_procs(1:NX, 1:NY, 1:N_procs), Fzc_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(out) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        integer i, j, k
        
        if (step==1) then  ! 1ステップ目のみ例外処理
            Ax0_procs(:, :, :) = Ax_procs(:, :, :)
            Ay0_procs(:, :, :) = Ay_procs(:, :, :)
            Az0_procs(:, :, :) = Az_procs(:, :, :)
        endif

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Ua_procs(i, j, k) = U_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i+1, j, k))/dX &
                                       +dt*(3.0d0*Ax_procs(i, j, k) - Ax0_procs(i, j, k))/2.0d0 &
                                       +dt*Bx_procs(i, j, k)

                    Va_procs(i, j, k) = V_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j+1, k))/dY &
                                       +dt*(3.0d0*Ay_procs(i, j, k) - Ay0_procs(i, j, k))/2.0d0 &
                                       +dt*By_procs(i, j, k)

                    Wa_procs(i, j, k) = W_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j, k+1))/dZ &
                                       +dt*(3.0d0*Az_procs(i, j, k) - Az0_procs(i, j, k))/2.0d0 &
                                       +dt*Bz_procs(i, j, k)
                enddo
            enddo
        enddo

        ! 周期境界
        call MPI_Boundary(Ua_procs)
        call MPI_Boundary(Va_procs)
        call MPI_Boundary(Wa_procs)
        call PBM(Ua_procs)
        call PBM(Va_procs)
        call PBM(Wa_procs)

        ! n-1ステップ目の保存
        Ax0_procs(:, :, :) = Ax_procs(:, :, :)
        Ay0_procs(:, :, :) = Ay_procs(:, :, :)
        Az0_procs(:, :, :) = Az_procs(:, :, :)

        ! 円柱上での外力初期化
        Fxc_procs(:, :, :) = 0.0d0
        Fyc_procs(:, :, :) = 0.0d0
        Fzc_procs(:, :, :) = 0.0d0

        ! 反復回数によらずHelmholtzの式を使うため
        Up_procs(:, :, :) = Ua_procs(:, :, :)
        Vp_procs(:, :, :) = Va_procs(:, :, :)
        Wp_procs(:, :, :) = Wa_procs(:, :, :)

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

    subroutine ibm_Helmholtz(Up_procs, Vp_procs, Wp_procs, X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs,&
                             Uc_procs, Vc_procs, Wc_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
        real(8), intent(in) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: X_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Y_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Z_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Xc_procs(1:NC, 1:NL, 1:N_procs), Yc_procs(1:NC, 1:NL, 1:N_procs), Zc_procs(1:NC, 1:NL, 1:N_procs)
        real(8), intent(in) :: Uc_procs(1:NC, 1:NL, 1:N_procs), Vc_procs(1:NC, 1:NL, 1:N_procs), Wc_procs(1:NC, 1:NL, 1:N_procs)
        real(8), intent(inout) :: Fxc_procs(1:NX, 1:NY, 1:N_procs), Fyc_procs(1:NX, 1:NY, 1:N_procs), Fzc_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(out) :: fxint_procs(1:NX, 1:NY, 1:N_procs), fyint_procs(1:NX, 1:NY, 1:N_procs), fzint_procs(1:NX, 1:NY, 1:N_procs)
        integer i, j, k, l, m, n
        real(8) Ub_procs(1:NC, 1:NL, 1:N_procs), Vb_procs(1:NC, 1:NL, 1:N_procs), Wb_procs(1:NC, 1:NL, 1:N_procs)
        real(8) fxtmp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), fytmp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), fztmp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        
        Ub_procs(:, :, :) = 0.0d0
        Vb_procs(:, :, :) = 0.0d0
        Wb_procs(:, :, :) = 0.0d0

        do n = 1, N_procs
            do m = 1, NL
                do l = 1, NC
                    ! 外力点周囲の3*3*3点の和をとる
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY), int(Yc_procs(l, m, n)/dY) + 2
                            do i = int(Xc_procs(l, m, n)/dX - 0.5d0), int(Xc_procs(l, m, n)/dX - 0.5d0) + 2
                                Ub_procs(l, m, n) = Ub_procs(l, m, n) &
                                            + Up_procs(i, j, k) * delta(X_procs(i, j, k)+0.5d0*dX, Y_procs(i, j, k), Z_procs(i, j, k), &
                                                                  Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY - 0.5d0), int(Yc_procs(l, m, n)/dY - 0.5d0) + 2
                            do i = int(Xc_procs(l, m, n)/dX), int(Xc_procs(l, m, n)/dX) + 2
                                Vb_procs(l, m, n) = Vb_procs(l, m, n) &
                                            + Vp_procs(i, j, k) * delta(X_procs(i, j, k), Y_procs(i, j, k)+0.5d0*dY, Z_procs(i, j, k), &
                                                                  Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY) , int(Yc_procs(l, m, n)/dY) + 2
                            do i = int(Xc_procs(l, m, n)/dX), int(Xc_procs(l, m, n)/dX) + 2
                                Wb_procs(l, m, n) = Wb_procs(l, m, n) &
                                            + Wp_procs(i, j, k) * delta(X_procs(i, j, k), Y_procs(i, j, k), Z_procs(i, j, k)+0.5d0*dZ, &
                                                                  Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! if (myrank == 0) write(*, *) sum(Ub_procs)/(N_procs*NL*NC)
        ! if (myrank == 0) write(*, *) sum(Vb_procs)/(N_procs*NL*NC)
        ! if (myrank == 0) write(*, *) sum(Wb_procs)/(N_procs*NL*NC)

        do n = 1, N_procs
            do m = 1, NL
                do l = 1, NC
                    Fxc_procs(l, m, n) = Fxc_procs(l, m, n) + (Uc_procs(l, m, n) - Ub_procs(l, m, n))/dt
                    Fyc_procs(l, m, n) = Fyc_procs(l, m, n) + (Vc_procs(l, m, n) - Vb_procs(l, m, n))/dt
                    Fzc_procs(l, m, n) = Fzc_procs(l, m, n) + (Wc_procs(l, m, n) - Wb_procs(l, m, n))/dt
                enddo
            enddo
        enddo

        fxtmp_procs(:, :, :) = 0.0d0
        fytmp_procs(:, :, :) = 0.0d0
        fztmp_procs(:, :, :) = 0.0d0
        do n = 1, N_procs
            do m = 1, NL
                do l = 1, NC
                    ! 格子点周辺の外力点の和をとる、つまり外力点周辺の格子点へ和をのこす
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY), int(Yc_procs(l, m, n)/dY) + 2
                            do i = int(Xc_procs(l, m, n)/dX - 0.5d0), int(Xc_procs(l, m, n)/dX - 0.5d0) + 2
                                fxtmp_procs(i, j, k) = fxtmp_procs(i, j, k) &
                                               + Fxc_procs(l, m, n) * delta(X_procs(i, j, k)+0.5d0*dX, Y_procs(i, j, k), Z_procs(i, j, k), &
                                                                      Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dV
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY - 0.5d0), int(Yc_procs(l, m, n)/dY - 0.5d0) + 2
                            do i = int(Xc_procs(l, m, n)/dX), int(Xc_procs(l, m, n)/dX) + 2
                                fytmp_procs(i, j, k) = fytmp_procs(i, j, k) &
                                               + Fyc_procs(l, m, n) * delta(X_procs(i, j, k), Y_procs(i, j, k)+0.5d0*dY, Z_procs(i, j, k), &
                                                                      Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dV
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY) , int(Yc_procs(l, m, n)/dY) + 2
                            do i = int(Xc_procs(l, m, n)/dX), int(Xc_procs(l, m, n)/dX) + 2
                                fztmp_procs(i, j, k) = fztmp_procs(i, j, k) &
                                               + Fzc_procs(l, m, n) * delta(X_procs(i, j, k), Y_procs(i, j, k), Z_procs(i, j, k)+0.5d0*dZ, &
                                                                      Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dV
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        fxint_procs(:, :, :) = fxtmp_procs(1:NX, 1:NY, 1:N_procs)
        fyint_procs(:, :, :) = fytmp_procs(1:NX, 1:NY, 1:N_procs)
        fzint_procs(:, :, :) = fztmp_procs(1:NX, 1:NY, 1:N_procs)
        ! 上書き
        fxtmp_procs(:, :, 1) = fxtmp_procs(:, :, 0)
        fytmp_procs(:, :, 1) = fytmp_procs(:, :, 0)
        fztmp_procs(:, :, 1) = fztmp_procs(:, :, 0)
        fxtmp_procs(:, :, N_procs) = fxtmp_procs(:, :, N_procs+1)
        fytmp_procs(:, :, N_procs) = fytmp_procs(:, :, N_procs+1)
        fztmp_procs(:, :, N_procs) = fztmp_procs(:, :, N_procs+1)
        ! 通信
        call MPI_Boundary(fxtmp_procs)
        call MPI_Boundary(fytmp_procs)
        call MPI_Boundary(fztmp_procs)
        ! 足し算
        fxint_procs(:, :, 1) = fxint_procs(:, :, 1) + fxtmp_procs(1:NX, 1:NY, 0)
        fyint_procs(:, :, 1) = fyint_procs(:, :, 1) + fytmp_procs(1:NX, 1:NY, 0)
        fzint_procs(:, :, 1) = fzint_procs(:, :, 1) + fztmp_procs(1:NX, 1:NY, 0)
        fxint_procs(:, :, N_procs) = fxint_procs(:, :, N_procs) + fxtmp_procs(1:NX, 1:NY, N_procs+1)
        fyint_procs(:, :, N_procs) = fyint_procs(:, :, N_procs) + fytmp_procs(1:NX, 1:NY, N_procs+1)
        fzint_procs(:, :, N_procs) = fzint_procs(:, :, N_procs) + fztmp_procs(1:NX, 1:NY, N_procs+1)

    end subroutine ibm_Helmholtz

    subroutine ibm_predict(Ua_procs, Va_procs, Wa_procs, fxint_procs, fyint_procs, fzint_procs, Bx_procs, By_procs, Bz_procs, Up_procs, Vp_procs, Wp_procs)
        real(8), intent(in) :: Ua_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Va_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wa_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: fxint_procs(1:NX, 1:NY, 1:N_procs), fyint_procs(1:NX, 1:NY, 1:N_procs), fzint_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(out) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k
        real(8) Q_procs(1:NX, 1:NY, 1:N_procs), LHS_procs(1:NX/2+1, 1:NY/procs, 1:NZ)
        integer NY_procs

        NY_procs = NY / procs

        ! 左辺の定数部分の計算
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    LHS_procs(i, j, k) = 1 + dt/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                                        + 2*(1-cos((myrank*NY_procs + j-1)*dY))/dY**2 &
                                                        + 2*(1-cos((k-1)*dZ))/dZ**2)
                enddo
            enddo
        enddo

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = Ua_procs(i, j, k) + dt*fxint_procs(i, j, k) - dt*Bx_procs(i, j, k)/2
                enddo
            enddo
        enddo
        call fft_solve(Q_procs, LHS_procs, Up_procs(1:NX, 1:NY, 1:N_procs))

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = Va_procs(i, j, k) + dt*fyint_procs(i, j, k) - dt*By_procs(i, j, k)/2
                enddo
            enddo
        enddo
        call fft_solve(Q_procs, LHS_procs, Vp_procs(1:NX, 1:NY, 1:N_procs))

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = Wa_procs(i, j, k) + dt*fzint_procs(i, j, k) - dt*Bz_procs(i, j, k)/2
                enddo
            enddo
        enddo
        call fft_solve(Q_procs, LHS_procs, Wp_procs(1:NX, 1:NY, 1:N_procs))

        ! Phiを求めるときに0番目の値も必要
        call MPI_Boundary(Up_procs)
        call MPI_Boundary(Vp_procs)
        call MPI_Boundary(Wp_procs)
        call PBM(Up_procs)
        call PBM(Vp_procs)
        call PBM(Wp_procs)
    end subroutine ibm_predict


    subroutine ave_init(U_ave_procs, V_ave_procs, W_ave_procs, U_rms_procs, V_rms_procs, W_rms_procs)
        real(8), allocatable :: U_ave_procs(:, :, :), V_ave_procs(:, :, :), W_ave_procs(:, :, :)
        real(8), allocatable :: U_rms_procs(:, :, :), V_rms_procs(:, :, :), W_rms_procs(:, :, :)
        allocate(U_ave_procs(1:NX, 1:NY, 1:N_procs), V_ave_procs(1:NX, 1:NY, 1:N_procs), W_ave_procs(1:NX, 1:NY, 1:N_procs))
        allocate(U_rms_procs(1:NX, 1:NY, 1:N_procs), V_rms_procs(1:NX, 1:NY, 1:N_procs), W_rms_procs(1:NX, 1:NY, 1:N_procs))
        U_ave_procs(:, :, :) = 0.0d0
        V_ave_procs(:, :, :) = 0.0d0
        W_ave_procs(:, :, :) = 0.0d0
        U_rms_procs(:, :, :) = 0.0d0
        V_rms_procs(:, :, :) = 0.0d0
        W_rms_procs(:, :, :) = 0.0d0
    end subroutine ave_init


    subroutine ibm_write(U_rms_procs, V_rms_procs, W_rms_procs)
        real(8) U_rms_procs(1:NX, 1:NY, 1:N_procs), V_rms_procs(1:NX, 1:NY, 1:N_procs), W_rms_procs(1:NX, 1:NY, 1:N_procs)
        real(8) U_var, V_var, W_var, U_var_sum, V_var_sum, W_var_sum
        real(8) var

        U_rms_procs(:, :, :) = sqrt(U_rms_procs(:, :, :)/Nstep)
        V_rms_procs(:, :, :) = sqrt(V_rms_procs(:, :, :)/Nstep)
        W_rms_procs(:, :, :) = sqrt(W_rms_procs(:, :, :)/Nstep)

        U_var = sum(U_rms_procs(:, :, :))/(NX*NY*NZ) * U_C
        V_var = sum(V_rms_procs(:, :, :))/(NX*NY*NZ) * U_C
        W_var = sum(W_rms_procs(:, :, :))/(NX*NY*NZ) * U_C

        call MPI_Allreduce(U_var, U_var_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_var, V_var_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_var, W_var_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        var = sqrt((U_var_sum**2 + V_var_sum**2 + W_var_sum**2)/3.0d0)

        if (myrank == 0) write(*, *) U_var_sum, V_var_sum, W_var_sum, var

    end subroutine ibm_write
end module ibm

program main
    use smac
    use fft
    use ibm
    implicit none
    real(8), allocatable :: U_procs(:, :, :), V_procs(:, :, :), W_procs(:, :, :)
    real(8), allocatable :: P_procs(:, :, :), Phi_procs(:, :, :)
    real(8), allocatable :: Ax_procs(:, :, :), Ay_procs(:, :, :), Az_procs(:, :, :)
    real(8), allocatable :: Bx_procs(:, :, :), By_procs(:, :, :), Bz_procs(:, :, :)
    real(8), allocatable :: Ax0_procs(:, :, :), Ay0_procs(:, :, :), Az0_procs(:, :, :)
    real(8), allocatable :: Bx0_procs(:, :, :), By0_procs(:, :, :), Bz0_procs(:, :, :)
    real(8), allocatable :: Up_procs(:, :, :), Vp_procs(:, :, :), Wp_procs(:, :, :)
    real(8), allocatable :: Fx_procs(:, :, :), Fy_procs(:, :, :), Fz_procs(:, :, :)
    ! ibmのために追加した変数
    real(8), allocatable :: X_procs(:, :, :), Y_procs(:, :, :), Z_procs(:, :, :)
    real(8), allocatable :: Xc_procs(:, :, :), Yc_procs(:, :, :), Zc_procs(:, :, :)
    real(8), allocatable :: Uc_procs(:, :, :), Vc_procs(:, :, :), Wc_procs(:, :, :)
    real(8), allocatable :: Ua_procs(:, :, :), Va_procs(:, :, :), Wa_procs(:, :, :)
    real(8), allocatable :: Fxc_procs(:, :, :), Fyc_procs(:, :, :), Fzc_procs(:, :, :)
    real(8), allocatable :: fxint_procs(:, :, :), fyint_procs(:, :, :), fzint_procs(:, :, :)
    ! 実験と定量的な評価をするために追加した変数
    real(8), allocatable :: U_ave_procs(:, :, :), V_ave_procs(:, :, :), W_ave_procs(:, :, :)
    real(8), allocatable :: U_rms_procs(:, :, :), V_rms_procs(:, :, :), W_rms_procs(:, :, :)
    real(8) Energy(0:NX)  ! エネルギーカスケードのために保存
    integer s
    integer step
    real(8) time1, time2
    Energy(:) = 0.0d0
    call mk_dir(dir)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

    ! write(*, *) 'bbb'
    ! call MPI_Barrier(MPI_COMM_WORLD)  ! これより前に処理差がつく作業(例えばwrite)をしないとエラーがでる。
    time1 = MPI_Wtime()   ! 計測開始

    call init(U_procs, V_procs, W_procs, P_procs, Phi_procs, &
              Ax_procs, Ay_procs, Az_procs, Bx_procs, By_procs, Bz_procs, &
              Ax0_procs, Ay0_procs, Az0_procs, Bx0_procs, By0_procs, Bz0_procs, &
              Up_procs, Vp_procs, Wp_procs, Fx_procs, Fy_procs, Fz_procs)
    if (input_step > 0) call input(U_procs, V_procs, W_procs, P_procs)
    if (method == 2 .or. method == 3) then
        call ibm_init(X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, Uc_procs, Vc_procs, Wc_procs, &
                      Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
    endif
    if (method == 3) call ave_init(U_ave_procs, V_ave_procs, W_ave_procs, U_rms_procs, V_rms_procs, W_rms_procs)


    do step = 1, Nstep
        call convection(U_procs, V_procs, W_procs, Ax_procs, Ay_procs, Az_procs)
        call viscous(U_procs, V_procs, W_procs, Bx_procs, By_procs, Bz_procs)

        if (method == 0) then
            call navier(U_procs, V_procs, W_procs, P_procs, Up_procs, Vp_procs, Wp_procs, Ax_procs, Ay_procs, Az_procs, &
                        Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, Bx0_procs, By0_procs, Bz0_procs, &
                        Fx_procs, Fy_procs, Fz_procs, step)
            call poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs, step)
            call march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
        endif
        
        if (method == 1) then
            call fft_navier(U_procs, V_procs, W_procs, P_procs, Up_procs, Vp_procs, Wp_procs, Ax_procs, Ay_procs, Az_procs, &
                            Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, Fx_procs, Fy_procs, Fz_procs, step)
            call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
            call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
        endif

        if (method == 2) then
            call ibm_preliminary(U_procs, V_procs, W_procs, P_procs, Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                                 Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, Up_procs, Vp_procs, Wp_procs, step)
            do s = 1, NS
                call ibm_Helmholtz(Up_procs, Vp_procs, Wp_procs, X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, &
                                   Uc_procs, Vc_procs, Wc_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
                call ibm_predict(Ua_procs, Va_procs, Wa_procs, fxint_procs, fyint_procs, fzint_procs, Bx_procs, By_procs, Bz_procs, Up_procs, Vp_procs, Wp_procs)
            enddo
            call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
            call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
        endif

        if (method == 3) then
            call ibm_preliminary(U_procs, V_procs, W_procs, P_procs, Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                                 Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, Up_procs, Vp_procs, Wp_procs, step)
            call ibm_Helmholtz(Up_procs, Vp_procs, Wp_procs, X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, &
                            Uc_procs, Vc_procs, Wc_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
            call ibm_predict(Ua_procs, Va_procs, Wa_procs, fxint_procs, fyint_procs, fzint_procs, Bx_procs, By_procs, Bz_procs, Up_procs, Vp_procs, Wp_procs)
            call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
            call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)

            U_ave_procs(:, :, :) = U_ave_procs(:, :, :) + U_procs(1:NX, 1:NY, 1:N_procs)
            V_ave_procs(:, :, :) = V_ave_procs(:, :, :) + V_procs(1:NX, 1:NY, 1:N_procs)
            W_ave_procs(:, :, :) = W_ave_procs(:, :, :) + W_procs(1:NX, 1:NY, 1:N_procs)
        endif

        if (mod(step, Gstep)==0) call get_data(U_procs, V_procs, W_procs, step)
        ! if (mod(step, 10)==0) call taylor_debag(U_procs, V_procs, W_procs, step)
        if (mod(step, output_step)==0) call output(U_procs, V_procs, W_procs, P_procs, step)
        call logging(U_procs, V_procs, W_procs, step, time1)


        call energy_sum(U_procs, V_procs, W_procs, Energy)  ! エネルギーの足し算
        if (mod(step, Estep)==0) call energy_reset(U_procs, V_procs, W_procs, step, Energy)
    enddo


    if (method == 3) then
        call input(U_procs, V_procs, W_procs, P_procs)
        U_ave_procs(:, :, :) = U_ave_procs(:, :, :)/Nstep
        V_ave_procs(:, :, :) = V_ave_procs(:, :, :)/Nstep
        W_ave_procs(:, :, :) = W_ave_procs(:, :, :)/Nstep

        do step = 1, Nstep
            call convection(U_procs, V_procs, W_procs, Ax_procs, Ay_procs, Az_procs)
            call viscous(U_procs, V_procs, W_procs, Bx_procs, By_procs, Bz_procs)
            call ibm_preliminary(U_procs, V_procs, W_procs, P_procs, Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                                 Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, Up_procs, Vp_procs, Wp_procs, step)
            call ibm_Helmholtz(Up_procs, Vp_procs, Wp_procs, X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, &
                                Uc_procs, Vc_procs, Wc_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
            call ibm_predict(Ua_procs, Va_procs, Wa_procs, fxint_procs, fyint_procs, fzint_procs, Bx_procs, By_procs, Bz_procs, Up_procs, Vp_procs, Wp_procs)
            call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
            call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)

            U_rms_procs(:, :, :) = U_rms_procs(:, :, :) + (U_procs(1:NX, 1:NY, 1:N_procs) - U_ave_procs(:, :, :))**2
            V_rms_procs(:, :, :) = V_rms_procs(:, :, :) + (V_procs(1:NX, 1:NY, 1:N_procs) - V_ave_procs(:, :, :))**2
            W_rms_procs(:, :, :) = W_rms_procs(:, :, :) + (W_procs(1:NX, 1:NY, 1:N_procs) - W_ave_procs(:, :, :))**2

            call logging(U_procs, V_procs, W_procs, step, time1)
        enddo

        call ibm_write(U_rms_procs, V_rms_procs, W_rms_procs)
    endif

    ! call MPI_Barrier(MPI_COMM_WORLD)
    time2 = MPI_Wtime()  ! 計測終了
    if (myrank == 0) write(*, *) 'time:', time2 - time1

    call MPI_Finalize(ierr)  ! 終わり
    
end program main