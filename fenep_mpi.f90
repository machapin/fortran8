module smac
    implicit none
    include 'mpif.h'
    real(8), parameter :: PI = acos(-1.0d0)
    ! ステップ数
    integer, parameter :: Nstep = 100  ! 7桁まで
    integer, parameter :: Tstep = 10 ! 時系列データを取得する間隔
    integer, parameter :: Lstep = 10  ! ログを表示する間隔
    integer, parameter :: Gstep = 100  ! Mapやエンストロフィやエネルギーを保存する間隔
    integer, parameter :: Ostep = 5000  ! 配列を保存する間隔
    ! パラメータ
    integer, parameter :: NX = 256, NY = NX, NZ = NX
    real(8), parameter :: dX = 2*PI/NX, dY = 2*PI/NY, dZ = 2*PI/NZ  ! 規格化長さは2*PI
    real(8), parameter :: dt = 0.002d0
    ! 手法
    integer, parameter :: input_type = 1  ! 0:初期条件なし、1:初期条件ありで1周、2:初期値ありで2周
    integer, parameter :: method = 2  ! 1:FFT、2:IBM
    integer, parameter :: ibm_type = 22  ! 11:円柱表面1つ、12:円柱表面4つ、21:円柱内部1つ、22:円柱内部4つ
    real(8), parameter :: DC = 2*PI / 4.0d0  ! 円柱の直径、円柱の中心点を確認すること！
    integer, parameter :: flow_type = 0  ! 0:外力なし(f=0)、3:テイラーグリーン外力、4:テイラーグリーン渦の減衰
    integer, parameter :: eigen_method = 3  ! 0:カルダノ、1:シルベスター、3:ビエト
    ! 無次元パラメータ
    ! real(8), parameter :: Re = 2.0d0
    real(8), parameter :: beta = 1.0d0
    real(8), parameter :: Wi = 5.0d0
    real(8), parameter :: Lp = 100.0d0
    ! 有次元パラメータ
    ! real(8), parameter :: L_C = 1.0d0  ! 長さが100なら100/2*PI
    ! real(8), parameter :: U_C = 1.0d0  ! 本来は乱流テイラーグリーン渦の平均流の速さ
    ! real(8), parameter :: nu = L_C*U_C/Re
    ! 実験と同様のパラメータ
    real(8), parameter :: f_C = 0.4d0  ! 円柱回転速度[rps]
    real(8), parameter :: D_C = 0.03d0  ! 円柱直径[m]
    real(8), parameter :: U_C = f_C * PI * D_C  ! 代表速度
    real(8), parameter :: L_C = D_C / DC  ! 代表長さ
    real(8), parameter :: nu = 1.0d-6  ! 動粘性係数
    real(8), parameter :: Re = U_C*L_C/nu  ! betaが絡んできそうだけどReを合わせたいので一旦これで
    ! その他のパラメータ
    real(8), parameter :: dX_C = dX*L_C, dY_C = dY*L_C, dZ_C = dZ*L_C
    real(8), parameter :: dt_C = dt*L_C/U_C
    ! その他
    character(64) :: input_dir = './'
    character(64) :: output_dir = './debug/'
    ! MPI用変数
    integer N_procs, NY_procs
    integer ierr, procs, myrank
    integer next_rank, former_rank
    integer req1s, req1r, req2s, req2r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s, sta1r, sta2s, sta2r
    integer seed, seeds(2)  ! 乱数

contains
    subroutine init(U_procs, V_procs, W_procs, P_procs, Phi_procs, &
                    Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, &
                    Bx_procs, By_procs, Bz_procs, Bx0_procs, By0_procs, Bz0_procs, &
                    Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, &
                    Up_procs, Vp_procs, Wp_procs, Fx_procs, Fy_procs, Fz_procs, C_procs, &
                    Cpx_procs, Cnx_procs, Cpy_procs, Cny_procs, Cpz_procs, Cnz_procs, Cx_procs)  ! はじめに実行
        real(8), allocatable :: U_procs(:, :, :), V_procs(:, :, :), W_procs(:, :, :)
        real(8), allocatable :: P_procs(:, :, :), Phi_procs(:, :, :)
        real(8), allocatable :: Ax_procs(:, :, :), Ay_procs(:, :, :), Az_procs(:, :, :)
        real(8), allocatable :: Ax0_procs(:, :, :), Ay0_procs(:, :, :), Az0_procs(:, :, :)
        real(8), allocatable :: Bx_procs(:, :, :), By_procs(:, :, :), Bz_procs(:, :, :)
        real(8), allocatable :: Bx0_procs(:, :, :), By0_procs(:, :, :), Bz0_procs(:, :, :)
        real(8), allocatable :: Tx_procs(:, :, :), Ty_procs(:, :, :), Tz_procs(:, :, :)
        real(8), allocatable :: Tx0_procs(:, :, :), Ty0_procs(:, :, :), Tz0_procs(:, :, :)
        real(8), allocatable :: Up_procs(:, :, :), Vp_procs(:, :, :), Wp_procs(:, :, :)
        real(8), allocatable :: Fx_procs(:, :, :), Fy_procs(:, :, :), Fz_procs(:, :, :)
        real(8), allocatable :: C_procs(:, :, :, :)
        real(8), allocatable :: Cpx_procs(:, :, :, :), Cnx_procs(:, :, :, :)
        real(8), allocatable :: Cpy_procs(:, :, :, :), Cny_procs(:, :, :, :)
        real(8), allocatable :: Cpz_procs(:, :, :, :), Cnz_procs(:, :, :, :)
        real(8), allocatable :: Cx_procs(:, :, :, :)
        integer i, j, k
        character(2) str
        write(str, '(I2.2)') ibm_type

        N_procs = NZ / procs
        NY_procs = NY / procs

        ! コメント
        if (mod(NZ, procs)/=0 .or. mod(NY, procs)/=0) stop 'NZ or NY is not divisible by procs'
        ! if (1.0d0*dt/dX>1.0d0/6 .or. 1.0d0*dt/dY>1.0d0/6 .or. 1.0d0*dt/dZ>1.0d0/6) stop 'CFL condition is not met.'
        if (myrank == 0) then
            write(*, '(a, F8.3, F8.3, F8.3)') 'CFL:', 1.0d0*dt/dX*6.0d0, 1.0d0*dt/dY*6.0d0, 1.0d0*dt/dZ*6.0d0
            write(*, '(a, F8.3)') 'L_C =', L_C
            write(*, '(a, F8.3)') 'U_C =', U_C
            write(*, '(a, E12.4)') 'nu =', nu
            write(*, '(a, F10.3)') 'Re  =', Re
            write(*, '(a, E12.4)') 'dX_C =', dX_C
            write(*, '(a, E12.4)') 'dt_C =', dt_C
            write(*, '(a, F8.3)') 'step /[s]=', 1.0d0/dt_C
        endif

        ! ファイルの初期化
        if (myrank == 0) then
            call rm_dir(trim(output_dir))
            call mk_dir(trim(output_dir))
            call mk_dir(trim(output_dir)//'Log/')
            call mk_dir(trim(output_dir)//'Map/')
            call mk_dir(trim(output_dir)//'Qti/')
            call mk_dir(trim(output_dir)//'Cylinder/')
            call mk_dir(trim(output_dir)//'Restart/')
        endif

        ! 乱数用
        call system_clock(seed)  ! 現在時刻をシード値として使用する
        seed = seed + myrank  ! プロセスごとにシード値を変える
        seeds(1:2) = [seed, 0]
        call random_seed(put=seeds)  ! 入力するシード値は配列

        ! 配列のallocate
        allocate(U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(P_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Phi_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Bx0_procs(1:NX, 1:NY, 1:N_procs), By0_procs(1:NX, 1:NY, 1:N_procs), Bz0_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Tx_procs(1:NX, 1:NY, 1:N_procs), Ty_procs(1:NX, 1:NY, 1:N_procs), Tz_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Tx0_procs(1:NX, 1:NY, 1:N_procs), Ty0_procs(1:NX, 1:NY, 1:N_procs), Tz0_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Fx_procs(1:NX, 1:NY, 1:N_procs), Fy_procs(1:NX, 1:NY, 1:N_procs), Fz_procs(1:NX, 1:NY, 1:N_procs))
        allocate(C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Cpx_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cnx_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Cpy_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cny_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Cpz_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cnz_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Cx_procs(6, 1:NX, 1:NY, 1:N_procs))

        ! 初期条件
        U_procs(:, :, :) = 0.0d0
        V_procs(:, :, :) = 0.0d0
        W_procs(:, :, :) = 0.0d0
        call random_number(U_procs)
        call random_number(V_procs)
        call random_number(W_procs)
        U_procs(:, :, :) = 0.001d0 * (U_procs(:, :, :)-0.5d0)
        V_procs(:, :, :) = 0.001d0 * (V_procs(:, :, :)-0.5d0)
        W_procs(:, :, :) = 0.001d0 * (W_procs(:, :, :)-0.5d0)
        P_procs(:, :, :) = 0.0d0
        Phi_procs(:, :, :) = 0.0d0
        Fx_procs(:, :, :) = 0.0d0
        Fy_procs(:, :, :) = 0.0d0
        Fz_procs(:, :, :) = 0.0d0
        call random_number(C_procs)
        C_procs(:, :, :, :) = 0.001d0 * (C_procs(:, :, :, :)-0.5d0)
        C_procs(1, :, :, :) = C_procs(1, :, :, :) + 1.0d0
        C_procs(4, :, :, :) = C_procs(4, :, :, :) + 1.0d0
        C_procs(6, :, :, :) = C_procs(6, :, :, :) + 1.0d0
        if (flow_type == 3) then
            do k = 1, N_procs
                do j = 1, NY
                    do i = 1, NX
                        Fx_procs(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)
                        Fy_procs(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)
                        U_procs(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)/(2.0d0 * nu)
                        V_procs(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)/(2.0d0 * nu)
                        ! Fx_procs(i, j, k) = -sin(i*dX) * cos((myrank*N_procs+k-0.5d0)*dZ)
                        ! Fz_procs(i, j, k) = cos((i-0.5d0)*dX) * sin((myrank*N_procs+k)*dZ)
                        ! Fx_procs(i, j, k) = -sin(2.0d0*(j-0.5d0)*dY)
                        ! Fy_procs(i, j, k) = sin(2.0d0*(i-0.5d0)*dX)
                        ! Fx_procs(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY) * cos((myrank*N_procs + k-0.5d0)*dZ)
                        ! Fy_procs(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY) * cos((myrank*N_procs + k-0.5d0)*dZ)
                    enddo
                enddo
            enddo
            ! デバッグモード
            W_procs(:, :, :) = 0.0d0
            C_procs(:, :, :, :) = 0.0d0
            C_procs(1, :, :, :) = 1.0d0
            C_procs(4, :, :, :) = 1.0d0
            C_procs(6, :, :, :) = 1.0d0
        endif
        if (flow_type == 4) then
            do k = 1, N_procs
                do j = 1, NY
                    do i = 1, NX
                        U_procs(i, j, k) = -sin(i*dX) * cos((j-0.5d0)*dY)
                        V_procs(i, j, k) = cos((i-0.5d0)*dX) * sin(j*dY)
                    enddo
                enddo
            enddo
        endif


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
        call C_MPI_Boundary(C_procs)

        ! 境界条件
        call PBM(U_procs)
        call PBM(V_procs)
        call PBM(W_procs)
        call PBM(P_procs)
        call PBM(Phi_procs)
        call C_PBM(C_procs)
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

    subroutine C_MPI_Boundary(C_procs)
        real(8), intent(inout) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        ! 次のランクへ送信
        call MPI_Isend(C_procs(1, 0, 0, N_procs), 6*(NX+2)*(NY+2), MPI_REAL8, next_rank, 1, MPI_COMM_WORLD, req1s, ierr)
        call MPI_Irecv(C_procs(1, 0, 0, 0), 6*(NX+2)*(NY+2), MPI_REAL8, former_rank, 1, MPI_COMM_WORLD, req1r, ierr)
  
        ! 前のランクへ送信
        call MPI_Isend(C_procs(1, 0, 0, 1), 6*(NX+2)*(NY+2), MPI_REAL8, former_rank, 2, MPI_COMM_WORLD, req2s, ierr)
        call MPI_Irecv(C_procs(1, 0, 0, N_procs+1), 6*(NX+2)*(NY+2), MPI_REAL8, next_rank, 2, MPI_COMM_WORLD, req2r, ierr)
  
        ! 待機
        call MPI_Wait(req1s, sta1s, ierr)
        call MPI_Wait(req1r, sta1r, ierr)
        call MPI_Wait(req2s, sta2s, ierr)
        call MPI_Wait(req2r, sta2r, ierr)
    end subroutine C_MPI_Boundary


    subroutine PBM(A_procs)  ! 一般的な周期境界条件
        real(8), intent(inout) :: A_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        A_procs(0, :, :) = A_procs(NX, :, :)
        A_procs(NX+1, :, :) = A_procs(1, :, :)
        A_procs(:, 0, :) = A_procs(:, NY, :)
        A_procs(:, NY+1, :) = A_procs(:, 1, :)
    end subroutine PBM
    
    subroutine C_PBM(C_procs)
        real(8), intent(inout) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        C_procs(:, 0, :, :) = C_procs(:, NX, :, :)
        C_procs(:, NX+1, :, :) = C_procs(:, 1, :, :)
        C_procs(:, :, 0, :) = C_procs(:, :, NY, :)
        C_procs(:, :, NY+1, :) = C_procs(:, :, 1, :)
    end subroutine C_PBM

    subroutine CpxCnx(C_procs, Cpx_procs, Cnx_procs)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Cpx_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cnx_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k, l, index
        real(8) Cptemp(6, 3), Cntemp(6, 3), Eigen(6), minlist(4)

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Cptemp(:, 1) = C_procs(:, i, j, k)*3/2 - C_procs(:, i+1, j, k)/2
                    Cntemp(:, 1) = C_procs(:, i, j, k)/2 + C_procs(:, i+1, j, k)/2
                    Cptemp(:, 2) = C_procs(:, i-1, j, k)/2 + C_procs(:, i, j, k)/2
                    Cntemp(:, 2) = -C_procs(:, i-1, j, k)/2 + C_procs(:, i, j, k)*3/2
                    Cptemp(:, 3) = C_procs(:, i-1, j, k)/4 + C_procs(:, i, j, k) - C_procs(:, i+1, j, k)/4
                    Cntemp(:, 3) = -C_procs(:, i-1, j, k)/4 + C_procs(:, i, j, k) + C_procs(:, i+1, j, k)/4

                    do l = 1, 3
                        if (eigen_method == 0) then
                            call Cardano(Cptemp(:, l), Eigen(1), Eigen(2), Eigen(3))  ! CardanoでCpnの95%ぐらい時間かかってる。
                            call Cardano(Cntemp(:, l), Eigen(4), Eigen(5), Eigen(6))
                            minlist(l) = minval(Eigen(:))
                        endif
                        if (eigen_method == 1) then
                            call Sylvester(Cptemp(:, l), Eigen(1))
                            call Sylvester(Cntemp(:, l), Eigen(2))
                            minlist(l) = minval(Eigen(1:2))
                        endif
                        if (eigen_method == 3) then
                            call Vieta(Cptemp(:, l), Eigen(1))
                            call Vieta(Cntemp(:, l), Eigen(2))
                            minlist(l) = minval(Eigen(1:2))
                        endif
                    enddo
                    minlist(4) = 0.0d0
                    index = maxloc(minlist, 1)

                    if (index == 4) then
                        Cpx_procs(:, i, j, k) = C_procs(:, i, j, k)
                        Cnx_procs(:, i, j, k) = C_procs(:, i, j, k)
                    else
                        if (eigen_method == 1) index = 3
                        Cpx_procs(:, i, j, k) = Cptemp(:, index)  ! 周期条件をそのまま使いたいため保存する場所を変更
                        Cnx_procs(:, i, j, k) = Cntemp(:, index)
                    endif
                enddo
            enddo
        enddo

        call C_MPI_Boundary(Cpx_procs)
        call C_MPI_Boundary(Cnx_procs)
        call C_PBM(Cpx_procs)
        call C_PBM(Cnx_procs)
    end subroutine CpxCnx

    subroutine CpyCny(C_procs, Cpy_procs, Cny_procs)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Cpy_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cny_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k, l, index
        real(8) Cptemp(6, 3), Cntemp(6, 3), Eigen(6), minlist(4)

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Cptemp(:, 1) = C_procs(:, i, j, k)*3/2 - C_procs(:, i, j+1, k)/2
                    Cntemp(:, 1) = C_procs(:, i, j, k)/2 + C_procs(:, i, j+1, k)/2
                    Cptemp(:, 2) = C_procs(:, i, j-1, k)/2 + C_procs(:, i, j, k)/2
                    Cntemp(:, 2) = -C_procs(:, i, j-1, k)/2 + C_procs(:, i, j, k)*3/2
                    Cptemp(:, 3) = C_procs(:, i, j-1, k)/4 + C_procs(:, i, j, k) - C_procs(:, i, j+1, k)/4
                    Cntemp(:, 3) = -C_procs(:, i, j-1, k)/4 + C_procs(:, i, j, k) + C_procs(:, i, j+1, k)/4

                    minlist(0:3) = 0.0d0
                    do l = 1, 3
                        if (eigen_method == 0) then
                            call Cardano(Cptemp(:, l), Eigen(1), Eigen(2), Eigen(3))  ! CardanoでCpnの95%ぐらい時間かかってる。
                            call Cardano(Cntemp(:, l), Eigen(4), Eigen(5), Eigen(6))
                            minlist(l) = minval(Eigen(:))
                        endif
                        if (eigen_method == 1) then
                            call Sylvester(Cptemp(:, l), Eigen(1))
                            call Sylvester(Cntemp(:, l), Eigen(2))
                            minlist(l) = minval(Eigen(1:2))
                        endif
                        if (eigen_method == 3) then
                            call Vieta(Cptemp(:, l), Eigen(1))
                            call Vieta(Cntemp(:, l), Eigen(2))
                            minlist(l) = minval(Eigen(1:2))
                        endif
                    enddo
                    minlist(4) = 0.0d0
                    index = maxloc(minlist, 1)

                    if (index == 4) then
                        Cpy_procs(:, i, j, k) = C_procs(:, i, j, k)
                        Cny_procs(:, i, j, k) = C_procs(:, i, j, k)
                    else
                        if (eigen_method == 1) index = 3
                        Cpy_procs(:, i, j, k) = Cptemp(:, index)
                        Cny_procs(:, i, j, k) = Cntemp(:, index)
                    endif
                enddo
            enddo
        enddo

        call C_MPI_Boundary(Cpy_procs)
        call C_MPI_Boundary(Cny_procs)
        call C_PBM(Cpy_procs)
        call C_PBM(Cny_procs)
    end subroutine CpyCny

    subroutine CpzCnz(C_procs, Cpz_procs, Cnz_procs)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Cpz_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cnz_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k, l, index
        real(8) Cptemp(6, 3), Cntemp(6, 3), Eigen(6), minlist(4)

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Cptemp(:, 1) = C_procs(:, i, j, k)*3/2 - C_procs(:, i, j, k+1)/2
                    Cntemp(:, 1) = C_procs(:, i, j, k)/2 + C_procs(:, i, j, k+1)/2
                    Cptemp(:, 2) = C_procs(:, i, j, k-1)/2 + C_procs(:, i, j, k)/2
                    Cntemp(:, 2) = -C_procs(:, i, j, k-1)/2 + C_procs(:, i, j, k)*3/2
                    Cptemp(:, 3) = C_procs(:, i, j, k-1)/4 + C_procs(:, i, j, k) - C_procs(:, i, j, k+1)/4
                    Cntemp(:, 3) = -C_procs(:, i, j, k-1)/4 + C_procs(:, i, j, k) + C_procs(:, i, j, k+1)/4

                    minlist(0:3) = 0.0d0
                    do l = 1, 3
                        if (eigen_method == 0) then
                            call Cardano(Cptemp(:, l), Eigen(1), Eigen(2), Eigen(3))  ! CardanoでCpnの95%ぐらい時間かかってる。
                            call Cardano(Cntemp(:, l), Eigen(4), Eigen(5), Eigen(6))
                            minlist(l) = minval(Eigen(:))
                        endif
                        if (eigen_method == 1) then
                            call Sylvester(Cptemp(:, l), Eigen(1))
                            call Sylvester(Cntemp(:, l), Eigen(2))
                            minlist(l) = minval(Eigen(1:2))
                        endif
                        if (eigen_method == 3) then
                            call Vieta(Cptemp(:, l), Eigen(1))
                            call Vieta(Cntemp(:, l), Eigen(2))
                            minlist(l) = minval(Eigen(1:2))
                        endif
                    enddo
                    minlist(4) = 0.0d0
                    index = maxloc(minlist, 1)

                    if (index == 4) then
                        Cpz_procs(:, i, j, k) = C_procs(:, i, j, k)
                        Cnz_procs(:, i, j, k) = C_procs(:, i, j, k)
                    else
                        if (eigen_method == 1) index = 3
                        Cpz_procs(:, i, j, k) = Cptemp(:, index)
                        Cnz_procs(:, i, j, k) = Cntemp(:, index)
                    endif
                enddo
            enddo
        enddo

        call C_MPI_Boundary(Cpz_procs)
        call C_MPI_Boundary(Cnz_procs)
        call C_PBM(Cpz_procs)
        call C_PBM(Cnz_procs)
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

            w0 = cmplx(1.0d0, 0.0d0, kind=8)
            w1 = cmplx(-0.5d0, sqrt(3.0d0)/2, kind=8)
            w2 = cmplx(-0.5d0, -sqrt(3.0d0)/2, kind=8)

            e0 = w0*u3**(1.0d0/3.0d0) + w0*v3**(1.0d0/3.0d0) - a2/3
            e1 = w1*u3**(1.0d0/3.0d0) + w2*v3**(1.0d0/3.0d0) - a2/3
            e2 = w2*u3**(1.0d0/3.0d0) + w1*v3**(1.0d0/3.0d0) - a2/3

            re0 = real(e0)
            re1 = real(e1)
            re2 = real(e2)
        endif
    end subroutine Cardano

    subroutine Vieta(a, x)
        real(8), intent(in) :: a(6)
        real(8), intent(out) :: x
        real(8) a0, a1, a2, p, q, t

        a2 = -(a(1) + a(4) + a(6))
        a1 = a(1)*a(4) + a(4)*a(6) + a(6)*a(1) - a(5)*a(5) - a(3)*a(3) - a(2)*a(2)
        a0 = -(a(1)*a(4)*a(6) + a(2)*a(5)*a(3) + a(2)*a(5)*a(3) - a(1)*a(5)*a(5) - a(2)*a(2)*a(6) - a(3)*a(4)*a(3))

        p = a1 - a2**2/3.0d0
        q = a0 - a1*a2/3.0d0 + 2.0d0*a2**3/27.0d0
        t = (q/2.0d0)**2 + (p/3.0d0)**3  ! 三次方程式の判別式

        if (t >= 0.0d0) then
            x = -100.0d0  ! 選ばれないようにする
        else
            x = 2.0d0*sqrt(-p/3.0d0)*cos(acos(3.0d0*q/(2.0d0*p)*sqrt(-3.0d0/p))/3.0d0 + 2.0d0*PI/3.0d0) - a2/3.0d0
        endif

    end subroutine Vieta

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

    subroutine Cstar(Cpx_procs, Cnx_procs, Cpy_procs, Cny_procs, Cpz_procs, Cnz_procs, U_procs, V_procs, W_procs, Cx_procs)
        real(8), intent(in) :: Cpx_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cnx_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Cpy_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cny_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Cpz_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), Cnz_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Cx_procs(6, 1:NX, 1:NY, 1:N_procs)
        integer i, j, k
        real(8) s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12
        real(8) u1, u2, v1, v2, w1, w2

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    u1 = U_procs(i-1, j, k)
                    u2 = U_procs(i, j, k)
                    v1 = V_procs(i, j-1, k)
                    v2 = V_procs(i, j, k)
                    w1 = W_procs(i, j, k-1)
                    w2 = W_procs(i, j, k)

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

                    Cx_procs(:, i, j, k) = s3*Cpx_procs(:, i, j, k) + s1*Cpx_procs(:, i+1, j, k) + s4*Cnx_procs(:, i-1, j, k) + s2*Cnx_procs(:, i, j, k) &
                                          +s7*Cpy_procs(:, i, j, k) + s5*Cpy_procs(:, i, j+1, k) + s8*Cny_procs(:, i, j-1, k) + s6*Cny_procs(:, i, j, k) &
                                          +s11*Cpz_procs(:, i, j, k) + s9*Cpz_procs(:, i, j, k+1) + s12*Cnz_procs(:, i, j, k-1) + s10*Cnz_procs(:, i, j, k)
                enddo
            enddo
        enddo

    end subroutine Cstar

    function f(C) result(y)  ! FENE-Pモデルの核
        real(8), intent(in) :: C(6) 
        real(8) y  ! y=f(C)
        y = (Lp**2 - 3)/(Lp**2 - C(1)-C(4)-C(6))
        ! y = 1.0d0  ! Oldroyd-Bモデル
    end function f

    subroutine Lyapunov(Cx_procs, U_procs, V_procs, W_procs, C_procs)
        real(8), intent(in) :: Cx_procs(6, 1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k
        integer row, colum, itr
        real(8) A0(3, 3), b(6)
        real(8) x(6), x1(6), x2(6), l(6), l1(6), l2(6), Jac(6, 6), r(6)
        real(8), parameter :: dc = 1.0d-5  ! 3, 4, 5あたりがいい。
        integer, parameter :: itrmax = 100
        real(8), parameter :: eps = 1.0d-5
        ! newton_itr = 0

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    ! 定数部分を計算
                    b(:) = Cx_procs(:, i, j, k)
                    b(1) = b(1) + dt/Wi
                    b(4) = b(4) + dt/Wi
                    b(6) = b(6) + dt/Wi

                    A0(1, 1) = -dt*(-U_procs(i-1, j, k) + U_procs(i, j, k))/dX
                    A0(1, 2) = -dt*(-V_procs(i-1, j-1, k) + V_procs(i+1, j-1, k) - V_procs(i-1, j, k) + V_procs(i+1, j, k))/(4*dX)
                    A0(1, 3) = -dt*(-W_procs(i-1, j, k-1) + W_procs(i+1, j, k-1) - W_procs(i-1, j, k) + W_procs(i+1, j, k))/(4*dX)
                    A0(2, 1) = -dt*(-U_procs(i-1, j-1, k) - U_procs(i, j-1, k) + U_procs(i-1, j+1, k) + U_procs(i, j+1, k))/(4*dY)
                    A0(2, 2) = -dt*(-V_procs(i, j-1, k) + V_procs(i, j, k))/dY
                    A0(2, 3) = -dt*(-W_procs(i, j-1, k-1) + W_procs(i, j+1, k-1) - W_procs(i, j-1, k) + W_procs(i, j+1, k))/(4*dY)
                    A0(3, 1) = -dt*(-U_procs(i-1, j, k-1) - U_procs(i, j, k-1) + U_procs(i-1, j, k+1) + U_procs(i, j, k+1))/(4*dZ)
                    A0(3, 2) = -dt*(-V_procs(i, j-1, k-1) - V_procs(i, j, k-1) + V_procs(i, j-1, k+1) + V_procs(i, j, k+1))/(4*dZ)
                    A0(3, 3) = -dt*(-W_procs(i, j, k-1) + W_procs(i, j, k))/dZ
                    A0 = transpose(A0)  ! 転置して関数に渡す

                    ! ニュートン法
                    x(:) = C_procs(:, i, j, k)  ! Cから探索開始
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
                        
                        call gauss(Jac, l, r)  ! ガウスの消去法
                        x(:) = x(:) - r(:)

                        if (sqrt(sum(r**2)/sum(x**2)) < eps) exit
                    enddo
                    ! newton_itr = newton_itr + itr
                    C_procs(:, i, j, k) = x(:)
                enddo
            enddo
        enddo

        call C_MPI_Boundary(C_procs)
        call C_PBM(C_procs)
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

    subroutine polymer_stress(C_procs, Tx_procs, Ty_procs, Tz_procs)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Tx_procs(1:NX, 1:NY, 1:N_procs), Ty_procs(1:NX, 1:NY, 1:N_procs), Tz_procs(1:NX, 1:NY, 1:N_procs)
        integer i, j, k
        real(8) C0_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), E(6)

        ! C0_procs(:, :, :, :) = 0.0d0
        E(:) = (/1.0, 0.0, 0.0, 1.0, 0.0, 1.0/)
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    C0_procs(:, i, j, k) = (f(C_procs(:, i, j, k))*C_procs(:, i, j, k) - E(:))
                enddo
            enddo
        enddo
        call C_MPI_Boundary(C0_procs)
        call C_PBM(C0_procs)

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Tx_procs(i, j, k) = (-C0_procs(1, i, j, k) + C0_procs(1, i+1, j, k))/dX &
                                       +(-C0_procs(2, i, j-1, k) - C0_procs(2, i+1, j-1, k) + C0_procs(2, i, j+1, k) + C0_procs(2, i+1, j+1, k))/(4*dY) &
                                       +(-C0_procs(3, i, j, k-1) - C0_procs(3, i+1, j, k-1) + C0_procs(3, i, j, k+1) + C0_procs(3, i+1, j, k+1))/(4*dZ)
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Ty_procs(i, j, k) = (-C0_procs(2, i-1, j, k) + C0_procs(2, i+1, j, k) - C0_procs(2, i-1, j+1, k) + C0_procs(2, i+1, j+1, k))/(4*dX) &
                                       +(-C0_procs(4, i, j, k) + C0_procs(4, i, j+1, k))/dY &
                                       +(-C0_procs(5, i, j, k-1) - C0_procs(5, i, j+1, k-1) + C0_procs(5, i, j, k+1) + C0_procs(5, i, j+1, k+1))/(4*dZ)
                enddo
            enddo
        enddo
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Tz_procs(i, j, k) = (-C0_procs(3, i-1, j, k) + C0_procs(3, i+1, j, k) - C0_procs(3, i-1, j, k+1) + C0_procs(3, i+1, j, k+1))/(4*dX) &
                                       +(-C0_procs(5, i, j-1, k) + C0_procs(5, i, j+1, k) - C0_procs(5, i, j-1, k+1) + C0_procs(5, i, j+1, k+1))/(4*dY) &
                                       +(-C0_procs(6, i, j, k) + C0_procs(6, i, j, k+1))/dZ
                enddo
            enddo
        enddo
        Tx_procs(:, :, :) = (1.0d0-beta)/(Re*Wi)*Tx_procs(:, :, :)
        Ty_procs(:, :, :) = (1.0d0-beta)/(Re*Wi)*Ty_procs(:, :, :)
        Tz_procs(:, :, :) = (1.0d0-beta)/(Re*Wi)*Tz_procs(:, :, :)
    end subroutine polymer_stress

    
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
        Bx_procs(:, :, :) = beta/Re*Bx_procs(:, :, :)
        By_procs(:, :, :) = beta/Re*By_procs(:, :, :)
        Bz_procs(:, :, :) = beta/Re*Bz_procs(:, :, :)

    end subroutine viscous


    subroutine get_data_binary(U_procs, V_procs, W_procs, C_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        real(8) U(0:NX+1, 0:NY+1, 0:NZ+1), V(0:NX+1, 0:NY+1, 0:NZ+1), W(0:NX+1, 0:NY+1, 0:NZ+1)
        real(8) C(6, 0:NX+1, 0:NY+1, 0:NZ+1)
        integer i, j, k
        real(8) D(3, 3), Omega(3), S(3), Qti, trC, E(3)
        character(8) str
        
        call MPI_Gather(U_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, U(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(V_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, V(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(W_procs(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, W(0, 0, 1), (NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(C_procs(1, 0, 0, 1), 6*(NX+2)*(NY+2)*N_procs, MPI_REAL8, C(1, 0, 0, 1), 6*(NX+2)*(NY+2)*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        if (myrank == 0) then
            U(:, :, 0) = U(:, :, NZ)
            U(:, :, NZ+1) = U(:, :, 1)
            V(:, :, 0) = V(:, :, NZ)
            V(:, :, NZ+1) = V(:, :, 1)
            W(:, :, 0) = W(:, :, NZ)
            W(:, :, NZ+1) = W(:, :, 1)
            C(:, :, :, 0) = C(:, :, :, NZ)
            C(:, :, :, NZ+1) = C(:, :, :, 1)

            write(str, '(I8.8)') step
            open(10, file=trim(output_dir)//'Map/z_'//str//'.bin', form='unformatted', status='replace', access='stream')
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
                    D(:, :) = D(:, :)*U_C/L_C

                    E(1) = (U(i-1, j, k)+U(i, j, k))/2*U_C
                    E(2) = (V(i, j-1, k)+V(i, j, k))/2*U_C
                    E(3) = (W(i, j, k-1)+W(i, j, k))/2*U_C
                    Omega(1) = (D(3, 2) - D(2, 3))
                    Omega(2) = (D(1, 3) - D(3, 1))
                    Omega(3) = (D(2, 1) - D(1, 2))
                    S(1) = (D(3, 2) + D(2, 3))
                    S(2) = (D(1, 3) + D(3, 1))
                    S(3) = (D(2, 1) + D(1, 2))
                    Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                    trC = C(1, i, j, k) + C(4, i, j, k) + C(6, i, j, k)
                    write(10) (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, E(1), E(2), E(3), Qti, &
                              Omega(1), Omega(2), Omega(3), sum(Omega**2)/2, &
                              trC, C(1, i, j, k), C(2, i, j, k), C(3, i, j, k), C(4, i, j, k), C(5, i, j, k), C(6, i, j, k)
                enddo
            enddo
            close(10)


            open(20, file=trim(output_dir)//'Map/y_'//str//'.bin', form='unformatted', status='replace', access='stream')
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
                    D(:, :) = D(:, :)*U_C/L_C
    
                    E(1) = (U(i-1, j, k)+U(i, j, k))/2*U_C
                    E(2) = (V(i, j-1, k)+V(i, j, k))/2*U_C
                    E(3) = (W(i, j, k-1)+W(i, j, k))/2*U_C
                    Omega(1) = (D(3, 2) - D(2, 3))
                    Omega(2) = (D(1, 3) - D(3, 1))
                    Omega(3) = (D(2, 1) - D(1, 2))
                    S(1) = (D(3, 2) + D(2, 3))
                    S(2) = (D(1, 3) + D(3, 1))
                    S(3) = (D(2, 1) + D(1, 2))
                    Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                    trC = C(1, i, j, k) + C(4, i, j, k) + C(6, i, j, k)
                    write(20) (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, E(1), E(2), E(3), Qti, &
                              Omega(1), Omega(2), Omega(3), sum(Omega**2)/2, &
                              trC, C(1, i, j, k), C(2, i, j, k), C(3, i, j, k), C(4, i, j, k), C(5, i, j, k), C(6, i, j, k)
                enddo
            enddo
            close(20)

            open(30, file=trim(output_dir)//'Map/x_'//str//'.bin', form='unformatted', status='replace', access='stream')
            i = NX/2  ! x-z平面
            do k = 1, NZ
                do j = 1, NY
                    D(1, 1) = (U(i, j, k) - U(i-1, j, k))/dX
                    D(1, 2) = (U(i, j+1, k) - U(i, j-1, k) + U(i-1, j+1, k) - U(i-1, j-1, k))/(4*dY)
                    D(1, 3) = (U(i, j, k+1) - U(i, j, k-1) + U(i-1, j, k+1) - U(i-1, j, k-1))/(4*dZ)
                    D(2, 1) = (V(i+1, j, k) - V(i-1, j, k) + V(i+1, j-1, k) - V(i-1, j-1, k))/(4*dX)
                    D(2, 2) = (V(i, j, k) - V(i, j-1, k))/dY
                    D(2, 3) = (V(i, j, k+1) - V(i, j, k-1) + V(i, j-1, k+1) - V(i, j-1, k-1))/(4*dZ)
                    D(3, 1) = (W(i+1, j, k) - W(i-1, j, k) + W(i+1, j, k-1) - W(i-1, j, k-1))/(4*dX)
                    D(3, 2) = (W(i, j+1, k) - W(i, j-1, k) + W(i, j+1, k-1) - W(i, j-1, k-1))/(4*dY)
                    D(3, 3) = (W(i, j, k) - W(i, j, k-1))/dZ
                    D(:, :) = D(:, :)*U_C/L_C
    
                    E(1) = (U(i-1, j, k)+U(i, j, k))/2*U_C
                    E(2) = (V(i, j-1, k)+V(i, j, k))/2*U_C
                    E(3) = (W(i, j, k-1)+W(i, j, k))/2*U_C
                    Omega(1) = (D(3, 2) - D(2, 3))
                    Omega(2) = (D(1, 3) - D(3, 1))
                    Omega(3) = (D(2, 1) - D(1, 2))
                    S(1) = (D(3, 2) + D(2, 3))
                    S(2) = (D(1, 3) + D(3, 1))
                    S(3) = (D(2, 1) + D(1, 2))
                    Qti = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)
                    trC = C(1, i, j, k) + C(4, i, j, k) + C(6, i, j, k)
                    write(30) (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, (k-0.5d0)*dZ_C, E(1), E(2), E(3), Qti, &
                              Omega(1), Omega(2), Omega(3), sum(Omega**2)/2, &
                              trC, C(1, i, j, k), C(2, i, j, k), C(3, i, j, k), C(4, i, j, k), C(5, i, j, k), C(6, i, j, k)
                enddo
            enddo
            close(30)
        endif
    end subroutine get_data_binary

    subroutine input_binary(U_procs, V_procs, W_procs, P_procs, C_procs, Ax0_procs, Ay0_procs, Az0_procs, Tx0_procs, Ty0_procs, Tz0_procs)
        real(8), intent(out) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1), C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(out) :: Tx0_procs(1:NX, 1:NY, 1:N_procs), Ty0_procs(1:NX, 1:NY, 1:N_procs), Tz0_procs(1:NX, 1:NY, 1:N_procs)
        integer i, j, k
        character(8) str
        write(str, '(I8.8)') myrank

        open(100, file=trim(input_dir)//'Restart/'//str//'.bin', form='unformatted')
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    read(100) U_procs(i, j, k), V_procs(i, j, k), W_procs(i, j, k), P_procs(i, j, k), C_procs(:, i, j, k), &
                              Ax0_procs(i, j, k), Ay0_procs(i, j, k), Az0_procs(i, j, k), Tx0_procs(i, j, k), Ty0_procs(i, j, k), Tz0_procs(i, j, k)
                enddo
            enddo
        enddo
        close(100)
        if (myrank == 0) write(*, *) 'input_binary done'

        call MPI_Boundary(U_procs)
        call MPI_Boundary(V_procs)
        call MPI_Boundary(W_procs)
        call MPI_Boundary(P_procs)
        call C_MPI_Boundary(C_procs)
        call PBM(U_procs)
        call PBM(V_procs)
        call PBM(W_procs)
        call PBM(P_procs)
        call C_PBM(C_procs)
    end subroutine input_binary

    subroutine output_binary(U_procs, V_procs, W_procs, P_procs, C_procs, Ax_procs, Ay_procs, Az_procs, Tx_procs, Ty_procs, Tz_procs)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1), C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Tx_procs(1:NX, 1:NY, 1:N_procs), Ty_procs(1:NX, 1:NY, 1:N_procs), Tz_procs(1:NX, 1:NY, 1:N_procs)
        integer i, j, k
        character(8) str
        write(str, '(I8.8)') myrank

        open(100, file=trim(output_dir)//'Restart/'//str//'.bin', form='unformatted', status='replace')  ! 上書き保存する
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    write(100) U_procs(i, j, k), V_procs(i, j, k), W_procs(i, j, k), P_procs(i, j, k), C_procs(:, i, j, k), &
                               Ax_procs(i, j, k), Ay_procs(i, j, k), Az_procs(i, j, k), Tx_procs(i, j, k), Ty_procs(i, j, k), Tz_procs(i, j, k)
                enddo
            enddo
        enddo
        close(100)
        if (myrank == 0) write(*, *) 'output_binary done'
    end subroutine output_binary

    subroutine vtk_binary(U_procs, V_procs, W_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        integer i, j, k
        real(8) D(3, 3), Omega(3), Scalar_procs(1:NX, 1:NY, 1:N_procs), Scalar(1:NX, 1:NY, 1:NZ)
        character(8) str
        character(8) key
        character(120) buffer
        character lf
        lf = char(10)  ! 改行コード

        if (step < 10000000) key = 'Qti'
        if (10000000 <= step .and. step < 20000000) key = 'QtiA'
        if (20000000 <= step .and. step < 30000000) key = 'QtiB'
        if (30000000 <= step) key = 'QtiC'

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
                    D(:, :) = D(:, :)*U_C/L_C

                    Omega(1) = (D(3, 2) - D(2, 3))
                    Omega(2) = (D(1, 3) - D(3, 1))
                    Omega(3) = (D(2, 1) - D(1, 2))
                    ! Scalar_procs(i, j, k) = sum(Omega**2)/2  ! Enstrophy
                    Scalar_procs(i, j, k) = D(2, 2)*D(3, 3) - D(3, 2)*D(2, 3) + D(1, 1)*D(2, 2) - D(2, 1)*D(1, 2) + D(1, 1)*D(3, 3) - D(3, 1)*D(1, 3)  ! Qti
                enddo
            enddo
        enddo
        call MPI_Gather(Scalar_procs(1, 1, 1), NX*NY*N_procs, MPI_REAL8, Scalar(1, 1, 1), NX*NY*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        
        if (myrank == 0) then
            write(str, '(I8.8)') mod(step, 10000000)/Gstep  ! 連番にするため
            open(10, file=trim(output_dir)//'Qti/'//trim(key)//'_'//str//'.vtk', status='replace', form='unformatted', &
                    action='write', access='stream', convert='big_endian')

            write(buffer,'(a)') '# vtk DataFile Version 3.0'//lf
            write(10) trim(buffer)
            write(buffer,'(a)') 'Qti'//lf
            write(10) trim(buffer)
            write(buffer,'(a)') 'BINARY'//lf
            write(10) trim(buffer)
            write(buffer,'(a)') 'DATASET STRUCTURED_POINTS'//lf
            write(10) trim(buffer)
            write(buffer,'(a, 3(1x, i4))') 'DIMENSIONS', NX, NY, NZ
            write(10) trim(buffer)
            write(buffer,'(a, 3(1x, i3))') lf//'ORIGIN', 0, 0, 0
            write(10) trim(buffer)
            write(buffer,'(a, 3(1x, f16.10))') lf//'SPACING', dX_C, dY_C, dZ_C
            write(10) trim(buffer)
            write(buffer,'(a, i10)') lf//'POINT_DATA ', NX*NY*NZ
            write(10) trim(buffer)
            write(buffer,'(a)') lf//'SCALARS '//'Qti'//' float'//lf  ! doubleだと倍精度
            write(10) trim(buffer)
            write(buffer,'(a)') 'LOOKUP_TABLE default'//lf
            write(10) trim(buffer)
            do k = 1, NZ
                do j = 1, NY
                    do i = 1, NX
                        write(10) real(Scalar(i, j, k))  ! 単精度で十分
                    enddo
                enddo
            enddo
            close(10)
        endif
    end subroutine vtk_binary

    subroutine mk_dir(outdir)
        implicit none
        character(*), intent(in) :: outdir
        character(256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dir

    subroutine rm_dir(outdir)
        implicit none
        character(*), intent(in) :: outdir
        character(256) command
        write(command, *) 'rm -r ', trim(outdir)
        call system(command)
    end subroutine rm_dir

end module smac

module fft
    use smac
    implicit none
    include 'fftw3.f'
    integer(8) plan1, plan2, plan3, plan4, plan5, plan6
    real(8) Re1(1:NX, 1:NY), Re5(1:Nstep)
    complex(8) Im1(1:NX/2+1, 1:NY), Im2(1:NZ), Im3(1:NZ), Im5(1:Nstep/2+1)
contains
    subroutine fft_init
        call dfftw_plan_dft_r2c_2d(plan1, NX, NY, Re1, Im1, FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(plan2, NZ, Im2, Im3, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(plan3, NZ, Im3, Im2, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(plan4, NX, NY, Im1, Re1, FFTW_ESTIMATE)
        call dfftw_plan_dft_r2c_1d(plan5, Nstep, Re5, Im5, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_1d(plan6, Nstep, Im5, Re5, FFTW_ESTIMATE)
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

    ! 時系列方向のフーリエ変換用
    subroutine fftr2c_1d(input, output)
        real(8), intent(in) :: input(1:Nstep)
        complex(8), intent(out) :: output(1:Nstep/2+1)
        Re5(:) = input(:)
        call dfftw_execute(plan5, Re5, Im5)
        output(:) = Im5(:)
    end subroutine fftr2c_1d

    subroutine fftc2r_1d(input, output)
        complex(8), intent(in) :: input(1:Nstep/2+1)
        real(8), intent(out) :: output(1:Nstep)
        Im5(:) = input(:)
        call dfftw_execute(plan6, Im5, Re5)
        output(:) = Re5(:)
    end subroutine fftc2r_1d

    subroutine fft_finalize
        call dfftw_destroy_plan(plan1)
        call dfftw_destroy_plan(plan2)
        call dfftw_destroy_plan(plan3)
        call dfftw_destroy_plan(plan4)
        call dfftw_destroy_plan(plan5)
        call dfftw_destroy_plan(plan6)
    end subroutine fft_finalize


    subroutine fft_solve(Q_procs, LHS_procs, Phi_procs)
        real(8), intent(in) :: Q_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: LHS_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        real(8), intent(out) :: Phi_procs(1:NX, 1:NY, 1:N_procs)
        complex(8) Q_hat_procs(1:NX/2+1, 1:NY, 1:N_procs), Phi_hat_procs(1:NX/2+1, 1:NY, 1:N_procs)
        complex(8) Q_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ), Phi_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        complex(8) Z_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        complex(8) T1_procs(1:NX/2+1, 1:NY_procs, 1:procs, 1:N_procs), T2_procs(1:NX/2+1, 1:NY_procs, 1:N_procs, 1:procs)
        integer i, j, k

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
                    if (abs(LHS_procs(i, j, k)) < 1.0d-14) then  ! マシン零以下
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

    end subroutine fft_solve


    subroutine fft_navier(U_procs, V_procs, W_procs, P_procs, Up_procs, Vp_procs, Wp_procs, Ax_procs, Ay_procs, Az_procs, &
                          Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                          Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, Fx_procs, Fy_procs, Fz_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Tx_procs(1:NX, 1:NY, 1:N_procs), Ty_procs(1:NX, 1:NY, 1:N_procs), Tz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Tx0_procs(1:NX, 1:NY, 1:N_procs), Ty0_procs(1:NX, 1:NY, 1:N_procs), Tz0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Fx_procs(1:NX, 1:NY, 1:N_procs), Fy_procs(1:NX, 1:NY, 1:N_procs), Fz_procs(1:NX, 1:NY, 1:N_procs)
        integer, intent(in) :: step
        real(8) Q_procs(1:NX, 1:NY, 1:N_procs)
        real(8) LHS_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        integer i, j, k

        if (step==1 .and. input_type==0) then  ! 1ステップ目のみ例外処理
            Ax0_procs(:, :, :) = Ax_procs(:, :, :)
            Ay0_procs(:, :, :) = Ay_procs(:, :, :)
            Az0_procs(:, :, :) = Az_procs(:, :, :)
            Tx0_procs(:, :, :) = Tx_procs(:, :, :)
            Ty0_procs(:, :, :) = Ty_procs(:, :, :)
            Tz0_procs(:, :, :) = Tz_procs(:, :, :)
        endif

        ! 左辺の定数部分の計算
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    LHS_procs(i, j, k) = 1.0d0 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
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
                                      +dt*(Tx_procs(i, j, k) + Tx0_procs(i, j, k))/2.0d0 &
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
                                      +dt*(Ty_procs(i, j, k) + Ty0_procs(i, j, k))/2.0d0 &
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
                                       +dt*(Tz_procs(i, j, k) + Tz0_procs(i, j, k))/2.0d0 &
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
        Tx0_procs(:, :, :) = Tx_procs(:, :, :)
        Ty0_procs(:, :, :) = Ty_procs(:, :, :)
        Tz0_procs(:, :, :) = Tz_procs(:, :, :)
    end subroutine fft_navier

    subroutine fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
        real(8), intent(in) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: Phi_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8) Q_procs(1:NX, 1:NY, 1:N_procs)
        real(8) LHS_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        integer i, j, k

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
                                      -dt*beta/Re/2.0d0*((Phi_procs(i-1, j, k)-2*Phi_procs(i, j, k)+Phi_procs(i+1, j, k))/dX**2 &
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
        complex(8), intent(out) :: Q_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        complex(8) Q_hat_procs(1:NX/2+1, 1:NY, 1:N_procs)
        complex(8) Z_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        complex(8) T1_procs(1:NX/2+1, 1:NY_procs, 1:procs, 1:N_procs), T2_procs(1:NX/2+1, 1:NY_procs, 1:N_procs, 1:procs)
        integer i, j, k

        ! x,y方向にfft
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

    end subroutine fft_forward

    subroutine fft_z(Q_procs, Q_hat_procs)
        real(8), intent(in) :: Q_procs(1:NX, 1:NY, 1:N_procs)
        complex(8), intent(out) :: Q_hat_procs(1:NX, 1:NY_procs, 1:NZ)
        complex(8) Z_procs(1:NX, 1:NY_procs, 1:NZ)
        complex(8) T1_procs(1:NX, 1:NY_procs, 1:procs, 1:N_procs), T2_procs(1:NX, 1:NY_procs, 1:N_procs, 1:procs)
        integer i, j

        ! 軸変換  z分割をy分割に変更
        T1_procs(:, :, :, :) = reshape(Q_procs, shape(T1_procs))
        T2_procs(:, :, :, :) = reshape(T1_procs, shape(T2_procs), order=[1, 2, 4, 3])
        call MPI_Alltoall(T2_procs(1, 1, 1, 1), NX*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, &
                          Z_procs(1, 1, 1), NX*NY_procs*N_procs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        ! z方向にfft
        do j = 1, NY_procs
            do i = 1, NX
                call fftc2c_1d_forward(Z_procs(i, j, :), Q_hat_procs(i, j, :))  ! c2cだからNZのままでok
            enddo
        enddo
        Q_hat_procs(:, :, :) = Q_hat_procs(:, :, :)/NZ

    end subroutine fft_z


    subroutine fft_backward(Phi_hat_hat_procs, Phi_procs)
        complex(8), intent(inout) :: Phi_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        real(8), intent(out) :: Phi_procs(1:NX, 1:NY, 1:N_procs)
        complex(8) Phi_hat_procs(1:NX/2+1, 1:NY, 1:N_procs)
        complex(8) Z_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        complex(8) T1_procs(1:NX/2+1, 1:NY_procs, 1:procs, 1:N_procs), T2_procs(1:NX/2+1, 1:NY_procs, 1:N_procs, 1:procs)
        integer i, j, k

        ! z方向に逆fft
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

    end subroutine fft_backward


    subroutine scale_vtk(U_procs, V_procs, W_procs, step)  ! 渦のスケール分解
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        complex(8) U_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ), V_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ), W_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        complex(8) U_hat_scale_procs(3, 1:NX/2+1, 1:NY_procs, 1:NZ),V_hat_scale_procs(3, 1:NX/2+1, 1:NY_procs, 1:NZ), W_hat_scale_procs(3, 1:NX/2+1, 1:NY_procs, 1:NZ)
        real(8) K_abs_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        real(8) U_tmp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_tmp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_tmp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k, index
        real(8) k_index(4)
        real(8) kx, ky, kz

        ! 速度場をフーリエ変換
        call fft_forward(U_procs(1:NX, 1:NY, 1:N_procs), U_hat_hat_procs)
        call fft_forward(V_procs(1:NX, 1:NY, 1:N_procs), V_hat_hat_procs)
        call fft_forward(W_procs(1:NX, 1:NY, 1:N_procs), W_hat_hat_procs)

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

        ! バンドパスフィルタ
        U_hat_scale_procs(:, :, :, :) = 0.0d0
        V_hat_scale_procs(:, :, :, :) = 0.0d0
        W_hat_scale_procs(:, :, :, :) = 0.0d0
        k_index = [1.0d0, 3.0d0, 9.0d0, 27.0d0]  ! 波数の範囲を選択
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    do index = 1, 3
                        if (k_index(index) <= K_abs_procs(i, j, k) .and. K_abs_procs(i, j, k) < k_index(index+1)) then
                            U_hat_scale_procs(index, i, j, k) = U_hat_hat_procs(i, j, k)
                            V_hat_scale_procs(index, i, j, k) = V_hat_hat_procs(i, j, k)
                            W_hat_scale_procs(index, i, j, k) = W_hat_hat_procs(i, j, k)
                        endif
                    enddo
                enddo
            enddo
        enddo
        do index = 1, 3
            call fft_backward(U_hat_scale_procs(index, 1:NX/2+1, 1:NY_procs, 1:NZ), U_tmp_procs(1:NX, 1:NY, 1:N_procs))
            call fft_backward(V_hat_scale_procs(index, 1:NX/2+1, 1:NY_procs, 1:NZ), V_tmp_procs(1:NX, 1:NY, 1:N_procs))
            call fft_backward(W_hat_scale_procs(index, 1:NX/2+1, 1:NY_procs, 1:NZ), W_tmp_procs(1:NX, 1:NY, 1:N_procs))
            call MPI_Boundary(U_tmp_procs)
            call MPI_Boundary(V_tmp_procs)
            call MPI_Boundary(W_tmp_procs)
            call PBM(U_tmp_procs)
            call PBM(V_tmp_procs)
            call PBM(W_tmp_procs)
            call vtk_binary(U_tmp_procs, V_tmp_procs, W_tmp_procs, step + index*10000000)
        enddo

    end subroutine scale_vtk
end module fft


module ibm
    use smac
    use fft
    implicit none
    ! IBM用パラメータ
    integer NC, NL
    real(8) dV
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
        integer index
        integer NL_lap  ! 円柱内部の殻数
        integer, allocatable :: NL_shell(:)  ! 殻に対する外力点の数

        if (mod(ibm_type, 10) == 1) NC = 1  ! 円柱の数
        if (mod(ibm_type, 10) == 2) NC = 4

        if (ibm_type/10 == 1) then  ! 外力点が円柱表面のみの場合
            NL = nint(PI*DC/dX)
            dV = PI*DC/NL * dX * dX
        endif
        
        if (ibm_type/10 == 2) then  ! 外力点が円柱内部にもある場合
            NL_lap = int(DC/(2*dX))  ! 円柱内部の殻数
            allocate(NL_shell(0:NL_lap))
            do i = 0, NL_lap
                NL_shell(i) = nint(PI*(DC-2*i*dX)/dX)
            enddo
            if (NL_shell(NL_lap) == 0) NL_shell(NL_lap) = 1  ! 中心に外力点を1つ置く
            NL = sum(NL_shell)
            dV = PI*(DC+dX)**2/4 / NL * dX
            if (myrank == 0) write(*, '(a, 100I5)') 'NL_shell:  ', NL_shell(:)
        endif

        ! コメント
        if (dX /= dY .or. dX /= dZ) stop 'dX, dY and dZ are not equal'
        if (myrank == 0) then
            write(*, '(a, I5)') 'NL   :', NL
            ! write(*, '(a, E12.4)') 'dX**3:', dX**3
            ! write(*, '(a, E12.4)') 'dV   :', dV
        endif

        ! 配列のallocate
        N_procs = NZ / procs
        allocate(X_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Y_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Z_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Xc_procs(1:NC, 1:NL, 1:N_procs), Yc_procs(1:NC, 1:NL, 1:N_procs), Zc_procs(1:NC, 1:NL, 1:N_procs))
        allocate(Uc_procs(1:NC, 1:NL, 1:N_procs), Vc_procs(1:NC, 1:NL, 1:N_procs), Wc_procs(1:NC, 1:NL, 1:N_procs))
        allocate(Ua_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Va_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wa_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
        allocate(Fxc_procs(1:NC, 1:NL, 1:N_procs), Fyc_procs(1:NC, 1:NL, 1:N_procs), Fzc_procs(1:NC, 1:NL, 1:N_procs))
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

        ! 外力点の定義
        Uc_procs(:, :, :) = 0.0d0
        Vc_procs(:, :, :) = 0.0d0
        Wc_procs(:, :, :) = 0.0d0
        if (ibm_type == 11) then
            do j = 1, NL
                Xc_procs(1, j, :) = 2*PI/4 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(1, j, :) = 2*PI/2 + DC/2 * sin(2*PI/NL*j)
            enddo
        endif
        if (ibm_type == 12) then
            do j = 1, NL
                Xc_procs(1, j, :) = 2*PI/4   + DC/2 * cos(2*PI/NL*j)
                Yc_procs(1, j, :) = 2*PI/4   + DC/2 * sin(2*PI/NL*j)
                Xc_procs(2, j, :) = 2*PI*3/4 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(2, j, :) = 2*PI/4   + DC/2 * sin(2*PI/NL*j)
                Xc_procs(3, j, :) = 2*PI/4   + DC/2 * cos(2*PI/NL*j)
                Yc_procs(3, j, :) = 2*PI*3/4 + DC/2 * sin(2*PI/NL*j)
                Xc_procs(4, j, :) = 2*PI*3/4 + DC/2 * cos(2*PI/NL*j)
                Yc_procs(4, j, :) = 2*PI*3/4 + DC/2 * sin(2*PI/NL*j)
                Uc_procs(1, j, :) = sin(2*PI/NL*j)
                Vc_procs(1, j, :) = -cos(2*PI/NL*j)
                Uc_procs(2, j, :) = -sin(2*PI/NL*j)
                Vc_procs(2, j, :) = cos(2*PI/NL*j)
                Uc_procs(3, j, :) = -sin(2*PI/NL*j)
                Vc_procs(3, j, :) = cos(2*PI/NL*j)
                Uc_procs(4, j, :) = sin(2*PI/NL*j)
                Vc_procs(4, j, :) = -cos(2*PI/NL*j)
            enddo
        endif
        if (ibm_type == 21) then
            index = 0
            do i = 0, NL_lap
                do j = 1, NL_shell(i)
                    index = index + 1
                    Xc_procs(1, index, :) = 2*PI/2 + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc_procs(1, index, :) = 2*PI/2 + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Uc_procs(1, index, :) =  (DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc_procs(1, index, :) = -(DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                enddo
            enddo
        endif
        if (ibm_type == 22) then
            index = 0
            do i = 0, NL_lap
                do j = 1, NL_shell(i)
                    index = index + 1
                    Xc_procs(1, index, :) = 2*PI/4   + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc_procs(1, index, :) = 2*PI/4   + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Xc_procs(2, index, :) = 2*PI*3/4 + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc_procs(2, index, :) = 2*PI/4   + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Xc_procs(3, index, :) = 2*PI/4   + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc_procs(3, index, :) = 2*PI*3/4 + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Xc_procs(4, index, :) = 2*PI*3/4 + (DC-2*i*dX)/2 * cos(2*PI/NL_shell(i)*j)
                    Yc_procs(4, index, :) = 2*PI*3/4 + (DC-2*i*dX)/2 * sin(2*PI/NL_shell(i)*j)
                    Uc_procs(1, index, :) =  (DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)  ! 正しい
                    Vc_procs(1, index, :) = -(DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                    Uc_procs(2, index, :) = -(DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc_procs(2, index, :) =  (DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                    Uc_procs(3, index, :) = -(DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc_procs(3, index, :) =  (DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                    Uc_procs(4, index, :) =  (DC-2*i*dX)/DC * sin(2*PI/NL_shell(i)*j)
                    Vc_procs(4, index, :) = -(DC-2*i*dX)/DC * cos(2*PI/NL_shell(i)*j)
                enddo
            enddo
        endif
        do k = 1, N_procs
            Zc_procs(:, :, k) = (myrank*N_procs + k-0.25d0)*dZ
        enddo
        ! Xc_procs(:, :, :) = Xc_procs(:, :, :) + 10d-10*dX
        ! Yc_procs(:, :, :) = Yc_procs(:, :, :) + 10d-10*dY
        ! Zc_procs(:, :, :) = Zc_procs(:, :, :) + 10d-10*dZ

    end subroutine ibm_init

    subroutine ibm_preliminary(U_procs, V_procs, W_procs, P_procs, Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                               Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, &
                               Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, Up_procs, Vp_procs, Wp_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: P_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: Ax_procs(1:NX, 1:NY, 1:N_procs), Ay_procs(1:NX, 1:NY, 1:N_procs), Az_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Ax0_procs(1:NX, 1:NY, 1:N_procs), Ay0_procs(1:NX, 1:NY, 1:N_procs), Az0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Tx_procs(1:NX, 1:NY, 1:N_procs), Ty_procs(1:NX, 1:NY, 1:N_procs), Tz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(inout) :: Tx0_procs(1:NX, 1:NY, 1:N_procs), Ty0_procs(1:NX, 1:NY, 1:N_procs), Tz0_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(out) :: Ua_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Va_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wa_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(inout) :: Fxc_procs(1:NC, 1:NL, 1:N_procs), Fyc_procs(1:NC, 1:NL, 1:N_procs), Fzc_procs(1:NC, 1:NL, 1:N_procs)
        real(8), intent(out) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        integer i, j, k
        
        if (step==1 .and. input_type==0) then  ! 1ステップ目のみ例外処理
            Ax0_procs(:, :, :) = Ax_procs(:, :, :)
            Ay0_procs(:, :, :) = Ay_procs(:, :, :)
            Az0_procs(:, :, :) = Az_procs(:, :, :)
            Tx0_procs(:, :, :) = Tx_procs(:, :, :)
            Ty0_procs(:, :, :) = Ty_procs(:, :, :)
            Tz0_procs(:, :, :) = Tz_procs(:, :, :)
        endif

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Ua_procs(i, j, k) = U_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i+1, j, k))/dX &
                                       +dt*(3.0d0*Ax_procs(i, j, k) - Ax0_procs(i, j, k))/2.0d0 &
                                       +dt*Bx_procs(i, j, k) &
                                       +dt*(Tx_procs(i, j, k) + Tx0_procs(i, j, k))/2.0d0

                    Va_procs(i, j, k) = V_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j+1, k))/dY &
                                       +dt*(3.0d0*Ay_procs(i, j, k) - Ay0_procs(i, j, k))/2.0d0 &
                                       +dt*By_procs(i, j, k) &
                                       +dt*(Ty_procs(i, j, k) + Ty0_procs(i, j, k))/2.0d0

                    Wa_procs(i, j, k) = W_procs(i, j, k) - dt*(-P_procs(i, j, k)+P_procs(i, j, k+1))/dZ &
                                       +dt*(3.0d0*Az_procs(i, j, k) - Az0_procs(i, j, k))/2.0d0 &
                                       +dt*Bz_procs(i, j, k) &
                                       +dt*(Tz_procs(i, j, k) + Tz0_procs(i, j, k))/2.0d0
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
        Tx0_procs(:, :, :) = Tx_procs(:, :, :)
        Ty0_procs(:, :, :) = Ty_procs(:, :, :)
        Tz0_procs(:, :, :) = Tz_procs(:, :, :)

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
        real(8), intent(inout) :: Fxc_procs(1:NC, 1:NL, 1:N_procs), Fyc_procs(1:NC, 1:NL, 1:N_procs), Fzc_procs(1:NC, 1:NL, 1:N_procs)
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
                                                                      Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY - 0.5d0), int(Yc_procs(l, m, n)/dY - 0.5d0) + 2
                            do i = int(Xc_procs(l, m, n)/dX), int(Xc_procs(l, m, n)/dX) + 2
                                fytmp_procs(i, j, k) = fytmp_procs(i, j, k) &
                                               + Fyc_procs(l, m, n) * delta(X_procs(i, j, k), Y_procs(i, j, k)+0.5d0*dY, Z_procs(i, j, k), &
                                                                      Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                    do k = n-1, n+1
                        do j = int(Yc_procs(l, m, n)/dY) , int(Yc_procs(l, m, n)/dY) + 2
                            do i = int(Xc_procs(l, m, n)/dX), int(Xc_procs(l, m, n)/dX) + 2
                                fztmp_procs(i, j, k) = fztmp_procs(i, j, k) &
                                               + Fzc_procs(l, m, n) * delta(X_procs(i, j, k), Y_procs(i, j, k), Z_procs(i, j, k)+0.5d0*dZ, &
                                                                      Xc_procs(l, m, n), Yc_procs(l, m, n), Zc_procs(l, m, n)) * dX**3
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! 通信
        call IBM_MPI_Boundary(fxtmp_procs, fxint_procs)
        call IBM_MPI_Boundary(fytmp_procs, fyint_procs)
        call IBM_MPI_Boundary(fztmp_procs, fzint_procs)

    end subroutine ibm_Helmholtz


    subroutine IBM_MPI_Boundary(ftmp_procs, fint_procs)
        real(8), intent(inout) :: ftmp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(out) :: fint_procs(1:NX, 1:NY, 1:N_procs)
        real(8) :: f1_procs(0:NX+1, 0:NY+1)
        real(8) :: f2_procs(0:NX+1, 0:NY+1)

        ! xy方向周期境界
        ftmp_procs(1, :, :) = ftmp_procs(1, :, :) + ftmp_procs(NX+1, :, :)
        ftmp_procs(NX, :, :) = ftmp_procs(NX, :, :) + ftmp_procs(0, :, :)
        ftmp_procs(:, 1, :) = ftmp_procs(:, 1, :) + ftmp_procs(:, NY+1, :)
        ftmp_procs(:, NY, :) = ftmp_procs(:, NY, :) + ftmp_procs(:, 0, :)
        fint_procs(:, :, :) = ftmp_procs(1:NX, 1:NY, 1:N_procs)

        ! z方向通信
        ! 次のランクへ送信
        call MPI_Isend(ftmp_procs(0, 0, N_procs+1), (NX+2)*(NY+2), MPI_REAL8, next_rank, 1, MPI_COMM_WORLD, req1s, ierr)
        call MPI_Irecv(f1_procs(0, 0), (NX+2)*(NY+2), MPI_REAL8, former_rank, 1, MPI_COMM_WORLD, req1r, ierr)  ! 一番上に足すやつ
        ! 前のランクへ送信
        call MPI_Isend(ftmp_procs(0, 0, 0), (NX+2)*(NY+2), MPI_REAL8, former_rank, 2, MPI_COMM_WORLD, req2s, ierr)
        call MPI_Irecv(f2_procs(0, 0), (NX+2)*(NY+2), MPI_REAL8, next_rank, 2, MPI_COMM_WORLD, req2r, ierr)  ! 一番下に足すやつ
        ! 待機
        call MPI_Wait(req1s, sta1s, ierr)
        call MPI_Wait(req1r, sta1r, ierr)
        call MPI_Wait(req2s, sta2s, ierr)
        call MPI_Wait(req2r, sta2r, ierr)
        fint_procs(:, :, 1) = fint_procs(:, :, 1) + f1_procs(1:NX, 1:NY)
        fint_procs(:, :, N_procs) = fint_procs(:, :, N_procs) + f2_procs(1:NX, 1:NY)

    end subroutine IBM_MPI_Boundary

    subroutine ibm_predict(Ua_procs, Va_procs, Wa_procs, fxint_procs, fyint_procs, fzint_procs, Bx_procs, By_procs, Bz_procs, Up_procs, Vp_procs, Wp_procs)
        real(8), intent(in) :: Ua_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Va_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wa_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: fxint_procs(1:NX, 1:NY, 1:N_procs), fyint_procs(1:NX, 1:NY, 1:N_procs), fzint_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: Bx_procs(1:NX, 1:NY, 1:N_procs), By_procs(1:NX, 1:NY, 1:N_procs), Bz_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(out) :: Up_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Vp_procs(0:NX+1, 0:NY+1, 0:N_procs+1), Wp_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        integer i, j, k
        real(8) Q_procs(1:NX, 1:NY, 1:N_procs), LHS_procs(1:NX/2+1, 1:NY_procs, 1:NZ)

        ! 左辺の定数部分の計算
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX/2+1
                    LHS_procs(i, j, k) = 1.0d0 + dt*beta/Re/2.0d0*(2*(1-cos((i-1)*dX))/dX**2 &
                                                             + 2*(1-cos((myrank*NY_procs + j-1)*dY))/dY**2 &
                                                             + 2*(1-cos((k-1)*dZ))/dZ**2)
                enddo
            enddo
        enddo

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = Ua_procs(i, j, k) + dt*fxint_procs(i, j, k) - dt*Bx_procs(i, j, k)/2.0d0
                enddo
            enddo
        enddo
        call fft_solve(Q_procs, LHS_procs, Up_procs(1:NX, 1:NY, 1:N_procs))

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = Va_procs(i, j, k) + dt*fyint_procs(i, j, k) - dt*By_procs(i, j, k)/2.0d0
                enddo
            enddo
        enddo
        call fft_solve(Q_procs, LHS_procs, Vp_procs(1:NX, 1:NY, 1:N_procs))

        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Q_procs(i, j, k) = Wa_procs(i, j, k) + dt*fzint_procs(i, j, k) - dt*Bz_procs(i, j, k)/2.0d0
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


    subroutine ave_init(U_ave_procs, V_ave_procs, W_ave_procs, U_rms_procs, V_rms_procs, W_rms_procs, Energy_map_procs, U_step_procs, V_step_procs, W_step_procs)
        real(8), allocatable :: U_ave_procs(:, :, :), V_ave_procs(:, :, :), W_ave_procs(:, :, :)
        real(8), allocatable :: U_rms_procs(:, :, :), V_rms_procs(:, :, :), W_rms_procs(:, :, :)
        real(8), allocatable :: Energy_map_procs(:, :, :)
        real(8), allocatable :: U_step_procs(:, :), V_step_procs(:, :), W_step_procs(:, :)

        allocate(U_ave_procs(1:NX, 1:NY, 1:N_procs), V_ave_procs(1:NX, 1:NY, 1:N_procs), W_ave_procs(1:NX, 1:NY, 1:N_procs))
        allocate(U_rms_procs(1:NX, 1:NY, 1:N_procs), V_rms_procs(1:NX, 1:NY, 1:N_procs), W_rms_procs(1:NX, 1:NY, 1:N_procs))
        allocate(Energy_map_procs(1:NX, 1:NY_procs, 0:NZ))
        allocate(U_step_procs(1:N_procs, 1:Nstep), V_step_procs(1:N_procs, 1:Nstep), W_step_procs(1:N_procs, 1:Nstep))
        ! あるx, yに対し、z座標における速度を全ステップで保存
        U_ave_procs(:, :, :) = 0.0d0
        V_ave_procs(:, :, :) = 0.0d0
        W_ave_procs(:, :, :) = 0.0d0
        U_rms_procs(:, :, :) = 0.0d0
        V_rms_procs(:, :, :) = 0.0d0
        W_rms_procs(:, :, :) = 0.0d0
        Energy_map_procs(:, :, :) = 0.0d0
        U_step_procs(:, :) = 0.0d0
        V_step_procs(:, :) = 0.0d0
        W_step_procs(:, :) = 0.0d0
    end subroutine ave_init

    subroutine ave_mean(U_ave_procs, V_ave_procs, W_ave_procs)
        real(8), intent(inout) :: U_ave_procs(1:NX, 1:NY, 1:N_procs), V_ave_procs(1:NX, 1:NY, 1:N_procs), W_ave_procs(1:NX, 1:NY, 1:N_procs)
        real(8) U_ave_mean(1:NX, 1:NY), V_ave_mean(1:NX, 1:NY), W_ave_mean(1:NX, 1:NY)
        real(8) U_ave_tmp(1:NX, 1:NY), V_ave_tmp(1:NX, 1:NY), W_ave_tmp(1:NX, 1:NY)
        integer i, j, k
        do j = 1, NY
            do i = 1, NX
                U_ave_mean(i, j) = sum(U_ave_procs(i, j, :))/NZ
                V_ave_mean(i, j) = sum(V_ave_procs(i, j, :))/NZ
                W_ave_mean(i, j) = sum(W_ave_procs(i, j, :))/NZ
            enddo
        enddo

        call MPI_Allreduce(U_ave_mean, U_ave_tmp, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_ave_mean, V_ave_tmp, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_ave_mean, W_ave_tmp, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        do k = 1, N_procs
            U_ave_procs(:, :, k) = U_ave_tmp(:, :)
            V_ave_procs(:, :, k) = V_ave_tmp(:, :)
            W_ave_procs(:, :, k) = W_ave_tmp(:, :)
        enddo
        
    end subroutine ave_mean


    subroutine ibm_var(U_rms_procs, V_rms_procs, W_rms_procs)
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

        if (myrank == 0) write(*, *) 'ibm_var:', var

    end subroutine ibm_var


    subroutine ibm_vtk(Xc_procs, Yc_procs, Zc_procs)
        real(8), intent(in) :: Xc_procs(1:NC, 1:NL, 1:N_procs), Yc_procs(1:NC, 1:NL, 1:N_procs), Zc_procs(1:NC, 1:NL, 1:N_procs)
        real(8) Xc(1:NC, 1:NL, 1:NZ), Yc(1:NC, 1:NL, 1:NZ), Zc(1:NC, 1:NL, 1:NZ)
        integer NR, l, m, n, o
        character(120) buffer
        character(8) str
        character lf
        lf = char(10)  ! 改行コード
        NR = nint(PI*DC/dX) + 1  ! 分割数+1

        call MPI_Gather(Xc_procs(1, 1, 1), NC*NL*N_procs, MPI_REAL8, Xc(1, 1, 1), NC*NL*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(Yc_procs(1, 1, 1), NC*NL*N_procs, MPI_REAL8, Yc(1, 1, 1), NC*NL*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(Zc_procs(1, 1, 1), NC*NL*N_procs, MPI_REAL8, Zc(1, 1, 1), NC*NL*N_procs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        if (myrank == 0) then
            do l = 1, NC
                write(str, '(I8.8)') l
                open(10, file=trim(output_dir)//'Cylinder/'//'Cylinder_'//str//'.vtk', status='replace', form='unformatted', &
                        action='write', access='stream', convert='big_endian')

                write(buffer,'(a)') '# vtk DataFile Version 3.0'//lf
                write(10) trim(buffer)
                write(buffer,'(a)') 'Cylinder Data'//lf  ! hogeはなんでもいい
                write(10) trim(buffer)
                write(buffer,'(a)') 'BINARY'//lf
                write(10) trim(buffer)
                ! 座標値をvtkファイルに書き出す
                write(buffer,'(a)') 'DATASET STRUCTURED_GRID'//lf
                write(10) trim(buffer)
                write(buffer,'(a, 3(1x, i4))') 'DIMENSIONS', NR, NR, NZ
                write(10) trim(buffer)
                write(buffer,'(a, i10, a)') lf//'POINTS ', NR*NR*NZ, 'float'//lf  ! doubleだと倍精度
                write(10) trim(buffer)
                do n = 1, NZ
                    do o = 1, NR  ! 3次元で円柱描くとき必要
                        do m = 1, NR-1  ! NLではなくNR-1まで
                            write(10) real(Xc(l, m, n)*L_C), real(Yc(l, m, n)*L_C), real(Zc(l, m, n)*L_C)
                        enddo
                        write(10) real(Xc(l, 1, n)*L_C), real(Yc(l, 1, n)*L_C), real(Zc(l, 1, n)*L_C)
                    enddo
                enddo

                ! 点データをvtkファイルに書き出す
                write(buffer, '(a, i10)') lf//'POINT_DATA ', NR*NR*NZ
                write(10) trim(buffer)
                write(buffer, '(a)') lf//'SCALARS Cylinder int 1'//lf
                write(10) trim(buffer)
                write(buffer, '(a)') 'LOOKUP_TABLE default'//lf
                write(10) trim(buffer)

                do n = 1, NR*NR*NZ
                    write(10) 1  ! n-1
                enddo
                close(10)
            enddo
        endif
    end subroutine ibm_vtk
end module ibm


module data
    use smac
    use fft
    use ibm
    implicit none
contains
    subroutine log_progress(U_procs, V_procs, W_procs, C_procs, step, time1)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1), time1
        integer, intent(in) :: step
        real(8) time0
        integer time2, hour, min, sec
        real(8) U_reg_procs(1:NX, 1:NY, 1:N_procs), V_reg_procs(1:NX, 1:NY, 1:N_procs), W_reg_procs(1:NX, 1:NY, 1:N_procs)
        real(8) K_energy
        real(8) trC_procs(1:NX, 1:NY, 1:N_procs)
        real(8) trC_tmp, trC_mean, trC_std, trC_max, trC_min
        real(8) U_tmp, V_tmp, W_tmp, U_tmp_sum, V_tmp_sum, W_tmp_sum
        integer i, j, k

        ! 運動エネルギー
        do k = 1, N_procs  ! レギュラー格子に直して計算
            do j = 1, NY
                do i = 1, NX
                    U_reg_procs(i, j, k) = (U_procs(i, j, k) + U_procs(i-1, j, k))*U_C/2.0d0
                    V_reg_procs(i, j, k) = (V_procs(i, j, k) + V_procs(i, j-1, k))*U_C/2.0d0
                    W_reg_procs(i, j, k) = (W_procs(i, j, k) + W_procs(i, j, k-1))*U_C/2.0d0
                enddo
            enddo
        enddo
        U_tmp = sum(U_reg_procs**2)/(NX*NY*NZ)
        V_tmp = sum(V_reg_procs**2)/(NX*NY*NZ)
        W_tmp = sum(W_reg_procs**2)/(NX*NY*NZ)
        call MPI_Allreduce(U_tmp, U_tmp_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_tmp, V_tmp_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_tmp, W_tmp_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        K_energy = (U_tmp_sum + V_tmp_sum + W_tmp_sum)/2.0d0
        if (K_energy*0.0d0 /= 0.0d0) then
            if (myrank == 0) write(*, *) 'step:', step, 'NaN value'
            stop 'NaN value'  ! NaNの判定
        endif

        ! CFL条件
        U_tmp = maxval(U_procs(1:NX, 1:NY, 1:N_procs))
        V_tmp = maxval(V_procs(1:NX, 1:NY, 1:N_procs))
        W_tmp = maxval(W_procs(1:NX, 1:NY, 1:N_procs))
        call MPI_Allreduce(U_tmp, U_tmp_sum, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_tmp, V_tmp_sum, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_tmp, W_tmp_sum, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

        ! 伸長長さ
        trC_procs(:, :, :) = C_procs(1, 1:NX, 1:NY, 1:N_procs)+C_procs(4, 1:NX, 1:NY, 1:N_procs)+C_procs(6, 1:NX, 1:NY, 1:N_procs)
        trC_tmp = sum(trC_procs)/(NX*NY*NZ)
        call MPI_Allreduce(trC_tmp, trC_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        trC_tmp = sum((trC_procs - trC_mean)**2)/(NX*NY*NZ)
        call MPI_Allreduce(trC_tmp, trC_std, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        trC_std = sqrt(trC_std)
        trC_tmp = maxval(trC_procs(1:NX, 1:NY, 1:N_procs))
        call MPI_Allreduce(trC_tmp, trC_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
        trC_tmp = minval(trC_procs(1:NX, 1:NY, 1:N_procs))
        call MPI_Allreduce(trC_tmp, trC_min, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
        
        if (myrank == 0) then
            time0 = MPI_Wtime()
            if (input_type /= 2) time2 = int((time0-time1)*(Nstep-step)/step)
            if (input_type == 2) time2 = int((time0-time1)*(2*Nstep-step)/step)
            hour = time2 / 3600
            min = mod(time2, 3600) / 60
            sec = mod(time2, 60)

            write(*, '(a, I7)', advance='no') 'step:', step
            write(*, '(a, e12.4)', advance='no') ' | K_energy:', K_energy
            write(*, '(a, F7.3)', advance='no') ' | CFL:', max(U_tmp_sum*dt/dX, V_tmp_sum*dt/dY, W_tmp_sum*dt/dZ)*6.0d0
            ! write(*, '(a, F8.3, a, F8.3)', advance='no') ' | trC:', trC_mean, ' +-', trC_std
            write(*, '(a, F8.3, a, F8.3)', advance='no') ' |', trC_min, ' < trC <', trC_max
            write(*, '(a, I3, a, I2, a, I2)', advance='no') ' | time_left:', hour, ':', min, ':', sec
            write(*, *) ''
        endif
    end subroutine log_progress


    subroutine all_time_dif(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, C_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: U_ave_procs(1:NX, 1:NY, 1:N_procs), V_ave_procs(1:NX, 1:NY, 1:N_procs), W_ave_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        integer i, j, k
        real(8) K_energy, K_staggered
        real(8) trC_procs(1:NX, 1:NY, 1:N_procs)
        real(8) Ep_procs(1:NX, 1:NY, 1:N_procs)
        real(8) U_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8) U_reg_procs(1:NX, 1:NY, 1:N_procs), V_reg_procs(1:NX, 1:NY, 1:N_procs), W_reg_procs(1:NX, 1:NY, 1:N_procs)
        real(8) U_grad_procs(1:NX, 1:NY, 1:N_procs), V_grad_procs(1:NX, 1:NY, 1:N_procs), W_grad_procs(1:NX, 1:NY, 1:N_procs)
        real(8) lambda, Re_lambda, lambda_epsilon
        real(8) tmp, tmp_sum, D(3, 3), S(3, 3)
        real(8) epsilon, L_kolmogorov, T_kolmogorov, U_kolmogorov
        real(8) trC_mean, trC_std
        real(8) U_tmp, V_tmp, W_tmp
        real(8) U_rms, V_rms, W_rms, U_all_rms
        real(8) U_grad, V_grad, W_grad, U_all_grad
        real(8) Ep_energy
        character(32) key
        if (sum(U_ave_procs**2) < 1.0d-10) key = 'same'  ! 平均を引かない
        if (sum(U_ave_procs**2) > 1.0d-10) key = 'diff'  ! 平均を引く
        
        ! 変動速度
        U_dif_procs(1:NX, 1:NY, 1:N_procs) = (U_procs(1:NX, 1:NY, 1:N_procs) - U_ave_procs(:, :, :)) * U_C
        V_dif_procs(1:NX, 1:NY, 1:N_procs) = (V_procs(1:NX, 1:NY, 1:N_procs) - V_ave_procs(:, :, :)) * U_C
        W_dif_procs(1:NX, 1:NY, 1:N_procs) = (W_procs(1:NX, 1:NY, 1:N_procs) - W_ave_procs(:, :, :)) * U_C
        call MPI_Boundary(U_dif_procs)
        call MPI_Boundary(V_dif_procs)
        call MPI_Boundary(W_dif_procs)
        call PBM(U_dif_procs)
        call PBM(V_dif_procs)
        call PBM(W_dif_procs)

        ! 変動速度の2乗
        do k = 1, N_procs  ! レギュラー格子に直して計算
            do j = 1, NY
                do i = 1, NX
                    U_reg_procs(i, j, k) = (U_dif_procs(i, j, k) + U_dif_procs(i-1, j, k))/2.0d0
                    V_reg_procs(i, j, k) = (V_dif_procs(i, j, k) + V_dif_procs(i, j-1, k))/2.0d0
                    W_reg_procs(i, j, k) = (W_dif_procs(i, j, k) + W_dif_procs(i, j, k-1))/2.0d0
                enddo
            enddo
        enddo
        U_tmp = sum(U_reg_procs**2)/(NX*NY*NZ)
        V_tmp = sum(V_reg_procs**2)/(NX*NY*NZ)
        W_tmp = sum(W_reg_procs**2)/(NX*NY*NZ)
        call MPI_Allreduce(U_tmp, U_rms, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_tmp, V_rms, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_tmp, W_rms, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! 変動速度の微分の2乗
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    U_grad_procs(i, j, k) = (U_dif_procs(i, j, k) - U_dif_procs(i-1, j, k))/dX_C
                    V_grad_procs(i, j, k) = (V_dif_procs(i, j, k) - V_dif_procs(i, j-1, k))/dY_C
                    W_grad_procs(i, j, k) = (W_dif_procs(i, j, k) - W_dif_procs(i, j, k-1))/dZ_C
                enddo
            enddo
        enddo
        U_tmp = sum(U_grad_procs**2)/(NX*NY*NZ)
        V_tmp = sum(V_grad_procs**2)/(NX*NY*NZ)
        W_tmp = sum(W_grad_procs**2)/(NX*NY*NZ)
        call MPI_Allreduce(U_tmp, U_grad, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_tmp, V_grad, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_tmp, W_grad, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! エネルギー散逸率
        tmp = 0.0d0
        do k = 1, N_procs  ! 変動速度で計算
            do j = 1, NY
                do i = 1, NX
                    D(1, 1) = (U_dif_procs(i, j, k) - U_dif_procs(i-1, j, k))/dX
                    D(1, 2) = (U_dif_procs(i, j+1, k) - U_dif_procs(i, j-1, k) + U_dif_procs(i-1, j+1, k) - U_dif_procs(i-1, j-1, k))/(4*dY)
                    D(1, 3) = (U_dif_procs(i, j, k+1) - U_dif_procs(i, j, k-1) + U_dif_procs(i-1, j, k+1) - U_dif_procs(i-1, j, k-1))/(4*dZ)
                    D(2, 1) = (V_dif_procs(i+1, j, k) - V_dif_procs(i-1, j, k) + V_dif_procs(i+1, j-1, k) - V_dif_procs(i-1, j-1, k))/(4*dX)
                    D(2, 2) = (V_dif_procs(i, j, k) - V_dif_procs(i, j-1, k))/dY
                    D(2, 3) = (V_dif_procs(i, j, k+1) - V_dif_procs(i, j, k-1) + V_dif_procs(i, j-1, k+1) - V_dif_procs(i, j-1, k-1))/(4*dZ)
                    D(3, 1) = (W_dif_procs(i+1, j, k) - W_dif_procs(i-1, j, k) + W_dif_procs(i+1, j, k-1) - W_dif_procs(i-1, j, k-1))/(4*dX)
                    D(3, 2) = (W_dif_procs(i, j+1, k) - W_dif_procs(i, j-1, k) + W_dif_procs(i, j+1, k-1) - W_dif_procs(i, j-1, k-1))/(4*dY)
                    D(3, 3) = (W_dif_procs(i, j, k) - W_dif_procs(i, j, k-1))/dZ
                    D(:, :) = D(:, :)/L_C
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

        ! エネルギー散逸率
        U_all_rms = (U_rms + V_rms + W_rms)/3.0d0
        U_all_grad = (U_grad + V_grad + W_grad)/3.0d0
        epsilon = nu/2.0d0*tmp_sum/(NX*NY*NZ)

        ! コルモゴロフ
        L_kolmogorov = (nu**3/epsilon)**0.25d0  ! 最小スケールの渦長さ
        T_kolmogorov = (nu/epsilon)**0.5d0  ! 最小スケールの渦時間
        U_kolmogorov = (nu*epsilon)**0.25d0  ! 最小スケールの渦速度

        ! テイラー長
        lambda = sqrt(U_all_rms/U_all_grad)  ! テイラー長(定義方法は色々ある)
        Re_lambda = sqrt(U_all_rms)*lambda/nu  ! テイラー長レイノルズ数
        lambda_epsilon = 10*nu*U_all_rms/lambda**2/epsilon  ! テイラー長とエネルギー散逸率の関係(デバッグ用)(本当は15)

        ! 変動速度の運動エネルギー
        K_energy = (U_rms + V_rms + W_rms)/2.0d0
        ! スタガード格子で計算
        tmp = sum(U_procs(1:NX, 1:NY, 1:N_procs)**2) + sum(V_procs(1:NX, 1:NY, 1:N_procs)**2) + sum(W_procs(1:NX, 1:NY, 1:N_procs)**2)
        tmp = tmp*U_C**2/(2.0d0*NX*NY*NZ)
        call MPI_Allreduce(tmp, K_staggered, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! 弾性エネルギー
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Ep_procs(i, j, k) = 0.5d0*(1.0d0-beta)/Re/Wi*(Lp**2.0d0-3.0d0)*log(f(C_procs(:, i, j, k)))
                enddo
            enddo
        enddo
        tmp = sum(Ep_procs(:, :, :))/(NX*NY*NZ)
        call MPI_Allreduce(tmp, Ep_energy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! 伸長長さ
        trC_procs(:, :, :) = C_procs(1, 1:NX, 1:NY, 1:N_procs)+C_procs(4, 1:NX, 1:NY, 1:N_procs)+C_procs(6, 1:NX, 1:NY, 1:N_procs)
        tmp = sum(trC_procs)/(NX*NY*NZ)
        call MPI_Allreduce(tmp, trC_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        tmp = sum((trC_procs - trC_mean)**2)/(NX*NY*NZ)
        call MPI_Allreduce(tmp, trC_std, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        trC_std = sqrt(trC_std)
        
        ! 時系列データの保存
        if (myrank == 0) then
            open(30, file = trim(output_dir)//'Log/time_'//trim(key)//'.d', position='append')
            write(30, *) step, step*dt_C, K_energy, K_staggered, lambda, Re_lambda, epsilon, L_kolmogorov, T_kolmogorov, U_kolmogorov, &
                            Ep_energy, trC_mean, trC_std
            close(30)
        endif

    end subroutine all_time_dif

    subroutine all_time_dif_map(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, C_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: U_ave_procs(1:NX, 1:NY, 1:N_procs), V_ave_procs(1:NX, 1:NY, 1:N_procs), W_ave_procs(1:NX, 1:NY, 1:N_procs)
        real(8), intent(in) :: C_procs(6, 0:NX+1, 0:NY+1, 0:N_procs+1)
        integer, intent(in) :: step
        integer i, j, k
        real(8) trC_procs(1:NX, 1:NY, 1:N_procs)
        real(8) Ep_procs(1:NX, 1:NY, 1:N_procs)
        real(8) U_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8) U_reg_procs(1:NX, 1:NY, 1:N_procs), V_reg_procs(1:NX, 1:NY, 1:N_procs), W_reg_procs(1:NX, 1:NY, 1:N_procs)
        real(8) U_grad_procs(1:NX, 1:NY, 1:N_procs), V_grad_procs(1:NX, 1:NY, 1:N_procs), W_grad_procs(1:NX, 1:NY, 1:N_procs)
        real(8) D(3, 3), S(3, 3)
        real(8) K_energy
        real(8) epsilon, L_kolmogorov, T_kolmogorov, U_kolmogorov
        real(8) lambda, Re_lambda, lambda_epsilon
        real(8) trC_mean, trC_std
        real(8) Ep_energy
        ! mapの定義
        real(8) K_energy_map(1:NX, 1:NY)
        real(8) tmp_map(1:NX, 1:NY), tmp_sum_map(1:NX, 1:NY)
        real(8) epsilon_map(1:NX, 1:NY), L_kolmogorov_map(1:NX, 1:NY), T_kolmogorov_map(1:NX, 1:NY), U_kolmogorov_map(1:NX, 1:NY)
        real(8) lambda_map(1:NX, 1:NY), Re_lambda_map(1:NX, 1:NY), lambda_epsilon_map(1:NX, 1:NY)
        real(8) U_tmp_map(1:NX, 1:NY), V_tmp_map(1:NX, 1:NY), W_tmp_map(1:NX, 1:NY)
        real(8) U_rms_map(1:NX, 1:NY), V_rms_map(1:NX, 1:NY), W_rms_map(1:NX, 1:NY), U_all_rms_map(1:NX, 1:NY)
        real(8) U_grad_map(1:NX, 1:NY), V_grad_map(1:NX, 1:NY), W_grad_map(1:NX, 1:NY), U_all_grad_map(1:NX, 1:NY)
        real(8) trC_mean_map(1:NX, 1:NY), trC_std_map(1:NX, 1:NY)
        real(8) Ep_energy_map(1:NX, 1:NY)
        real(8) mask_map(1:NX, 1:NY)
        character(8) str
        character(32) key
        if (sum(U_ave_procs**2) < 1.0d-10) key = 'same'  ! 平均を引かない
        if (sum(U_ave_procs**2) > 1.0d-10) key = 'diff'  ! 平均を引く
        
        ! 変動速度
        U_dif_procs(1:NX, 1:NY, 1:N_procs) = (U_procs(1:NX, 1:NY, 1:N_procs) - U_ave_procs(:, :, :)) * U_C
        V_dif_procs(1:NX, 1:NY, 1:N_procs) = (V_procs(1:NX, 1:NY, 1:N_procs) - V_ave_procs(:, :, :)) * U_C
        W_dif_procs(1:NX, 1:NY, 1:N_procs) = (W_procs(1:NX, 1:NY, 1:N_procs) - W_ave_procs(:, :, :)) * U_C
        call MPI_Boundary(U_dif_procs)
        call MPI_Boundary(V_dif_procs)
        call MPI_Boundary(W_dif_procs)
        call PBM(U_dif_procs)
        call PBM(V_dif_procs)
        call PBM(W_dif_procs)

        ! 変動速度の2乗
        do k = 1, N_procs  ! レギュラー格子に直して計算
            do j = 1, NY
                do i = 1, NX
                    U_reg_procs(i, j, k) = (U_dif_procs(i, j, k) + U_dif_procs(i-1, j, k))/2.0d0
                    V_reg_procs(i, j, k) = (V_dif_procs(i, j, k) + V_dif_procs(i, j-1, k))/2.0d0
                    W_reg_procs(i, j, k) = (W_dif_procs(i, j, k) + W_dif_procs(i, j, k-1))/2.0d0
                enddo
            enddo
        enddo
        do j = 1, NY
            do i = 1, NX
                U_tmp_map(i, j) = sum(U_reg_procs(i, j, :)**2)/NZ
                V_tmp_map(i, j) = sum(V_reg_procs(i, j, :)**2)/NZ
                W_tmp_map(i, j) = sum(W_reg_procs(i, j, :)**2)/NZ
            enddo
        enddo
        call MPI_Allreduce(U_tmp_map, U_rms_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_tmp_map, V_rms_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_tmp_map, W_rms_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! 変動速度の微分の2乗
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    U_grad_procs(i, j, k) = (U_dif_procs(i, j, k) - U_dif_procs(i-1, j, k))/dX_C
                    V_grad_procs(i, j, k) = (V_dif_procs(i, j, k) - V_dif_procs(i, j-1, k))/dY_C
                    W_grad_procs(i, j, k) = (W_dif_procs(i, j, k) - W_dif_procs(i, j, k-1))/dZ_C
                enddo
            enddo
        enddo
        do j = 1, NY
            do i = 1, NX
                U_tmp_map(i, j) = sum(U_grad_procs(i, j, :)**2)/NZ
                V_tmp_map(i, j) = sum(V_grad_procs(i, j, :)**2)/NZ
                W_tmp_map(i, j) = sum(W_grad_procs(i, j, :)**2)/NZ
            enddo
        enddo
        call MPI_Allreduce(U_tmp_map, U_grad_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_tmp_map, V_grad_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_tmp_map, W_grad_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! エネルギー散逸率
        tmp_map(:, :) = 0.0d0
        do k = 1, N_procs  ! 変動速度で計算
            do j = 1, NY
                do i = 1, NX
                    D(1, 1) = (U_dif_procs(i, j, k) - U_dif_procs(i-1, j, k))/dX
                    D(1, 2) = (U_dif_procs(i, j+1, k) - U_dif_procs(i, j-1, k) + U_dif_procs(i-1, j+1, k) - U_dif_procs(i-1, j-1, k))/(4*dY)
                    D(1, 3) = (U_dif_procs(i, j, k+1) - U_dif_procs(i, j, k-1) + U_dif_procs(i-1, j, k+1) - U_dif_procs(i-1, j, k-1))/(4*dZ)
                    D(2, 1) = (V_dif_procs(i+1, j, k) - V_dif_procs(i-1, j, k) + V_dif_procs(i+1, j-1, k) - V_dif_procs(i-1, j-1, k))/(4*dX)
                    D(2, 2) = (V_dif_procs(i, j, k) - V_dif_procs(i, j-1, k))/dY
                    D(2, 3) = (V_dif_procs(i, j, k+1) - V_dif_procs(i, j, k-1) + V_dif_procs(i, j-1, k+1) - V_dif_procs(i, j-1, k-1))/(4*dZ)
                    D(3, 1) = (W_dif_procs(i+1, j, k) - W_dif_procs(i-1, j, k) + W_dif_procs(i+1, j, k-1) - W_dif_procs(i-1, j, k-1))/(4*dX)
                    D(3, 2) = (W_dif_procs(i, j+1, k) - W_dif_procs(i, j-1, k) + W_dif_procs(i, j+1, k-1) - W_dif_procs(i, j-1, k-1))/(4*dY)
                    D(3, 3) = (W_dif_procs(i, j, k) - W_dif_procs(i, j, k-1))/dZ
                    D(:, :) = D(:, :)/L_C
                    S(1, 1) = (D(1, 1) + D(1, 1))
                    S(1, 2) = (D(1, 2) + D(2, 1))
                    S(1, 3) = (D(1, 3) + D(3, 1))
                    S(2, 1) = (D(2, 1) + D(1, 2))
                    S(2, 2) = (D(2, 2) + D(2, 2))
                    S(2, 3) = (D(2, 3) + D(3, 2))
                    S(3, 1) = (D(3, 1) + D(1, 3))
                    S(3, 2) = (D(3, 2) + D(2, 3))
                    S(3, 3) = (D(3, 3) + D(3, 3))
                    tmp_map(i, j) = tmp_map(i, j) + sum(S**2)
                enddo
            enddo
        enddo
        call MPI_Allreduce(tmp_map, tmp_sum_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        do j = 1, NY
            do i = 1, NX
                ! エネルギー散逸率
                U_all_rms_map(i, j) = (U_rms_map(i, j) + V_rms_map(i, j) + W_rms_map(i, j))/3.0d0
                U_all_grad_map(i, j) = (U_grad_map(i, j) + V_grad_map(i, j) + W_grad_map(i, j))/3.0d0
                epsilon_map(i, j) = nu/2.0d0*tmp_sum_map(i, j)/NZ

                ! コルモゴロフ
                L_kolmogorov_map(i, j) = (nu**3/epsilon_map(i, j))**0.25d0  ! 最小スケールの渦長さ
                T_kolmogorov_map(i, j) = (nu/epsilon_map(i, j))**0.5d0  ! 最小スケールの渦時間
                U_kolmogorov_map(i, j) = (nu*epsilon_map(i, j))**0.25d0  ! 最小スケールの渦速度

                ! テイラー長
                lambda_map(i, j) = sqrt(U_all_rms_map(i, j)/U_all_grad_map(i, j))  ! テイラー長(定義方法は色々ある)
                Re_lambda_map(i, j) = sqrt(U_all_rms_map(i, j))*lambda_map(i, j)/nu  ! テイラー長レイノルズ数
                lambda_epsilon_map(i, j) = 10*nu*U_all_rms_map(i, j)/lambda_map(i, j)**2/epsilon_map(i, j) ! テイラー長とエネルギー散逸率の関係(デバッグ用)(本当は15)

                ! 変動速度の運動エネルギー
                K_energy_map(i, j) = (U_rms_map(i, j) + V_rms_map(i, j) + W_rms_map(i, j))/2.0d0
            enddo
        enddo

        ! 弾性エネルギー
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    Ep_procs(i, j, k) = 0.5d0*(1.0d0-beta)/Re/Wi*(Lp**2.0d0-3.0d0)*log(f(C_procs(:, i, j, k)))
                enddo
            enddo
        enddo
        do j = 1, NY
            do i = 1, NX
                tmp_map(i, j) = sum(Ep_procs(i, j, :))/NZ
            enddo
        enddo
        call MPI_Allreduce(tmp_map, Ep_energy_map, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! 伸長長さ
        trC_procs(:, :, :) = C_procs(1, 1:NX, 1:NY, 1:N_procs)+C_procs(4, 1:NX, 1:NY, 1:N_procs)+C_procs(6, 1:NX, 1:NY, 1:N_procs)
        do j = 1, NY
            do i = 1, NX
                tmp_map(i, j) = sum(trC_procs(i, j, :))/NZ
            enddo
        enddo
        call MPI_Allreduce(tmp_map, trC_mean_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        do j = 1, NY
            do i = 1, NX
                tmp_map(i, j) = sum((trC_procs(i, j, :) - trC_mean_map(i, j))**2)/NZ  ! z方向の分散
            enddo
        enddo
        call MPI_Allreduce(tmp_map, trC_std_map, NX*NY, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        trC_std_map = sqrt(trC_std_map)

        ! マスクする場所
        mask_map(:, :) = 1.0d0
        do j = 1, NY
            do i = 1, NX
                if (((i-0.5d0)*dX - 2*PI  /4)**2 + ((j-0.5d0)*dX - 2*PI  /4)**2 < (DC*1.5d0/2.0d0)**2) mask_map(i, j) = 0.0d0 ! 円柱半径の1.5倍
                if (((i-0.5d0)*dX - 2*PI*3/4)**2 + ((j-0.5d0)*dX - 2*PI  /4)**2 < (DC*1.5d0/2.0d0)**2) mask_map(i, j) = 0.0d0
                if (((i-0.5d0)*dX - 2*PI  /4)**2 + ((j-0.5d0)*dX - 2*PI*3/4)**2 < (DC*1.5d0/2.0d0)**2) mask_map(i, j) = 0.0d0
                if (((i-0.5d0)*dX - 2*PI*3/4)**2 + ((j-0.5d0)*dX - 2*PI*3/4)**2 < (DC*1.5d0/2.0d0)**2) mask_map(i, j) = 0.0d0
            enddo
        enddo
        ! マスク後の値
        K_energy = sum(K_energy_map*mask_map)/sum(mask_map)
        lambda = sum(lambda_map*mask_map)/sum(mask_map)
        Re_lambda = sum(Re_lambda_map*mask_map)/sum(mask_map)
        epsilon = sum(epsilon_map*mask_map)/sum(mask_map)
        lambda_epsilon = sum(lambda_epsilon_map*mask_map)/sum(mask_map)
        L_kolmogorov = sum(L_kolmogorov_map*mask_map)/sum(mask_map)
        T_kolmogorov = sum(T_kolmogorov_map*mask_map)/sum(mask_map)
        U_kolmogorov = sum(U_kolmogorov_map*mask_map)/sum(mask_map)
        Ep_energy = sum(Ep_energy_map*mask_map)/sum(mask_map)
        trC_mean = sum(trC_mean_map*mask_map)/sum(mask_map)
        trC_std = sum(trC_std_map*mask_map)/sum(mask_map)

        ! Gstepごとにmapとして保存
        if (myrank == 0 .and. mod(step, Gstep)==0) then
            write(str, '(I8.8)') step
            open(10, file=trim(output_dir)//'Map/'//trim(key)//'_'//str//'.bin', form='unformatted', status='replace', access='stream')
            do j = 1, NY
                do i = 1, NX
                    write(10) (i-0.5d0)*dX_C, (j-0.5d0)*dY_C, K_energy_map(i, j), lambda_map(i, j), Re_lambda_map(i, j), epsilon_map(i, j), &
                              L_kolmogorov_map(i, j), T_kolmogorov_map(i, j), U_kolmogorov_map(i, j), &
                              Ep_energy_map(i, j), trC_mean_map(i, j), trC_std_map(i, j), mask_map(i, j), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
                enddo
            enddo
            close(10)
        endif
        
        ! 時系列データの保存  ! MapからK_staggeredは計算出来ない？出来てもする必要ないか
        if (myrank == 0) then
            open(30, file = trim(output_dir)//'Log/time_'//trim(key)//'_mask.d', position='append')
            write(30, *) step, step*dt_C, K_energy, 0.0d0, lambda, Re_lambda, epsilon, L_kolmogorov, T_kolmogorov, U_kolmogorov, &
                            Ep_energy, trC_mean, trC_std
            close(30)
        endif

    end subroutine all_time_dif_map


    subroutine fft_energy_dif(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, Energy, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: U_ave_procs(1:NX, 1:NY, 1:N_procs), V_ave_procs(1:NX, 1:NY, 1:N_procs), W_ave_procs(1:NX, 1:NY, 1:N_procs)
        integer, intent(in) :: step
        real(8), intent(inout) :: Energy(0:NX)
        real(8) U_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8) U_reg_procs(1:NX, 1:NY, 1:N_procs), V_reg_procs(1:NX, 1:NY, 1:N_procs), W_reg_procs(1:NX, 1:NY, 1:N_procs)
        complex(8) U_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ), V_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ), W_hat_hat_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        real(8) E_tmp_procs(1:NX/2+1, 1:NY_procs, 1:NZ), K_abs_procs(1:NX/2+1, 1:NY_procs, 1:NZ)
        real(8) Energy_procs(0:NX), Energy_sum(0:NX)
        integer i, j, k, index
        real(8) kx, ky, kz
        character(8) str
        character(32) key
        
        write(str, '(I8.8)') step
        if (sum(U_ave_procs**2) < 1.0d-10) key = 'same'  ! 平均を引かない
        if (sum(U_ave_procs**2) > 1.0d-10) key = 'diff'  ! 平均を引く

        ! 変動速度
        U_dif_procs(1:NX, 1:NY, 1:N_procs) = (U_procs(1:NX, 1:NY, 1:N_procs) - U_ave_procs(:, :, :)) * U_C
        V_dif_procs(1:NX, 1:NY, 1:N_procs) = (V_procs(1:NX, 1:NY, 1:N_procs) - V_ave_procs(:, :, :)) * U_C
        W_dif_procs(1:NX, 1:NY, 1:N_procs) = (W_procs(1:NX, 1:NY, 1:N_procs) - W_ave_procs(:, :, :)) * U_C
        call MPI_Boundary(U_dif_procs)
        call MPI_Boundary(V_dif_procs)
        call MPI_Boundary(W_dif_procs)
        call PBM(U_dif_procs)
        call PBM(V_dif_procs)
        call PBM(W_dif_procs)

        ! レギュラー格子での速度
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    U_reg_procs(i, j, k) = (U_dif_procs(i, j, k) + U_dif_procs(i-1, j, k))/2.0d0
                    V_reg_procs(i, j, k) = (V_dif_procs(i, j, k) + V_dif_procs(i, j-1, k))/2.0d0
                    W_reg_procs(i, j, k) = (W_dif_procs(i, j, k) + W_dif_procs(i, j, k-1))/2.0d0
                enddo
            enddo
        enddo

        ! 速度場をフーリエ変換  ! レギュラー格子での変動速度を用いる場合
        call fft_forward(U_reg_procs(1:NX, 1:NY, 1:N_procs), U_hat_hat_procs)
        call fft_forward(V_reg_procs(1:NX, 1:NY, 1:N_procs), V_hat_hat_procs)
        call fft_forward(W_reg_procs(1:NX, 1:NY, 1:N_procs), W_hat_hat_procs)

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
        call MPI_Allreduce(Energy_procs(0), Energy_sum(0), (NX+1), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        Energy(:) = Energy(:) + Energy_sum(:)

        if (mod(step, Gstep) == 0) then
            Energy(:) = Energy(:)/(dble(Gstep)/dble(Tstep))  ! 今までに足してきた回数で割る

            if (myrank == 0) then
                open(30, file = trim(output_dir)//'Log/energy_'//trim(key)//'_'//str//'.d', status='replace')
                do i = 0, NX-1
                    write(30, *) i, Energy(i)
                enddo
                close(30)
            endif

            Energy(:) = 0.0d0
        endif

    end subroutine fft_energy_dif


    subroutine fft_energy_dif_map(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, Energy_map_procs, step)
        real(8), intent(in) :: U_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8), intent(in) :: U_ave_procs(1:NX, 1:NY, 1:N_procs), V_ave_procs(1:NX, 1:NY, 1:N_procs), W_ave_procs(1:NX, 1:NY, 1:N_procs)
        integer, intent(in) :: step
        real(8), intent(inout) :: Energy_map_procs(1:NX, 1:NY_procs, 0:NZ)
        real(8) U_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), V_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1), W_dif_procs(0:NX+1, 0:NY+1, 0:N_procs+1)
        real(8) U_reg_procs(1:NX, 1:NY, 1:N_procs), V_reg_procs(1:NX, 1:NY, 1:N_procs), W_reg_procs(1:NX, 1:NY, 1:N_procs)
        complex(8) U_hat_procs(1:NX, 1:NY_procs, 1:NZ), V_hat_procs(1:NX, 1:NY_procs, 1:NZ), W_hat_procs(1:NX, 1:NY_procs, 1:NZ)
        real(8) E_tmp_procs(1:NX, 1:NY_procs, 1:NZ), K_abs_procs(1:NX, 1:NY_procs, 1:NZ)
        real(8) Energy_procs(1:NX, 1:NY_procs, 0:NZ)
        real(8) mask_map_procs(1:NX, 1:NY_procs)
        real(8) Energy_tmp(0:NZ), Energy(0:NZ)
        integer i, j, k, index
        ! real(8) kx, ky, kz  ! デバッグ用
        character(8) str
        character(32) key
        write(str, '(I8.8)') step
        if (sum(U_ave_procs**2) < 1.0d-10) key = 'same'  ! 平均を引かない
        if (sum(U_ave_procs**2) > 1.0d-10) key = 'diff'  ! 平均を引く

        ! 変動速度
        U_dif_procs(1:NX, 1:NY, 1:N_procs) = (U_procs(1:NX, 1:NY, 1:N_procs) - U_ave_procs(:, :, :)) * U_C
        V_dif_procs(1:NX, 1:NY, 1:N_procs) = (V_procs(1:NX, 1:NY, 1:N_procs) - V_ave_procs(:, :, :)) * U_C
        W_dif_procs(1:NX, 1:NY, 1:N_procs) = (W_procs(1:NX, 1:NY, 1:N_procs) - W_ave_procs(:, :, :)) * U_C
        call MPI_Boundary(U_dif_procs)
        call MPI_Boundary(V_dif_procs)
        call MPI_Boundary(W_dif_procs)
        call PBM(U_dif_procs)
        call PBM(V_dif_procs)
        call PBM(W_dif_procs)

        ! レギュラー格子での速度
        do k = 1, N_procs
            do j = 1, NY
                do i = 1, NX
                    U_reg_procs(i, j, k) = (U_dif_procs(i, j, k) + U_dif_procs(i-1, j, k))/2.0d0
                    V_reg_procs(i, j, k) = (V_dif_procs(i, j, k) + V_dif_procs(i, j-1, k))/2.0d0
                    W_reg_procs(i, j, k) = (W_dif_procs(i, j, k) + W_dif_procs(i, j, k-1))/2.0d0
                enddo
            enddo
        enddo

        ! 速度場をフーリエ変換  ! レギュラー格子での変動速度を用いる場合
        call fft_z(U_reg_procs(1:NX, 1:NY, 1:N_procs), U_hat_procs)
        call fft_z(V_reg_procs(1:NX, 1:NY, 1:N_procs), V_hat_procs)
        call fft_z(W_reg_procs(1:NX, 1:NY, 1:N_procs), W_hat_procs)

        ! 配列要素に対するエネルギー
        do k = 1, NZ
            do j = 1, NY_procs
                do i = 1, NX
                    E_tmp_procs(i, j, k) = abs(U_hat_procs(i, j, k))**2.0d0 &
                                         + abs(V_hat_procs(i, j, k))**2.0d0 &
                                         + abs(W_hat_procs(i, j, k))**2.0d0
                enddo
            enddo
        enddo

        ! 配列要素に対する波数
        do k = 1, NZ
            K_abs_procs(:, :, k) = sign(1.0d0, dble(NZ + 3)/2.0d0 - dble(k)) * (- dble(abs(NZ/2 + 1 - k)) + dble(NZ)/2.0d0)
        enddo

        ! 波数を四捨五入し、対応する整数の波数にエネルギーを足し合わせる
        Energy_procs(:, :, :) = 0.0d0
        do k = 1, NZ/2+1
            do j = 1, NY_procs
                do i = 1, NX
                    index = nint(K_abs_procs(i, j, k))
                    if (k==1 .or. k==NZ/2+1) then
                        Energy_procs(i, j, index) = Energy_procs(i, j, index) + E_tmp_procs(i, j, k)/2.0d0
                    else
                        Energy_procs(i, j, index) = Energy_procs(i, j, index) + E_tmp_procs(i, j, k)/2.0d0*2.0d0
                    endif
                enddo
            enddo
        enddo
        Energy_map_procs(:, :, :) = Energy_map_procs(:, :, :) + Energy_procs(:, :, :)

        ! デバッグ用
        ! ピッタリあってる!!
        ! kx = sum(U_reg_procs**2) + sum(V_reg_procs**2) + sum(W_reg_procs**2)
        ! kx = kx/(2*NX*NY*NZ)
        ! call MPI_Allreduce(kx, ky, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)  ! 変動速度を二乗して足し合わせた値
        ! kx = sum(Energy_procs)/dble(NX*NY)
        ! call MPI_Allreduce(kx, kz, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)  ! z方向のみにfftした速度場を用いて求めたエネルギーの合計
        ! if (myrank == 0) write(*, '(a, e12.4)') 'Energy Check = 1?', ky/kz


        if (mod(step, Gstep) == 0) then
            Energy_map_procs(:, :, :) = Energy_map_procs(:, :, :)/(dble(Gstep)/dble(Tstep))  ! 今までに足してきた回数で割る

            ! 空間平均したエネルギーを求める
            mask_map_procs(:, :) = 1.0d0  ! マスクする場所
            do j = 1, N_procs
                do i = 1, NX
                    if (((i-0.5d0)*dX - 2*PI  /4)**2 + ((myrank*NY_procs + j-0.5d0)*dX - 2*PI  /4)**2 < (DC*1.5d0/2.0d0)**2) mask_map_procs(i, j) = 0.0d0 ! 円柱半径の1.5倍
                    if (((i-0.5d0)*dX - 2*PI*3/4)**2 + ((myrank*NY_procs + j-0.5d0)*dX - 2*PI  /4)**2 < (DC*1.5d0/2.0d0)**2) mask_map_procs(i, j) = 0.0d0
                    if (((i-0.5d0)*dX - 2*PI  /4)**2 + ((myrank*NY_procs + j-0.5d0)*dX - 2*PI*3/4)**2 < (DC*1.5d0/2.0d0)**2) mask_map_procs(i, j) = 0.0d0
                    if (((i-0.5d0)*dX - 2*PI*3/4)**2 + ((myrank*NY_procs + j-0.5d0)*dX - 2*PI*3/4)**2 < (DC*1.5d0/2.0d0)**2) mask_map_procs(i, j) = 0.0d0
                enddo
            enddo
            do k = 0, NZ
                Energy_tmp(k) = sum(Energy_map_procs(:, :, k)*mask_map_procs(:, :))
            enddo
            call MPI_Allreduce(Energy_tmp(0), Energy(0), (NZ+1), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
            Energy(:) = Energy(:)/sum(mask_map_procs)  ! 足したx,y座標数で割る

            ! データの保存
            if (myrank == 0) then
                ! 特定のxyにおけるz方向波数のエネルギー
                open(40, file = trim(output_dir)//'Log/energy_'//trim(key)//'_z_'//str//'.d', status='replace')
                do k = 0, NZ-1
                    write(40, *) k, Energy_map_procs(1, 1, k)
                enddo
                close(40)
                ! 空間平均したエネルギー
                open(30, file = trim(output_dir)//'Log/energy_'//trim(key)//'_mask_'//str//'.d', status='replace')
                do k = 0, NZ-1
                    write(30, *) k, Energy(k)
                enddo
                close(30)
            endif

            Energy_map_procs(:, :, :) = 0.0d0
        endif

    end subroutine fft_energy_dif_map

    subroutine power_spectrum(U_step_procs, V_step_procs, W_step_procs)
        real(8), intent(in) :: U_step_procs(1:N_procs, 1:Nstep), V_step_procs(1:N_procs, 1:Nstep), W_step_procs(1:N_procs, 1:Nstep)
        complex(8) U_hat_procs(1:N_procs, 1:Nstep/2+1), V_hat_procs(1:N_procs, 1:Nstep/2+1), W_hat_procs(1:N_procs, 1:Nstep/2+1)
        ! real(8) U_debug_procs(1:N_procs, 1:Nstep), V_debug_procs(1:N_procs, 1:Nstep), W_debug_procs(1:N_procs, 1:Nstep)
        real(8) U_power_sum(1:Nstep/2+1), V_power_sum(1:Nstep/2+1), W_power_sum(1:Nstep/2+1)
        real(8) U_power_procs(1:Nstep/2+1), V_power_procs(1:Nstep/2+1), W_power_procs(1:Nstep/2+1)
        integer i, k

        do k = 1, N_procs
            call fftr2c_1d(U_step_procs(k, :)*U_C, U_hat_procs(k, :))
            call fftr2c_1d(V_step_procs(k, :)*U_C, V_hat_procs(k, :))
            call fftr2c_1d(W_step_procs(k, :)*U_C, W_hat_procs(k, :))
        enddo
        U_hat_procs(:, :) = U_hat_procs(:, :)/Nstep
        V_hat_procs(:, :) = V_hat_procs(:, :)/Nstep
        W_hat_procs(:, :) = W_hat_procs(:, :)/Nstep
        
        do i = 1, Nstep/2+1
            U_power_procs(i) = sum(abs(U_hat_procs(:, i))**2.0d0)/2.0d0
            V_power_procs(i) = sum(abs(V_hat_procs(:, i))**2.0d0)/2.0d0
            W_power_procs(i) = sum(abs(W_hat_procs(:, i))**2.0d0)/2.0d0
        enddo
        call MPI_Allreduce(U_power_procs(:), U_power_sum(:), Nstep/2+1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(V_power_procs(:), V_power_sum(:), Nstep/2+1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(W_power_procs(:), W_power_sum(:), Nstep/2+1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


        if (myrank == 0) then
            open(50, file = trim(output_dir)//'Log/step_velocity.d', status='replace')
            do i = 1, Nstep
                write(50, *) i, i*dt_C, U_step_procs(1, i), V_step_procs(1, i), W_step_procs(1, i)  ! x,yは中心、zは一番下の点の、時系列を保存
            enddo
            close(50)

            open(10, file = trim(output_dir)//'Log/Spectrum_x.d', status='replace')
            do i = 1, Nstep/2+1
                write(10, *) i-1, U_power_sum(i)  ! x,y中心のz方向に平均したuに対するパワースベクトル
            enddo
            close(10)

            open(20, file = trim(output_dir)//'Log/Spectrum_y.d', status='replace')
            do i = 1, Nstep/2+1
                write(20, *) i-1, V_power_sum(i)
            enddo
            close(20)

            open(30, file = trim(output_dir)//'Log/Spectrum_z.d', status='replace')
            do i = 1, Nstep/2+1
                write(30, *) i-1, W_power_sum(i)
            enddo
            close(30)
        endif

        
        ! デバッグ用
        ! 正しくFFTできている
        ! do k = 1, N_procs
        !     call fftc2r_1d(U_hat_procs(k, :), U_debug_procs(k, :))
        !     call fftc2r_1d(V_hat_procs(k, :), V_debug_procs(k, :))
        !     call fftc2r_1d(W_hat_procs(k, :), W_debug_procs(k, :))
        ! enddo
        ! if (myrank == 0) write(*, *) 'power_spectrum err =', sum((U_step_procs*U_C - U_debug_procs)**2)/(N_procs*Nstep)

    end subroutine power_spectrum

end module data

program main
    use smac
    use fft
    use ibm
    use data
    implicit none
    real(8), allocatable :: U_procs(:, :, :), V_procs(:, :, :), W_procs(:, :, :)
    real(8), allocatable :: P_procs(:, :, :), Phi_procs(:, :, :)
    real(8), allocatable :: Ax_procs(:, :, :), Ay_procs(:, :, :), Az_procs(:, :, :)
    real(8), allocatable :: Ax0_procs(:, :, :), Ay0_procs(:, :, :), Az0_procs(:, :, :)
    real(8), allocatable :: Bx_procs(:, :, :), By_procs(:, :, :), Bz_procs(:, :, :)
    real(8), allocatable :: Bx0_procs(:, :, :), By0_procs(:, :, :), Bz0_procs(:, :, :)
    real(8), allocatable :: Tx_procs(:, :, :), Ty_procs(:, :, :), Tz_procs(:, :, :)
    real(8), allocatable :: Tx0_procs(:, :, :), Ty0_procs(:, :, :), Tz0_procs(:, :, :)
    real(8), allocatable :: Up_procs(:, :, :), Vp_procs(:, :, :), Wp_procs(:, :, :)
    real(8), allocatable :: Fx_procs(:, :, :), Fy_procs(:, :, :), Fz_procs(:, :, :)
    real(8), allocatable :: C_procs(:, :, :, :)
    real(8), allocatable :: Cpx_procs(:, :, :, :), Cnx_procs(:, :, :, :)
    real(8), allocatable :: Cpy_procs(:, :, :, :), Cny_procs(:, :, :, :)
    real(8), allocatable :: Cpz_procs(:, :, :, :), Cnz_procs(:, :, :, :)
    real(8), allocatable :: Cx_procs(:, :, :, :)
    real(8), allocatable :: X_procs(:, :, :), Y_procs(:, :, :), Z_procs(:, :, :)
    real(8), allocatable :: Xc_procs(:, :, :), Yc_procs(:, :, :), Zc_procs(:, :, :)
    real(8), allocatable :: Uc_procs(:, :, :), Vc_procs(:, :, :), Wc_procs(:, :, :)
    real(8), allocatable :: Ua_procs(:, :, :), Va_procs(:, :, :), Wa_procs(:, :, :)
    real(8), allocatable :: Fxc_procs(:, :, :), Fyc_procs(:, :, :), Fzc_procs(:, :, :)
    real(8), allocatable :: fxint_procs(:, :, :), fyint_procs(:, :, :), fzint_procs(:, :, :)
    real(8), allocatable :: U_ave_procs(:, :, :), V_ave_procs(:, :, :), W_ave_procs(:, :, :)
    real(8), allocatable :: U_rms_procs(:, :, :), V_rms_procs(:, :, :), W_rms_procs(:, :, :)
    real(8), allocatable :: Energy_map_procs(:, :, :)
    real(8), allocatable :: U_step_procs(:, :), V_step_procs(:, :), W_step_procs(:, :)
    real(8) :: Energy(0:NX) = 0.0d0
    integer step
    real(8) time0, time1, time2, time3, time4, time5, time6, time23, time34, time45
    time23 = 0.0d0
    time34 = 0.0d0
    time45 = 0.0d0

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

    time0 = MPI_Wtime()

    call fft_init
    call init(U_procs, V_procs, W_procs, P_procs, Phi_procs, &
              Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, &
              Bx_procs, By_procs, Bz_procs, Bx0_procs, By0_procs, Bz0_procs, &
              Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, &
              Up_procs, Vp_procs, Wp_procs, Fx_procs, Fy_procs, Fz_procs, C_procs, &
              Cpx_procs, Cnx_procs, Cpy_procs, Cny_procs, Cpz_procs, Cnz_procs, Cx_procs)
    if (method == 2) then
        call ibm_init(X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, Uc_procs, Vc_procs, Wc_procs, &
                      Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
        if (input_type == 0) call ibm_vtk(Xc_procs, Yc_procs, Zc_procs)
    endif
    call ave_init(U_ave_procs, V_ave_procs, W_ave_procs, U_rms_procs, V_rms_procs, W_rms_procs, Energy_map_procs, U_step_procs, V_step_procs, W_step_procs)
    if (input_type > 0) call input_binary(U_procs, V_procs, W_procs, P_procs, C_procs, Ax0_procs, Ay0_procs, Az0_procs, Tx0_procs, Ty0_procs, Tz0_procs)
    if (input_type > 0 .and. 100*Gstep < Nstep) call vtk_binary(U_procs, V_procs, W_procs, 0)  ! 100個以上出力するなら初期値も出力する
    if (input_type > 0 .and. 100*Gstep < Nstep) call scale_vtk(U_procs, V_procs, W_procs, 0)

    time1 = MPI_Wtime()

    do step = 1, Nstep
        time2 = MPI_Wtime()
        if (beta == 1.0d0) then
            C_procs(:, :, :, :) = 0.0d0
            Tx_procs(:, :, :) = 0.0d0
            Ty_procs(:, :, :) = 0.0d0
            Tz_procs(:, :, :) = 0.0d0
        else
            call CpxCnx(C_procs, Cpx_procs, Cnx_procs)
            call CpyCny(C_procs, Cpy_procs, Cny_procs)
            call CpzCnz(C_procs, Cpz_procs, Cnz_procs)
            call Cstar(Cpx_procs, Cnx_procs, Cpy_procs, Cny_procs, Cpz_procs, Cnz_procs, U_procs, V_procs, W_procs, Cx_procs)
            call Lyapunov(Cx_procs, U_procs, V_procs, W_procs, C_procs)
            call polymer_stress(C_procs, Tx_procs, Ty_procs, Tz_procs)
        endif
        time3 = MPI_Wtime()
        
        call convection(U_procs, V_procs, W_procs, Ax_procs, Ay_procs, Az_procs)
        call viscous(U_procs, V_procs, W_procs, Bx_procs, By_procs, Bz_procs)

        if (method == 1) then
            call fft_navier(U_procs, V_procs, W_procs, P_procs, Up_procs, Vp_procs, Wp_procs, Ax_procs, Ay_procs, Az_procs, &
                            Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                            Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, Fx_procs, Fy_procs, Fz_procs, step)
            call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
            call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
        endif

        if (method == 2) then
            call ibm_preliminary(U_procs, V_procs, W_procs, P_procs, Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                                 Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, &
                                 Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, Up_procs, Vp_procs, Wp_procs, step)
            call ibm_Helmholtz(Up_procs, Vp_procs, Wp_procs, X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, &
                                Uc_procs, Vc_procs, Wc_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
            call ibm_predict(Ua_procs, Va_procs, Wa_procs, fxint_procs, fyint_procs, fzint_procs, Bx_procs, By_procs, Bz_procs, Up_procs, Vp_procs, Wp_procs)
            call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
            call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
        endif
        time4 = MPI_Wtime()
        
        U_ave_procs(:, :, :) = U_ave_procs(:, :, :) + U_procs(1:NX, 1:NY, 1:N_procs)
        V_ave_procs(:, :, :) = V_ave_procs(:, :, :) + V_procs(1:NX, 1:NY, 1:N_procs)
        W_ave_procs(:, :, :) = W_ave_procs(:, :, :) + W_procs(1:NX, 1:NY, 1:N_procs)

        ! x,y中心のz方向の格子点上の、時間に対するu,v,wを保存
        U_step_procs(:, step) = (U_procs(NX/2, NY/2, 1:N_procs) + U_procs(NX/2-1, NY/2, 1:N_procs))/2.0d0
        V_step_procs(:, step) = (V_procs(NX/2, NY/2, 1:N_procs) + V_procs(NX/2, NY/2-1, 1:N_procs))/2.0d0
        W_step_procs(:, step) = (W_procs(NX/2, NY/2, 1:N_procs) + W_procs(NX/2, NY/2, 0:N_procs-1))/2.0d0

        if (mod(step, Tstep) == 0) call all_time_dif(U_procs, V_procs, W_procs, U_rms_procs, V_rms_procs, W_rms_procs, C_procs, step)  ! U_rms_procsとかは0.0d0として代入
        if (mod(step, Tstep) == 0) call all_time_dif_map(U_procs, V_procs, W_procs, U_rms_procs, V_rms_procs, W_rms_procs, C_procs, step)
        if (mod(step, Tstep) == 0) call fft_energy_dif(U_procs, V_procs, W_procs, U_rms_procs, V_rms_procs, W_rms_procs, Energy, step)
        if (mod(step, Tstep) == 0) call fft_energy_dif_map(U_procs, V_procs, W_procs, U_rms_procs, V_rms_procs, W_rms_procs, Energy_map_procs, step)
        if (mod(step, Lstep) == 0) call log_progress(U_procs, V_procs, W_procs, C_procs, step, time1)

        if (mod(step, Gstep) == 0) then  ! 0step目に上書き保存はやめた
            call get_data_binary(U_procs, V_procs, W_procs, C_procs, step)
            call vtk_binary(U_procs, V_procs, W_procs, step)
            call scale_vtk(U_procs, V_procs, W_procs, step)
        endif
        if (mod(step, Ostep) == 0) call output_binary(U_procs, V_procs, W_procs, P_procs, C_procs, Ax_procs, Ay_procs, Az_procs, Tx_procs, Ty_procs, Tz_procs)
        time5 = MPI_Wtime()
        time23 = time23 + time3-time2
        time34 = time34 + time4-time3
        time45 = time45 + time5-time4
    enddo

    call power_spectrum(U_step_procs, V_step_procs, W_step_procs)


    if (input_type == 2) then
        U_ave_procs(:, :, :) = U_ave_procs(:, :, :)/Nstep
        V_ave_procs(:, :, :) = V_ave_procs(:, :, :)/Nstep
        W_ave_procs(:, :, :) = W_ave_procs(:, :, :)/Nstep
        call ave_mean(U_ave_procs, V_ave_procs, W_ave_procs)  ! z方向に対して平均
        ! U_ave_procs(:, :, :) = 0.0d0
        ! V_ave_procs(:, :, :) = 0.0d0
        ! W_ave_procs(:, :, :) = 0.0d0

        call input_binary(U_procs, V_procs, W_procs, P_procs, C_procs, Ax0_procs, Ay0_procs, Az0_procs, Tx0_procs, Ty0_procs, Tz0_procs)

        do step = 1, Nstep
            time2 = MPI_Wtime()
            if (beta == 1.0d0) then
                C_procs(:, :, :, :) = 0.0d0
                Tx_procs(:, :, :) = 0.0d0
                Ty_procs(:, :, :) = 0.0d0
                Tz_procs(:, :, :) = 0.0d0
            else
                call CpxCnx(C_procs, Cpx_procs, Cnx_procs)
                call CpyCny(C_procs, Cpy_procs, Cny_procs)
                call CpzCnz(C_procs, Cpz_procs, Cnz_procs)
                call Cstar(Cpx_procs, Cnx_procs, Cpy_procs, Cny_procs, Cpz_procs, Cnz_procs, U_procs, V_procs, W_procs, Cx_procs)
                call Lyapunov(Cx_procs, U_procs, V_procs, W_procs, C_procs)
                call polymer_stress(C_procs, Tx_procs, Ty_procs, Tz_procs)
            endif
            time3 = MPI_Wtime()

            call convection(U_procs, V_procs, W_procs, Ax_procs, Ay_procs, Az_procs)
            call viscous(U_procs, V_procs, W_procs, Bx_procs, By_procs, Bz_procs)

            if (method == 1) then
                call fft_navier(U_procs, V_procs, W_procs, P_procs, Up_procs, Vp_procs, Wp_procs, Ax_procs, Ay_procs, Az_procs, &
                                Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                                Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, Fx_procs, Fy_procs, Fz_procs, step)
                call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
                call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
            endif

            if (method == 2) then
                call ibm_preliminary(U_procs, V_procs, W_procs, P_procs, Ax_procs, Ay_procs, Az_procs, Ax0_procs, Ay0_procs, Az0_procs, Bx_procs, By_procs, Bz_procs, &
                                    Tx_procs, Ty_procs, Tz_procs, Tx0_procs, Ty0_procs, Tz0_procs, &
                                    Ua_procs, Va_procs, Wa_procs, Fxc_procs, Fyc_procs, Fzc_procs, Up_procs, Vp_procs, Wp_procs, step)
                call ibm_Helmholtz(Up_procs, Vp_procs, Wp_procs, X_procs, Y_procs, Z_procs, Xc_procs, Yc_procs, Zc_procs, &
                                    Uc_procs, Vc_procs, Wc_procs, Fxc_procs, Fyc_procs, Fzc_procs, fxint_procs, fyint_procs, fzint_procs)
                call ibm_predict(Ua_procs, Va_procs, Wa_procs, fxint_procs, fyint_procs, fzint_procs, Bx_procs, By_procs, Bz_procs, Up_procs, Vp_procs, Wp_procs)
                call fft_poisson(Up_procs, Vp_procs, Wp_procs, Phi_procs)
                call fft_march(Up_procs, Vp_procs, Wp_procs, U_procs, V_procs, W_procs, Phi_procs, P_procs)
            endif
            time4 = MPI_Wtime()

            U_rms_procs(:, :, :) = U_rms_procs(:, :, :) + (U_procs(1:NX, 1:NY, 1:N_procs) - U_ave_procs(:, :, :))**2
            V_rms_procs(:, :, :) = V_rms_procs(:, :, :) + (V_procs(1:NX, 1:NY, 1:N_procs) - V_ave_procs(:, :, :))**2
            W_rms_procs(:, :, :) = W_rms_procs(:, :, :) + (W_procs(1:NX, 1:NY, 1:N_procs) - W_ave_procs(:, :, :))**2

            if (mod(step, Tstep) == 0) call all_time_dif(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, C_procs, step)
            if (mod(step, Tstep) == 0) call all_time_dif_map(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, C_procs, step)
            if (mod(step, Tstep) == 0) call fft_energy_dif(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, Energy, step)
            if (mod(step, Tstep) == 0) call fft_energy_dif_map(U_procs, V_procs, W_procs, U_ave_procs, V_ave_procs, W_ave_procs, Energy_map_procs, step)
            if (mod(step, Lstep) == 0) call log_progress(U_procs, V_procs, W_procs, C_procs, step + Nstep, time1)
            time5 = MPI_Wtime()
            time23 = time23 + time3-time2
            time34 = time34 + time4-time3
            time45 = time45 + time5-time4
        enddo

        ! call ibm_var(U_rms_procs, V_rms_procs, W_rms_procs)  ! 今はいらない
        ! call output_binary(U_procs, V_procs, W_procs, P_procs, C_procs, Ax_procs, Ay_procs, Az_procs, Tx_procs, Ty_procs, Tz_procs)  ! 1周目で保存
    endif

    ! call MPI_Barrier(MPI_COMM_WORLD)
    time6 = MPI_Wtime()  ! 計測終了
    if (myrank == 0) then
        write(*, '(a9, F10.3, a3)') 'Total   :', (time6-time0), '[s]'
        write(*, '(a9, F10.3, a3, F8.3, a3)') 'init    :', time1-time0, '[s]', (time1-time0)/(time6-time0)*100, '[%]'
        write(*, '(a9, F10.3, a3, F8.3, a3)') 'polymer :', time23, '[s]', time23/(time6-time0)*100, '[%]'
        write(*, '(a9, F10.3, a3, F8.3, a3)') 'newton  :', time34, '[s]', time34/(time6-time0)*100, '[%]'
        write(*, '(a9, F10.3, a3, F8.3, a3)') 'output  :', time45, '[s]', time45/(time6-time0)*100, '[%]'
    endif

    call fft_finalize
    call MPI_Finalize(ierr)  ! 終わり
    
end program main