module smac
    implicit none
    integer :: i, j, l, itr, istep
    integer, parameter :: Nstep = 1
    integer, parameter :: NX = 256, NY = NX
    real(8), parameter :: PI = acos(-1.0d0)
    real(8), parameter :: Xmax = 2*PI, Ymax = Xmax
    real(8), parameter :: dX = Xmax/NX, dY = Ymax/NY
    real(8), parameter :: dt = 0.02d0
    real(8), parameter :: nu = 0.01d0
    real(8), parameter :: BETA = 2.0d0 / (1 + sin(PI*dX))
    real(8), parameter :: eps = 1.0d-5
    integer, parameter :: itrmax = 1000

    real(8) :: Ax0(1:NX, 1:NY), Ay0(1:NX, 1:NY)
    real(8) :: E(NX, NY), Q(NX, NY)  ! ポアソン方程式の右辺
    real(8) :: BYM, BXM, BXP, BYP, B0
    real(8) :: er, er0

contains
    subroutine init(U, V, P, Phi)  ! はじめに実行
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1), P(0:NX+1, 0:NY+1), Phi(0:NX+1, 0:NY+1)
        write(*, *) 'BETA =', BETA
        ! 初期条件
        U(:, :) = 0.0d0
        V(:, :) = 0.0d0
        P(:, :) = 0.0d0
        Phi(:, :) = 0.0d0
    end subroutine init

    subroutine Phi_ibm(Phi)  ! 圧力の周期境界条件
        real(8), intent(inout) :: Phi(0:NX+1, 0:NY+1)
        Phi(0, :) = Phi(NX, :)
        Phi(NX+1, :) = Phi(1, :)
        Phi(:, 0) = Phi(:, NY)
        Phi(:, NY+1) = Phi(:, 1)
    end subroutine Phi_ibm

    subroutine UV_ibm(U, V)  ! 速度の周期境界条件
        real(8), intent(inout) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        U(0, :) = U(NX, :)
        U(NX+1, :) = U(1, :)
        V(0, :) = V(NX, :)
        V(NX+1, :) = V(1, :)
        U(:, 0) = U(:, NY)
        U(:, NY+1) = U(:, 1)
        V(:, 0) = V(:, NY)
        V(:, NY+1) = V(:, 1)
    end subroutine UV_ibm


    subroutine convection(U, V, Ax, Ay)  ! 発散型
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)  ! U(:, :)だと何故かエラーが発生する。見にくいけど我慢。
        real(8), intent(out) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)
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

    subroutine navier_inkai(U, V, P, Up, Vp, Ax, Ay)
        real(8), intent(in) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1), P(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(in) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)


        if (istep==1) then  ! 1ステップ目のみ例外処理
            Ax0(:, :) = Ax(:, :)
            Ay0(:, :) = Ay(:, :)
        endif

        BXM = -dt*nu/2.0d0 / dX**2
        BXP = -dt*nu/2.0d0 / dX**2
        BYM = -dt*nu/2.0d0 / dY**2
        BYP = -dt*nu/2.0d0 / dY**2
        B0 = BYM + BXM + BXP + BYP - 1.0d0

        Up(:, :) = 0.0d0  ! 一旦全消去
        Vp(:, :) = 0.0d0

        ! 速度場Upを予測
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = U(i, j) &
                         -dt*(-P(i, j)+P(i+1, j))/dX &
                         +dt*(3.0d0*Ax(i, j) - Ax0(i, j))/2.0d0 &
                         +dt*nu/2.0d0*((U(i-1, j)-2*U(i, j)+U(i+1, j))/dX**2+ (U(i, j-1)-2*U(i, j)+U(i, j+1))/dY**2)
            enddo
        enddo
        do itr = 1, itrmax
            do j = 1, NY
                do l = 1, 2  ! 奇数と偶数に分けて計算する
                    do i = l, NX, 2
                        E(i, j) = BYM*Up(i, j-1) + BXM*Up(i-1, j) - B0*Up(i, j) &
                                 +BXP*Up(i+1, j) + BYP*Up(i, j+1) - Q(i, j)
                        Up(i, j) = Up(i, j) + E(i, j) / B0
                    enddo
                enddo
            enddo
            Up(0, :) = Up(NX, :)
            Up(NX+1, :) = Up(1, :)
            Up(:, 0) = Up(:, NY)
            Up(:, NY+1) = Up(:, 1)
            
            er = sum(E**2)
            er0 = sum(Q**2)
            if (er0 == 0.0d0) then  ! Qが全て0のとき計算できない
                er0 = sum(Q**2) + 1.0d-16
            endif
            if (sqrt(er / er0) < eps) then
                write(*, *) 'Up:', 'istep =', istep, 'itr =', itr, 'er', er, 'er0', er0 
                exit
            endif
        enddo

        ! 同様にVpについても求める
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = V(i, j) &
                         -dt*(-P(i, j)+P(i, j+1))/dY &
                         +dt*(3.0d0*Ay(i, j) - Ay0(i, j))/2.0d0 &
                         +dt*nu/2.0d0*((V(i-1, j)-2*V(i, j)+V(i+1, j))/dX**2 + (V(i, j-1)-2*V(i, j)+V(i, j+1))/dY**2)
            enddo
        enddo
        do itr = 1, itrmax
            do j = 1, NY
                do l = 1, 2  ! 奇数と偶数に分けて計算する
                    do i = l, NX, 2
                        E(i, j) = BYM*Vp(i, j-1) + BXM*Vp(i-1, j) - B0*Vp(i, j) &
                                 +BXP*Vp(i+1, j) + BYP*Vp(i, j+1) - Q(i, j)
                        Vp(i, j) = Vp(i, j) + E(i, j) / B0
                    enddo
                enddo
            enddo
            Vp(0, :) = Vp(NX, :)
            Vp(NX+1, :) = Vp(1, :)
            Vp(:, 0) = Vp(:, NY)
            Vp(:, NY+1) = Vp(:, 1)

            er = sum(E**2)
            er0 = sum(Q**2)
            if (er0 == 0.0d0) then  ! Qが全て0のとき計算できない
                er0 = sum(Q**2) + 1.0d-16
            endif
            if (sqrt(er / er0) < eps) then
                write(*, *) 'Vp:', 'istep =', istep, 'itr =', itr, 'er', er, 'er0', er0 
                exit
            endif 
        enddo

        ! n-1ステップ目の保存
        Ax0(:, :) = Ax(:, :)
        Ay0(:, :) = Ay(:, :)
    end subroutine navier_inkai

    subroutine poisson(Up, Vp, Phi)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
        real(8), intent(out) :: Phi(0:NX+1, 0:NY+1)
        ! ポアソン方程式で用いる定数の計算(今回は等間隔な格子)
        BXM = 1.0d0 / dX**2
        BXP = 1.0d0 / dX**2
        BYM = 1.0d0 / dY**2
        BYP = 1.0d0 / dY**2
        B0 = BYM + BXM + BXP + BYP
        ! 右辺の計算
        do j = 1, NY
            do i = 1, NX
                Q(i, j) = ((-Up(i-1, j) + Up(i, j)) / dX &
                          +(-Vp(i, j-1) + Vp(i, j)) / dY) / dt
            enddo
        enddo
        Phi(:, :) = 0.0d0
        ! SOR法
        do itr = 1, itrmax
            do j = 1, NY
                do l = 1, 2  ! 奇数と偶数に分けて計算する
                    do i = l, NX, 2
                        E(i, j) = BYM*Phi(i, j-1) + BXM*Phi(i-1, j) - B0*Phi(i, j) &
                                 +BXP*Phi(i+1, j) + BYP*Phi(i, j+1) - Q(i, j)
                        Phi(i, j) = Phi(i, j) + BETA * E(i, j) / B0
                    enddo
                enddo
            enddo
            Phi(0, :) = Phi(NX, :)  ! これいるくね？9/20
            Phi(NX+1, :) = Phi(1, :)
            Phi(:, 0) = Phi(:, NY)
            Phi(:, NY+1) = Phi(:, 1)
            
            er = sum(E**2)
            er0 = sum(Q**2)
            if (er0 == 0.0d0) then  ! Qが全て0のとき計算できない
                er0 = sum(Q**2) + 1.0d-16
            endif
            if (sqrt(er / er0) < eps) then
                write(*, *) 'Ph:', 'istep =', istep, 'itr =', itr, 'er', er, 'er0', er0 
                exit
            endif
            ! if (mod(itr, 10000)==0) then
            !     write(*, *) 'istep =', istep, 'itr =', itr, 'err =', sqrt(er / er0)
            ! endif
        enddo
    end subroutine poisson

    subroutine march_inkai(Up, Vp, U, V, Phi, P)
        real(8), intent(in) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1), Phi(0:NX+1, 0:NY+1)
        real(8), intent(out) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
        real(8), intent(inout) :: P(0:NX+1, 0:NY+1)
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
                P(i, j) = P(i, j) + Phi(i, j) - nu/2.0d0*((Phi(i-1, j)-2*Phi(i, j)+Phi(i+1, j))/dX**2 &
                                                         +(Phi(i, j-1)-2*Phi(i, j)+Phi(i, j+1))/dY**2)
            enddo
        enddo
        P(0, :) = P(NX, :)
        P(NX+1, :) = P(1, :)
        P(:, 0) = P(:, NY)
        P(:, NY+1) = P(:, 1)
    end subroutine march_inkai

end module smac

program main
    use smac
    implicit none
    real(8) :: U(0:NX+1, 0:NY+1), V(0:NX+1, 0:NY+1)
    real(8) :: Ax(1:NX, 1:NY), Ay(1:NX, 1:NY)  ! 対流項の計算
    real(8) :: Up(0:NX+1, 0:NY+1), Vp(0:NX+1, 0:NY+1)
    real(8) :: P(0:NX+1, 0:NY+1)
    real(8) :: Phi(0:NX+1, 0:NY+1)

    call init(U, V, P, Phi)
    ! 境界条件（今回は全て周期境界）
    call Phi_ibm(Phi)  ! 周期境界
    call UV_ibm(U, V)

    do istep = 1, Nstep
        ! 対流項、粘性項の計算
        call convection(U, V, Ax, Ay)
        ! call viscous(U, V, Bx, By)
        call navier_inkai(U, V, P, Up, Vp, Ax, Ay)

        call poisson(Up, Vp, Phi)
        ! call Phi_ibm(Phi) ! 圧力の周期境界条件

        call march_inkai(Up, Vp, U, V, Phi, P)
        call UV_ibm(U, V)  ! 速度の周期境界条件
    enddo
    
end program main