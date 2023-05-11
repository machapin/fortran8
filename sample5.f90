program main
    implicit none
    integer, parameter :: N = 3
    integer :: info
    real(8) :: A(N, N), w(N), work(3*N-1)  ! 3*N-1はこれ以上大きい値なら大丈夫
    integer :: lwork = 3*N-1
    
    ! Aを適当に初期化する
    A = reshape([2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0], [N, N])
    
    ! DSYEVを呼び出してAの固有値と固有ベクトルを計算する
    ! https://wwwnucl.ph.tsukuba.ac.jp/~hinohara/compphys2-18/doc/compphys2-8.pdf
    call dsyev('N', 'U', N, A, N, w, work, lwork, info)
    
    if (info /= 0) then
      print *, "DSYEV failed with info = ", info
    else
    ! 結果を表示する
    print *, "Eigenvalues:"
    print *, w
    ! print *, "Eigenvectors:"  ! 'N'なら固有ベクトルは計算されない
    ! print *, A
    end if
    
end program main

! ifort sample5.f90 -mkl
! gfortran sample5.f90 -llapack

! frtpx sample5.f90 -SSL2
! jxsub single.sh