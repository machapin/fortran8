! mpiifort sample4.f90
! qsub jobscript.sh

program main
    implicit none
    include 'mpif.h'
    integer, parameter :: NX = 4, NY = NX, NZ = NX
    real(8), allocatable :: U_procs(:, :, :)
    character(2) chmyrank
    integer N_procs
    integer ierr, procs, myrank
    integer seed
    integer seeds(1)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

    N_procs = NZ / procs
    allocate(U_procs(0:NX+1, 0:NY+1, 0:N_procs+1))

    call system_clock(seed)  ! 現在時刻をシード値として使用する
    seed = seed + myrank  ! プロセスごとにシード値を変える
    seeds(1) = seed
    call random_seed(put=seeds)  ! 入力するシード値は配列
    
    call random_number(U_procs)

    write(chmyrank, '(I2.2)') myrank
    open(10, file = './data/rank'//chmyrank//'.d')
    write(10, '(100e12.4)') U_procs(1:NX, 1:NY, 1:N_procs)
    close(10)


    call MPI_Finalize(ierr)
    
end program main
