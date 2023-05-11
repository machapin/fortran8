program main
    implicit none
    ! integer, parameter :: n = 100
    integer, parameter :: nx = 100, ny = 50 ! 配列のサイズ
    real(8) :: data1(nx, ny), data2(nx, ny) ! 2次元配列
    integer :: i, j
    character(len=20) :: filename1
    character(len=20) :: filename2

    do i = 1, nx
        do j = 1, ny
            data1(i,j) = dble(i*j) ! i*jを倍精度浮動小数点数に変換して代入
            data2(i,j) = dble(i*j)*10.0d0 ! i*jを倍精度浮動小数点数に変換して代入
        end do
    end do

    filename1 = "data/data.bin"
    open(10, file=filename1, access='stream', form='unformatted')
    do j = 1, ny
        do i = 1, nx
            write(10) data1(i, j), data2(i, j)
        end do
    end do
    close(10)

    filename2 = "data/data.d"
    open(20, file=filename2)
    do j = 1, ny
        do i = 1, nx
            write(20, *) data1(i, j), data2(i, j)
        end do
    end do
    close(20)
  
    
end program main