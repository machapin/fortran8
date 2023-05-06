program main
    implicit none
    character(64) :: dir='./'
    integer :: flow_type = 0
    
    if (flow_type == 0) then
       dir = trim(dir) // 'couette'
    end if
    if (flow_type == 1) then
       dir = trim(dir) // 'poiseuille'
    end if
    dir = trim(dir)//'aaaa'
    write(*, *) dir
    
end program main