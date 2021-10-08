program main

    implicit none

    integer :: io
    !character :: a
    character , dimension (15) :: a
    character(len=1500) ::line
    integer :: dimen
    integer :: i
    integer :: j
    integer :: blnk
    integer :: value_
    integer :: b
    integer :: end_
    real , allocatable :: c (:,:)

    open(newunit=io, file="../distance_matrix", status="old", action="read") 

    read(io, "(I3)") dimen
    print *, dimen

    allocate(c(dimen, dimen))

    do i = 1, dimen
        read(io, "(A)") line
        print *, line(:dimen*5)
        end_ = lnblnk(line) + 1
        do j = i + 1, dimen
            blnk = index(line, ' ')

            read(line(:blnk), '(I6)') value_ 
            print *, value_
            c(i, j) = value_
            line(:) = line(blnk+1:)
        end do
        !print *, line(b:b)
    end do
    !read(io, "(A)") line(2)
    !read(io, "(13I5)") cost(:)

    !i = 1
    !do while (i < dimen)
    !   !read(io, "(13F5.0)") c(i, i+1:dimen) 
    !   !print *, c(i, :)
    !end do

    
    close(io)

    do i = 1, dimen
    print *, c(1, i)
    end do

end program
