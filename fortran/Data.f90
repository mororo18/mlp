module rData
    implicit none

contains 

    subroutine load_matrix(c)
      real , allocatable :: c (:,:)

      integer :: io
      character(len=1500) ::line
      integer :: dimen
      integer :: i
      integer :: j
      integer :: blnk
      integer :: value_

      open(newunit=io, file="../distance_matrix", status="old", action="read") 
      read(io, "(I3)") dimen

      allocate(c(dimen, dimen))

      do i = 1, dimen
          read(io, "(A)") line
          do j = i + 1, dimen
              blnk = index(line, ' ')
              read(line(:blnk), '(I6)') value_ 
              c(i, j) = value_
              c(j, i) = value_
              line(:) = line(blnk+1:)
          end do
      end do

    end subroutine

end module rData
