#include "preproc.h"
module data
    implicit none

contains 

    subroutine load_matrix(c, rnd)
      use types
      real(typeReal), allocatable, intent(out):: c (:,:)
      integer, allocatable :: rnd (:)

      integer :: io
      character(len=5500) ::line
      integer :: dimen
      integer :: i
      integer :: j
      integer :: blnk
      integer :: value_
      integer :: rnd_count
      integer :: rnd_value

      open(newunit=io, file="../distance_matrix", status="old", action="read") 
      read(io, "(I5)") dimen

      allocate(c(dimen, dimen))

      do i = 1, dimen
      !! "(A)" -> Read w=len(variable) characters as a string.
          c(i,i) = real(0.0, typeReal)
          read(io, "(A)") line
          do j = i + 1, dimen
              blnk = index(line, ' ')
              read(line(:blnk), '(I8)') value_ 
              c(i, j) = real(value_, typeReal)
              c(j, i) = real(value_, typeReal)
              line(:) = line(blnk+1:)

          end do
      end do

      read(io, *)
      read(io, *)
      read(io, '(I8)') rnd_count

      allocate(rnd(rnd_count))

      print *, "qntd rnd  ", rnd_count

      do i = 1, rnd_count
          read(io, '(I6)') rnd(i)
      end do

      close(io)

    end subroutine

end module data
