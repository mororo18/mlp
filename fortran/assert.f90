module assertion
      implicit none
      interface assert
          procedure assert_dbl, assert_flt,  assert_int, assert_lgc 
      end interface
  contains

    subroutine error(filename, line)
      integer :: line
      character(len=*) :: filename
      character(len=128) :: buff

      write(buff, "(I4)") line
      print *, "Assertion failed. File " // filename // ":" // trim(buff)  
      call exit(0)
    end subroutine

    subroutine assert_dbl(a, b, filename, line)
      real(8) :: a, b
      integer :: line
      character(len=*) :: filename
      character(len=128) :: buff

      if (ABS(a-b) .gt. EPSILON(1.0)) then
          call error(filename, line)
      endif

  end subroutine


    subroutine assert_flt(a, b, filename, line)
      real(4) :: a, b
      integer :: line
      character(len=*) :: filename
      character(len=128) :: buff

      if (ABS(a-b) .gt. EPSILON(1.0)) then
          call error(filename, line)
      endif

      end subroutine

    subroutine assert_int(a, b, filename, line)
      integer :: a, b
      integer :: line
      character(len=*) :: filename
      character(len=128) :: buff

      if (a .ne. b) then
          call error(filename, line)
      endif

  end subroutine


    subroutine assert_lgc(a, b, filename, line)
      logical :: a, b
      integer :: line
      character(len=*) :: filename
      character(len=128) :: buff

      if (a .neqv. b) then
          call error(filename, line)
      endif

  end subroutine

end module
