module data_types
    implicit none

    type tInfo
        integer :: dimen
        real, allocatable :: cost (:,:)
        integer :: T = 1
        integer :: C = 2
        integer :: W = 3
        
        integer :: REINSERTION  = 1
        integer :: OR_OPT_2     = 2
        integer :: OR_OPT_3     = 3 
        integer :: SWAP         = 4
        integer :: TWO_OPT      = 5
    end type

    type tSolution
        integer, allocatable :: s (:)
        real, allocatable :: seq (:,:,:)
        real :: cost
    end type

end module

subroutine print_matrix(c)
    implicit none
    real, allocatable :: c(:,:)

    integer, allocatable :: nxm(:)
    integer :: dimen
    integer :: i
    integer :: j

    nxm = shape(c)
    dimen = nxm(1)

    do i=1, dimen
        do j=1, dimen 
            write(*, "(F6.1, a)", advance="no") c(i,j),  ''
        end do
        print *, new_line('A')
    end do

end subroutine

subroutine print_info(info)
    use data_types
    type(tInfo) :: info

    interface
        subroutine print_matrix(c)
            real, allocatable :: c(:,:)
        end subroutine
    end interface

    print *, info%dimen
    call print_matrix(info%cost)

end subroutine

subroutine subseq_load(sol, info)
    use data_types
    type(tSolution) :: sol
    type(tInfo) :: info
    integer :: i
    integer :: j
    integer :: k
    integer :: j_prev

    do i = 1, info%dimen+1
        k = 1 - i - merge(0, 1, i /= 1)

        sol%seq(i, i, info%T) = 0.0
        sol%seq(i, i, info%C) = 0.0
        sol%seq(i, i, info%W) = merge(1.0, 0.0, i /= 1)
        do j = i+1, info%dimen+1
            j_prev = j-1

            sol%seq(i,j,info%T) = info%cost(sol%s(j_prev), sol%s(j)) + sol%seq(i,j_prev,info%T)
            sol%seq(i,j,info%C) = sol%seq(i,j,info%T) + sol%seq(i,j_prev,info%C)
            sol%seq(i,j,info%W) = j+k

        end do
    end do
end subroutine

subroutine sort(cL, lmt, r, info)
    use data_types
    implicit none

    type(tInfo) :: info
    integer, dimension(info%dimen-1), intent(out) :: cL
    integer :: lmt
    integer :: r

    integer :: i
    integer :: j
    integer :: tmp

    do i=1, lmt
        do j=1, lmt-i
            if (info%cost(r, cL(j)) > info%cost(r, cL(j+1))) then
                tmp = cL(j)
                cL(j) = cL(j+1)
                cL(j+1) = tmp
            end if
        end do
    end do

end subroutine

subroutine arr_shift(arr, src, tgt, sz)
    implicit none 
    integer, dimension(*), intent(out) :: arr
    integer :: src, tgt
    integer :: i
    integer :: sz
    integer :: diff

    if (src < tgt) then
        diff = tgt - src
        do i=src+sz-1, src, -1
            arr(i+diff) = arr(i)
        end do
    else 
        do i=1, sz
            arr(tgt+i-1) = arr(i+src-1)
        end do
    end if

end subroutine

function construction(alpha, info) result(ret)
    use data_types
    implicit none

    real :: alpha
    type(tInfo) :: info
    integer, dimension(info%dimen+1) :: ret 

    integer :: i
    integer :: j
    integer, dimension(info%dimen+1) :: s
    integer, dimension(info%dimen-1) :: cL
    integer :: cL_size
    integer :: cN
    integer :: r
    integer :: rng
    integer :: index_
    real :: RND

    cL_size = info%dimen-1

    ! init cL
    cL = (/ (I, I=2, info%dimen) /)
    print *, cl
   !do i=2, info%dimen
   !    cL(i-1) = i
   !end do

    s(1) = 1
    r = 1
    do j=2, info%dimen
        call sort(cL, cL_size, r, info)

        rng = ceiling(cL_size * alpha)
        call random_number(RND)
        RND = merge(RND+0.0000000001, RND, RND < 0.0000000001)
        index_ = ceiling(rng * RND)

        cN = cL(index_)

        ! memmove on cL
        call arr_shift(cL, index_+1, index_, cL_size-index_)
        cL_size = cL_size-1

        s(j) = cN
        r = cN
    end do

    s(info%dimen+1) = 1

    ret = s(:)

end function

subroutine search_swap(solut, info, ret) 
    use data_types
    implicit none
    type(tSolution) :: solut
    type(tInfo) :: info
    logical, intent(out) :: ret
end subroutine

subroutine search_two_opt(solut, info, ret) 
    use data_types
    implicit none
    type(tSolution) :: solut
    type(tInfo) :: info
    logical, intent(out) :: ret
end subroutine

subroutine search_reinsertion(solut, info, opt, ret) 
    use data_types
    implicit none
    type(tSolution) :: solut
    type(tInfo) :: info
    integer :: opt
    logical, intent(out) :: ret
end subroutine

subroutine RVND(sol, info)
    use data_types
    implicit none

    type(tSolution) :: sol
    type(tInfo) :: info

    real :: rnd

    integer, dimension(5) :: neighbd_list
    integer :: nl_size
    integer :: neighbd

    integer :: index_
    logical :: improve_flag

    do while (nl_size > 0)
        call random_number(rnd)
        rnd = merge(rnd+0.0000000001, rnd, rnd < 0.0000000001)
        index_ = ceiling(rnd*nl_size)
        neighbd = neighbd_list(index_)

        improve_flag = .false.

        if (neighbd == info%REINSERTION) then
            call search_reinsertion(sol, info, info%REINSERTION, improve_flag)
        else if (neighbd == info%OR_OPT_2) then
            call search_reinsertion(sol, info, info%OR_OPT_2, improve_flag)
        else if (neighbd == info%OR_OPT_3) then
            call search_reinsertion(sol, info, info%OR_OPT_3, improve_flag)
        else if (neighbd == info%SWAP) then
            call search_swap(sol, info, improve_flag)
        else if (neighbd == info%TWO_OPT) then
            call search_two_opt(sol, info, improve_flag)
        endif

        if (improve_flag) then
            neighbd_list(:) = (/ info%REINSERTION, info%OR_OPT_2, info%OR_OPT_3, info%SWAP, info%TWO_OPT /)
        else
            call arr_shift(neighbd_list, index_+1, index_, nl_size-index_)
            nl_size = nl_size-1
        endif
        
    end do
end subroutine

function GILS_RVND(Imax, Iils, R, info) result(ret)
    use data_types

    implicit none

    !! parameters
    integer :: Imax
    integer :: Iils
    real, dimension(26) :: R
    type(tInfo) :: info
    type(tSolution) :: ret

    !! solutions
    type(tSolution) :: sol_best
    type(tSolution) :: sol_partial
    type(tSolution) :: sol_crnt

    !! aux variables
    integer :: i
    integer :: index_
    integer :: R_size = 26
    real :: alpha
    real :: rnd
    integer :: iterILS

    interface
        function construction(alpha, info) result (ret)
            use data_types
            implicit none
            real :: alpha
            type(tInfo) :: info
            integer, dimension(info%dimen+1) :: ret 
        end function
    end interface

    do i=1, Imax
        call random_number(rnd)
        rnd = merge(rnd+0.0000000001, rnd, rnd < 0.0000000001)
        index_ = ceiling(rnd*R_size)
        alpha = R(index_)
        sol_crnt%s = construction(alpha, info)
        sol_partial = sol_crnt

        iterILS = 0
        do while (iterILS < Iils)
            call RVND(sol_crnt, info)

            if (sol_crnt%cost < sol_partial%cost - EPSILON(1.0)) then
                sol_partial = sol_crnt
                sol_partial%cost = sol_partial%cost - EPSILON(1.0)
                iterILS = 0
            endif

            call subseq_load(sol_crnt, info)

            iterILS = iterILS + 1
        end do
        
        call subseq_load(sol_partial, info)

        if (sol_partial%cost < sol_best%cost) then
            sol_best = sol_partial
        endif
    
    end do

end function

program main
    use rData
    use data_types

    implicit none

    real , allocatable :: c (:,:)
    integer , allocatable :: s (:)
    integer, allocatable :: dimensions (:)
    integer :: dimen
    type(tInfo) :: info
    type(tSolution) :: sol
    integer :: i

    interface
        function construction(alpha, info) result (ret)
            use data_types
            implicit none
            real :: alpha
            type(tInfo) :: info
            integer, dimension(info%dimen+1) :: ret 
        end function
    end interface

    call load_matrix(info%cost)

    !print *, info% c
    dimensions = shape(info%cost(:,:))
    info%dimen = dimensions(1)

    call print_info(info)

    allocate(sol%s(info%dimen+1))
    allocate(sol%seq(info%dimen+1, info%dimen+1, 3))
    do i=1, info%dimen
        sol%s(i) = i
    end do
    sol%s(info%dimen+1) = 1

    !print *, sol%s(:)

    call subseq_load(sol, info)

    print *, sol%seq(1, info%dimen+1, info%C)

    sol%s = construction(0.2, info)

    print *, sol%s

end program