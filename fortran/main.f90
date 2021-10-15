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

        real :: fmax = 3.4028235E+38
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
    !print *, cl
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

subroutine swap(s, i, j)
    implicit none

    integer, dimension(*), intent(out) :: s 
    integer :: i
    integer :: j

    integer :: tmp

    tmp = s(i)
    s(i) = s(j)
    s(j) = tmp

end subroutine

subroutine reverse(s, i, j)
    implicit none

    integer, dimension(*), intent(out) :: s 
    integer :: i
    integer :: j

    integer :: index_

    do index_=i, int((i+j)/2)
        call swap(s, i, j)
    end do

end subroutine

subroutine reinsert(s, i, j, pos)
    implicit none

    integer, dimension(*), intent(out) :: s 
    integer :: i
    integer :: j
    integer :: pos

    integer, allocatable :: sub (:)

    integer :: sz
    allocate(sub(j-i+1))

    sub(:) = s(i:j)
     
    if (i < pos) then
        sz = pos-j-1
        call arr_shift(s, j+1, i, sz)
        s(i+sz: pos-1) = sub(:)
    else
        sz = i-pos
        call arr_shift(s, pos, j+1-sz, sz)
        s(pos:pos+j-i) = sub(:)
    endif

end subroutine

subroutine search_swap(solut, info, ret) 
    use data_types
    implicit none
    type(tSolution) :: solut
    type(tInfo) :: info
    logical, intent(out) :: ret

    integer :: i
    integer :: j
    integer :: i_prev
    integer :: i_next
    integer :: j_prev
    integer :: j_next

    integer :: I_best
    integer :: J_best

    real :: cost_best 
    real :: cost_new
    real :: cost_concat_1
    real :: cost_concat_2 
    real :: cost_concat_3 
    real :: cost_concat_4

    cost_best = info%fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0
    cost_concat_4 = 0.0

    do i=2, info%dimen-1
        i_prev = i-1
        i_next = i+1

        cost_concat_1 =                 solut%seq(1, i_prev, info%T) + info%cost(solut%s(i_prev), solut%s(i_next))
        cost_concat_2 = cost_concat_1 + solut%seq(i, i_next, info%T) + info%cost(solut%s(i), solut%s(i_next+1))

        cost_new = solut%seq(1, i_prev, info%C)                                                    + &
                solut%seq(i, i_next, info%W)               * (cost_concat_1) + info%cost(solut%s(i_next), solut%s(i))  + &
                solut%seq(i_next+1, info%dimen+1, info%W)   * (cost_concat_2) + solut%seq(i_next+1, info%dimen+1, info%C) 

        if (cost_new < cost_best) then
            cost_best = cost_new - EPSILON(1.0)
            I_best = i
            J_best = i_next
        endif

        do j=i_next+1, info%dimen

            j_prev = j-1
            j_next = j+1


            cost_concat_1 =                 solut%seq(1, i_prev, info%T)       + info%cost(solut%s(i_prev), solut%s(j))
            cost_concat_2 = cost_concat_1                           + info%cost(solut%s(j), solut%s(i_next))
            cost_concat_3 = cost_concat_2 + solut%seq(i_next, j_prev, info%T)  + info%cost(solut%s(j_prev), solut%s(i))
            cost_concat_4 = cost_concat_3                           + info%cost(solut%s(i), solut%s(j_next))


            cost_new = solut%seq(1, i_prev, info%C)                                                 + &     ! 1st subseq
                    cost_concat_1 + &                                                           ! concat 2nd subseq (single node)
                    solut%seq(i_next, j_prev, info%W)      * cost_concat_2 + solut%seq(i_next, j_prev, info%C) + &    ! concat 3rd subseq
                    cost_concat_3 + &                                                           ! concat 4th subseq (single node)
                    solut%seq(j_next, info%dimen+1, info%W) * cost_concat_4 + solut%seq(j_next, info%dimen+1, info%C)   ! concat 5th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new - EPSILON(1.0);
                I_best = i;
                J_best = j;
            endif

        end do
        
    end do

    if (cost_best < solut%seq(1, info%dimen+1, info%C) - EPSILON(1.0)) then
        call swap(solut%s, I_best, J_best)
        call subseq_load(solut, info)
        ret = .true.
    endif

    ret = .false.

end subroutine

subroutine search_two_opt(solut, info, ret) 
    use data_types
    implicit none
    type(tSolution) :: solut
    type(tInfo) :: info
    logical, intent(out) :: ret

    integer :: i
    integer :: j
    integer :: i_prev
    integer :: j_next

    integer :: I_best
    integer :: J_best

    real :: cost_best 
    real :: cost_new
    real :: cost_concat_1
    real :: cost_concat_2 
    real :: cost_concat_3 
    real :: cost_concat_4
    real :: rev_seq_cost

    cost_best = info%fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0

    do i=2, info%dimen-1
        i_prev = i-1
        rev_seq_cost = solut%seq(i, i+1, info%T)

        do j=i+2, info%dimen
            j_next = j+1

            rev_seq_cost = rev_seq_cost + info%cost(solut%s(j-1), solut%s(j)) * (solut%seq(i, j, info%W)-1.0)

            cost_concat_1 =                 solut%seq(1, i_prev, info%T)   + info%cost(solut%s(j), solut%s(i_prev))
            cost_concat_2 = cost_concat_1 + solut%seq(i, j, info%T)        + info%cost(solut%s(j_next), solut%s(i))

            cost_new = solut%seq(1, i_prev, info%C)                                                        + & !   1st subseq
                    solut%seq(i, j, info%W)                * cost_concat_1 + rev_seq_cost                  + & ! concat 2nd subseq (reversed seq)
                    solut%seq(j_next, info%dimen+1, info%W) * cost_concat_2 + solut%seq(j_next, info%dimen+1, info%C)       ! concat 3rd subseq

            if (cost_new < cost_best) then
                cost_best = cost_new - EPSILON(1.0)
                I_best = i
                J_best = j
            endif

        end do

    end do

    if (cost_best < solut%seq(1, info%dimen+1, info%C)) then
        call reverse(solut%s, I_best, J_best)
        call subseq_load(solut, info)
        ret = .true.
    endif

    ret = .false.

end subroutine

subroutine search_reinsertion(solut, info, opt, ret) 
    use data_types
    implicit none
    type(tSolution) :: solut
    type(tInfo) :: info
    integer :: opt
    logical, intent(out) :: ret

    integer :: i
    integer :: j
    integer :: k
    integer :: i_prev
    integer :: j_next
    integer :: k_next

    integer :: I_best
    integer :: J_best
    integer :: POS_best

    real :: cost_best 
    real :: cost_new
    real :: cost_concat_1
    real :: cost_concat_2 
    real :: cost_concat_3 
    real :: cost_concat_4

    cost_best = info%fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0

    do i=2, info%dimen-opt+1
        j = opt+i-1
        i_prev = i-1
        j_next = j+1

        do k=1, i_prev-1

            cost_concat_1 =                 solut%seq(1, k, info%T)            + info%cost(solut%s(k), solut%s(i))
            cost_concat_2 = cost_concat_1 + solut%seq(i, j, info%T)            + info%cost(solut%s(j), solut%s(k_next))
            cost_concat_3 = cost_concat_2 + solut%seq(k_next, i_prev, info%T)  + info%cost(solut%s(i_prev), solut%s(j_next))

            cost_new = solut%seq(1, k, info%C)                                                             + & !       1st subseq
                    solut%seq(i, j, info%W)                * cost_concat_1 + solut%seq(i, j, info%C)                  + & ! concat 2nd subseq (reinserted seq)
                    solut%seq(k_next, i_prev, info%W)      * cost_concat_2 + solut%seq(k_next, i_prev, info%C)        + & ! concat 3rd subseq
                    solut%seq(j_next, info%dimen+1, info%W) * cost_concat_3 + solut%seq(j_next, info%dimen+1, info%C)       ! concat 4th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new - EPSILON(1.0)
                I_best = i
                J_best = j
                POS_best = k
            endif

        end do

        do k=i+opt, info%dimen-opt-1
            k_next = k+1

            cost_concat_1 =                 solut%seq(1, i_prev, info%T)   + info%cost(solut%s(i_prev), solut%s(j_next))
            cost_concat_2 = cost_concat_1 + solut%seq(j_next, k, info%T)   + info%cost(solut%s(k), solut%s(i))
            cost_concat_3 = cost_concat_2 + solut%seq(i, j, info%T)        + info%cost(solut%s(j), solut%s(k_next))

            cost_new = solut%seq(1, i_prev, info%C)                                                        + & !       1st subseq
                    solut%seq(j_next, k, info%W)           * cost_concat_1 + solut%seq(j_next, k, info%C)             + & ! concat 2nd subseq
                    solut%seq(i, j, info%W)                * cost_concat_2 + solut%seq(i, j, info%C)                  + & ! concat 3rd subseq (reinserted seq)
                    solut%seq(k_next, info%dimen+1, info%W) * cost_concat_3 + solut%seq(k_next, info%dimen+1, info%C)       ! concat 4th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new - EPSILON(1.0)
                I_best = i
                J_best = j
                POS_best = k
            endif
        end do
    end do

    if (cost_best < solut%seq(1, info%dimen+1, info%C)) then
        call reinsert(solut%s, I_best, J_best, POS_best+1)
        call subseq_load(solut, info)
        ret = .true.
    endif

    ret = .false.

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

    call reinsert(sol%s, 8, 12, 2)
    print *, sol%s

end program
