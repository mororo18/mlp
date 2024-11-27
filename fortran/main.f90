#include "preproc.h"

subroutine print_matrix(c)
    use types
    implicit none
    real(typeReal), allocatable :: c(:,:)

    integer, allocatable :: nxm(:)
    integer :: dimen
    integer :: i
    integer :: j

    nxm = shape(c)
    dimen = nxm(1)

    do i=1, 2
        do j=1, dimen 
            write(*, "(F8.3, a)", advance="no") c(i,j),  ' '
        end do
        print *, new_line('A')
        print *, ''
    end do

end subroutine

subroutine print_info(info)
    use types
    type(tInfo) :: info

    interface
        subroutine print_matrix(c)
            real(typeReal), allocatable :: c(:,:)
        end subroutine
    end interface

    print *, info%dimen
    call print_matrix(info%cost)

end subroutine

subroutine subseq_load(sol, info)
    use types
    type(tSolution) :: sol
    type(tInfo) :: info
    integer :: i
    integer :: j
    integer :: k
    integer :: j_prev

    do i = 1, info%dimen+1
        k = 1 - i - merge(0, 1, i .ne. 1)

        sol%seq( i, i)%T = 0.0
        sol%seq( i, i)%C = 0.0
        sol%seq( i, i)%W = merge(1.0, 0.0, i .ne. 1)
        do j = i+1, info%dimen+1

            j_prev = j-1

            sol%seq( j, i)%T = info%cost(sol%s(j_prev), sol%s(j)) + sol%seq( j_prev,i)%T
            sol%seq(j,i)%C = sol%seq( j, i)%T + sol%seq(j_prev,i)%C
            sol%seq(j,i)%W = j+k

        end do
    end do

    sol%cost = sol%seq( info%dimen+1, 1)%C

end subroutine

subroutine sort(arr, len, r, info)
    use types
    implicit none 

    integer, intent(inout) :: arr(len)
    integer :: len
    integer :: r
    type(tInfo), intent(in) :: info
    call quicksort(arr, 1, len, info, r)
end subroutine sort


function partition(arr, left, right, info, r) result(ret)
    use types
    implicit none 

    integer, intent(inout) :: arr(*)
    integer, intent(in) :: left, right, r
    type(tInfo), intent(in) :: info
    integer :: pivot, i, j, temp
    integer :: ret

    pivot = arr(right)
    i = left - 1

    do j = left, right - 1
        if (info%cost(r, arr(j)) < info%cost(r, pivot)) then
            i = i + 1
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
        end if
    end do
    temp = arr(i + 1)
    arr(i + 1) = arr(right)
    arr(right) = temp

    ret = i + 1
end function partition

recursive subroutine quicksort(arr, left, right, info, r)
    use types
    implicit none 

    integer, intent(inout) :: arr(*)
    integer, intent(in) :: left, right, r
    type(tInfo), intent(in) :: info
    integer :: pivot

    interface
        function partition(arr, left, right, info, r) result(ret)
            use types
            implicit none 

            integer, intent(inout) :: arr(*)
            integer, intent(in) :: left, right, r
            type(tInfo), intent(in) :: info
            integer :: ret
        end function
    end interface

    if (left < right) then
        pivot = partition(arr, left, right, info, r)
        call quicksort(arr, left, pivot - 1, info, r)
        call quicksort(arr, pivot + 1, right, info, r)
    end if
end subroutine quicksort

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
    use types
    implicit none

    real(typeReal) :: alpha
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
    real(typeReal) :: RND

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

        index_ = info%rnd(info%rnd_index) + 1
        info%rnd_index = info%rnd_index + 1

        cN = cL(index_)

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
    integer :: bound

    bound = j
    do index_=i, int((i+j)/2)
        call swap(s, index_, bound)
        bound = bound - 1
    end do

end subroutine

subroutine reinsert(s, i, j, pos)
    implicit none

    integer, dimension(*), intent(out) :: s 
    integer :: i
    integer :: j
    integer :: pos

    integer, dimension(j-i+1) :: sub
    integer :: sz

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

subroutine is_feasible(solut, ret)
    use types
    implicit none
    type(tSolution) :: solut
    logical, intent(out) :: ret
    integer :: i
    integer :: n
    logical, dimension(solut%s_size-1) :: check
    n = solut%s_size-1
    check = (/ (.false.,i=1,n) /)

    do i=1, n
        check(solut%s(i)) = .true.
    end do

    ret = .true.
    do i=1, n
        if (check(i) .eqv. .false.) then
            ret = .false.
        endif
    end do

end subroutine

subroutine search_swap(solut, info, ret) 
#ifdef DEBUG
    use assertion
#endif
    use types
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

    real(typeReal) :: cost_best 
    real(typeReal) :: cost_new
    real(typeReal) :: cost_concat_1
    real(typeReal) :: cost_concat_2 
    real(typeReal) :: cost_concat_3 
    real(typeReal) :: cost_concat_4

    logical :: fsb

    cost_best = info%fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0
    cost_concat_4 = 0.0

    do i=2, info%dimen-1
        i_prev = i-1
        i_next = i+1

        cost_concat_1 =                 solut%seq( i_prev,1)%T + info%cost(solut%s(i_prev), solut%s(i_next))
        cost_concat_2 = cost_concat_1 + solut%seq( i_next,i)%T + info%cost(solut%s(i), solut%s(i_next+1))

        cost_new = solut%seq( i_prev,1)%C                                                                       + &
                solut%seq( i_next,i)%W               * (cost_concat_1) + info%cost(solut%s(i_next), solut%s(i)) + &
                solut%seq( info%dimen+1,i_next+1)%W  * (cost_concat_2) + solut%seq( info%dimen+1,i_next+1)%C 

        if (cost_new < cost_best) then
            cost_best = cost_new 
            I_best = i
            J_best = i_next
        endif

        do j=i_next+1, info%dimen

            j_prev = j-1
            j_next = j+1


            cost_concat_1 =                 solut%seq( i_prev,1)%T       + info%cost(solut%s(i_prev), solut%s(j))
            cost_concat_2 = cost_concat_1                           + info%cost(solut%s(j), solut%s(i_next))
            cost_concat_3 = cost_concat_2 + solut%seq( j_prev,i_next)%T  + info%cost(solut%s(j_prev), solut%s(i))
            cost_concat_4 = cost_concat_3                           + info%cost(solut%s(i), solut%s(j_next))


            cost_new = solut%seq( i_prev,1)%C                                                 + &     ! 1st subseq
                    cost_concat_1 + &                                                           ! concat 2nd subseq (single node)
                    solut%seq( j_prev,i_next)%W       * cost_concat_2 + solut%seq( j_prev,i_next)%C + &    ! concat 3rd subseq
                    cost_concat_3 + &                                                           ! concat 4th subseq (single node)
                    solut%seq( info%dimen+1,j_next)%W * cost_concat_4 + solut%seq( info%dimen+1,j_next)%C   ! concat 5th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new 
                I_best = i
                J_best = j
            endif

        end do
        
    end do

    if (cost_best < solut%cost - EPSILON(1.0)) then
        !print *, solut%cost, EPSILON(1.0)
        call swap(solut%s, I_best, J_best)
        call subseq_load(solut, info)
#ifdef DEBUG
        fsb = .false.
        call is_feasible(solut, fsb)
        !print *, cost_best, solut%cost
        ASSERT(fsb, .true.)
        ASSERT(cost_best, solut%cost)
#endif
        ret = .true.
    else
        ret = .false.
    endif

end subroutine

subroutine search_two_opt(solut, info, ret) 
#ifdef DEBUG
    use assertion
#endif
    use types
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

    real(typeReal) :: cost_best 
    real(typeReal) :: cost_new
    real(typeReal) :: cost_concat_1
    real(typeReal) :: cost_concat_2 
    real(typeReal) :: cost_concat_3 
    real(typeReal) :: cost_concat_4
    real(typeReal) :: rev_seq_cost

    logical :: fsb

    cost_best = info%fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0

    do i=2, info%dimen-1
        i_prev = i-1
        rev_seq_cost = solut%seq( i+1,i)%T

        do j=i+2, info%dimen
            j_next = j+1

            rev_seq_cost = rev_seq_cost + info%cost(solut%s(j-1), solut%s(j)) * ( solut%seq( j,i)%W-real(1.0, typeReal))

            cost_concat_1 =                 solut%seq( i_prev,1)%T   + info%cost(solut%s(j), solut%s(i_prev))
            cost_concat_2 = cost_concat_1 + solut%seq( j,i)%T        + info%cost(solut%s(j_next), solut%s(i))

            cost_new = solut%seq( i_prev,1)%C                                                        + & !   1st subseq
                    solut%seq( j,i)%W                * cost_concat_1 + rev_seq_cost                  + & ! concat 2nd subseq (reversed seq)
                    solut%seq( info%dimen+1,j_next)%W * cost_concat_2 + solut%seq( info%dimen+1,j_next)%C       ! concat 3rd subseq

            if (cost_new < cost_best ) then
                cost_best = cost_new 
                I_best = i
                J_best = j
            endif

        end do

    end do

    if (cost_best < solut%cost) then
        call reverse(solut%s, I_best, J_best)
        call subseq_load(solut, info)
#ifdef DEBUG
        fsb = .false.
        call is_feasible(solut, fsb)
        ASSERT(fsb, .true.)
        ASSERT(cost_best, solut%cost)
#endif
        ret = .true.
    else
        ret = .false.
    endif

end subroutine

subroutine search_reinsertion(solut, info, opt, ret) 
#ifdef DEBUG
    use assertion
#endif
    use types
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

    real(typeReal) :: cost_best 
    real(typeReal) :: cost_new
    real(typeReal) :: cost_concat_1
    real(typeReal) :: cost_concat_2 
    real(typeReal) :: cost_concat_3 
    real(typeReal) :: cost_concat_4

    logical :: fsb

    info%reinsert_call = info%reinsert_call + 1

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
            k_next = k+1

            cost_concat_1 =                 solut%seq( k,1)%T           + info%cost(solut%s(k), solut%s(i))
            cost_concat_2 = cost_concat_1 + solut%seq( j,i)%T           + info%cost(solut%s(j), solut%s(k_next))
            cost_concat_3 = cost_concat_2 + solut%seq( i_prev,k_next)%T + info%cost(solut%s(i_prev), solut%s(j_next))

            cost_new = solut%seq( k,1)%C                                                                   + & !       1st subseq
                    solut%seq( j,i)%W                 * cost_concat_1 + solut%seq( j,i)%C                  + & ! concat 2nd subseq (reinserted seq)
                    solut%seq( i_prev,k_next)%W       * cost_concat_2 + solut%seq( i_prev,k_next)%C        + & ! concat 3rd subseq
                    solut%seq( info%dimen+1,j_next)%W * cost_concat_3 + solut%seq( info%dimen+1,j_next)%C      ! concat 4th subseq

            if (cost_new .lt. cost_best ) then
                cost_best = cost_new
                I_best = i
                J_best = j
                POS_best = k
            endif

        end do

        do k=i+opt, info%dimen
            k_next = k+1

            cost_concat_1 =                 solut%seq( i_prev,1)%T + info%cost(solut%s(i_prev), solut%s(j_next))
            cost_concat_2 = cost_concat_1 + solut%seq( k,j_next)%T + info%cost(solut%s(k), solut%s(i))
            cost_concat_3 = cost_concat_2 + solut%seq( j,i)%T      + info%cost(solut%s(j), solut%s(k_next))

            cost_new = solut%seq( i_prev,1)%C                                                              + & !       1st subseq
                    solut%seq( k,j_next)%W            * cost_concat_1 + solut%seq( k,j_next)%C             + & ! concat 2nd subseq
                    solut%seq( j,i)%W                 * cost_concat_2 + solut%seq( j,i)%C                  + & ! concat 3rd subseq (reinserted seq)
                    solut%seq( info%dimen+1,k_next)%W * cost_concat_3 + solut%seq( info%dimen+1,k_next)%C      ! concat 4th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new 
                I_best = i
                J_best = j
                POS_best = k
            endif
        end do
    end do

    if (cost_best < solut%cost) then
        call reinsert(solut%s, I_best, J_best, POS_best+1)
        call subseq_load(solut, info)

#ifdef DEBUG
        fsb = .false.
        call is_feasible(solut, fsb)
        ASSERT(fsb, .true.)
        ASSERT(cost_best, solut%cost)
#endif
        ret = .true.
    else
        ret = .false.
    endif
    !call exit(0)

end subroutine

subroutine RVND(sol, info, it)
#ifdef DEBUG
    use assertion
#endif
    use types
    implicit none

    type(tSolution) :: sol
    type(tInfo) :: info
    integer, intent(out) :: it

    real(typeReal) :: rnd

    integer, dimension(5) :: neighbd_list
    integer :: nl_size
    integer :: neighbd

    integer :: index_
    integer :: i
    integer :: total
    logical :: improve_flag
    logical :: fsb
    logical :: yes

    integer, dimension(5) :: improv

    nl_size = 5
    neighbd_list = (/ info%SWAP, info%TWO_OPT, info%REINSERTION, info%OR_OPT_2, info%OR_OPT_3 /)
    yes = .false.
    total = 0
    do while (nl_size > 0)
        index_ = info%rnd(info%rnd_index) + 1
        info%rnd_index = info%rnd_index + 1

#ifdef DEBUG
        ASSERT(index_ > nl_size, .false.)
#endif

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

#ifdef DEBUG
        ! feasibility check
        fsb = .false.
        call is_feasible(sol, fsb)
        ASSERT(fsb, .true.)
#endif

        if (improve_flag .eqv. .true.) then
            neighbd_list(1) = info%SWAP
            neighbd_list(2) = info%TWO_OPT
            neighbd_list(3) = info%REINSERTION
            neighbd_list(4) = info%OR_OPT_2
            neighbd_list(5) = info%OR_OPT_3
            nl_size = 5
            !print *, "IMPROV", sol%cost
            improv(neighbd) = improv(neighbd) + 1
            yes = .true.

        else
            call arr_shift(neighbd_list, index_+1, index_, nl_size-index_)
            nl_size = nl_size-1
            total = total + 1
        endif
        !print *, nl_size

        it = it + 1
        
    end do

end subroutine

subroutine perturb(solut_crnt, solut_part, info)! result(ret)
    use types

    implicit none
    type(tSolution) :: solut_crnt
    type(tSolution) :: solut_part
    type(tInfo) :: info

    type(tSolution) :: ret

    !integer, dimension(info%dimen+1) :: s
    integer :: A_start, A_end
    integer :: B_start, B_end
    integer :: size_max
    integer :: size_min
    integer :: max_
    real(typeReal) :: rnd

    A_start = 1
    A_end = 1
    B_start = 1
    B_end = 1

    solut_crnt%s(:) = solut_part%s(:)

    size_max = (info%dimen+1) / 10
    size_max = merge(size_max, 2, size_max >= 2)
    size_min = 2

    do while ((A_start <= B_start .and. B_start <= A_end) .or. (B_start <= A_start .and. A_start <= B_end))
       A_start = info%rnd(info%rnd_index) + 1
       info%rnd_index = info%rnd_index + 1
       A_end = A_start + info%rnd(info%rnd_index) 
       info%rnd_index = info%rnd_index + 1

       B_start = info%rnd(info%rnd_index) + 1
       info%rnd_index = info%rnd_index + 1
       B_end = B_start + info%rnd(info%rnd_index) 
       info%rnd_index = info%rnd_index + 1

    end do

    if (A_start < B_start) then
        call reinsert (solut_crnt%s, B_start, B_end - 1, A_end)
        call reinsert (solut_crnt%s, A_start, A_end - 1, B_end)
    else
        call reinsert (solut_crnt%s, A_start, A_end - 1, B_end)
        call reinsert (solut_crnt%s, B_start, B_end - 1, A_end)
    end if

end subroutine

subroutine solut_init(solut, info)
    use types
    implicit none

    type(tSolution), intent(out) :: solut
    type(tInfo) :: info

    allocate(solut%s(info%dimen+1))
    allocate(solut%seq(info%dimen+1, info%dimen+1))
    solut%s_size = info%dimen+1
    !allocate(solut%seq(3, info%dimen+1, info%dimen+1))

end subroutine

function GILS_RVND(Imax, Iils, R, info) result(ret)
#ifdef DEBUG
    use assertion
#endif
    use types

    implicit none

    !! parameters
    integer :: Imax
    integer :: Iils
    real(typeReal), dimension(26) :: R
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
    real(typeReal) :: alpha
    real(typeReal) :: rnd
    integer :: iterILS
    logical :: fsb

    integer :: it
    integer :: Tit

    interface
        function construction(alpha, info) result (ret)
            use types
            implicit none
            real(typeReal) :: alpha
            type(tInfo) :: info
            integer, dimension(info%dimen+1) :: ret 
        end function
    end interface

    call solut_init(sol_best, info)
    call solut_init(sol_partial, info)
    call solut_init(sol_crnt, info)
    sol_best%cost = info%fmax

    Tit = 0
    do i=1, Imax

        index_ = info%rnd(info%rnd_index) + 1
        info%rnd_index = info%rnd_index + 1

        alpha = R(index_)
        print *, "[+] Local Search ", i
        sol_crnt%s = construction(alpha, info)

#ifdef DEBUG
        ! feasibility check
        fsb = .false.
        call is_feasible(sol_crnt, fsb)
        ASSERT(fsb, .true.)
#endif

        call subseq_load(sol_crnt, info)
        sol_partial = sol_crnt
        print *, "        [+] Constructing Inital Solution..", sol_crnt%cost

        print *, "        [+] Looking for the best Neighbor.."
        iterILS = 0
        it = 0

        do while (iterILS < Iils)
            call RVND(sol_crnt, info, it)

            if (sol_crnt%cost < sol_partial%cost - EPSILON(real(1.0, typeReal))) then
                sol_partial%s(:) = sol_crnt%s(:)
                sol_partial%cost = sol_crnt%cost 
                iterILS = 0
            endif

            call perturb(sol_crnt, sol_partial, info)
            call subseq_load(sol_crnt, info)
#ifdef DEBUG
            fsb = .false.
            call is_feasible(sol_crnt, fsb)
            ASSERT(fsb, .true.)
#endif

            iterILS = iterILS + 1
        end do
        Tit = Tit + it

        if (sol_partial%cost < sol_best%cost) then
            sol_best = sol_partial
        endif

        print *, "        Current best solution cost:", sol_best%cost
    
    end do

    print *, sol_best%s
    !print *, "RVND it", Tit

    ret = sol_best

end function

program main
#ifdef DEBUG
    use assertion
#endif
    use data
    use types

    implicit none

    integer, allocatable :: dimensions (:)
    integer, allocatable :: rnd (:)
    type(tInfo) :: info
    type(tSolution) :: sol
    real(typeReal), dimension(26) :: r
    integer :: Iils
    integer :: Imax
    integer :: i
    real(typeReal) :: opa

    logical :: viabilidade

    INTEGER :: begin, end_, rate

    interface
        function construction(alpha, info) result (ret)
            use types
            implicit none
            real(typeReal) :: alpha
            type(tInfo) :: info
            integer, dimension(info%dimen+1) :: ret 
        end function
        function GILS_RVND(Imax, Iils, R, info) result(ret)
            use types

            implicit none

            !! parameters
            integer :: Imax 
            integer :: Iils
            real(typeReal), dimension(26) :: R
            type(tInfo) :: info
            type(tSolution) :: ret
        end function
        subroutine print_matrix(c)
            implicit none
            real(typeReal), allocatable :: c(:,:)

            integer, allocatable :: nxm(:)
            integer :: dimen
            integer :: i
            integer :: j
        end subroutine

    end interface

    call load_matrix(info%cost, info%rnd)

    dimensions = shape(info%cost(:,:))
    info%dimen = dimensions(1)

    R = (/ 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, &
    0.20, 0.21, 0.22, 0.23, 0.24, 0.25 /)
    Iils = min(100, info%dimen)
    Imax = 10

    CALL SYSTEM_CLOCK(begin, rate)
    sol = GILS_RVND(Imax, Iils, R, info)
    CALL SYSTEM_CLOCK(end_)
    print *, "COST: ", sol%cost
    print *, "TIME: ", real(end_ - begin) / real(rate)
end program
