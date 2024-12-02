module data_types
    implicit none

    type tRnd
        integer, allocatable :: rnd(:)
        integer :: rnd_index = 1
    end type

    type tInfo
        integer :: dimen
        real(8), allocatable :: cost (:,:)
        integer :: T = 1
        integer :: C = 2
        integer :: W = 3
        
        integer :: REINSERTION  = 1
        integer :: OR_OPT_2     = 2
        integer :: OR_OPT_3     = 3 
        integer :: SWAP         = 4
        integer :: TWO_OPT      = 5

        real(8) :: fmax = 1073E12
        
        integer :: reinsert_call = 0

        integer, allocatable :: rnd(:)
        integer :: rnd_index = 1
    end type

    type tSolution
        integer, allocatable :: s (:)
        real(8), allocatable :: seq (:,:,:)
        real(8) :: cost
    end type

end module

subroutine print_matrix(c)
    implicit none
    real(8), allocatable :: c(:,:)

    integer, allocatable :: nxm(:)
    integer :: dimen
    integer :: i
    integer :: j

    nxm = shape(c)
    dimen = nxm(1)

    do i=1, 2
        do j=1, dimen 
            write(*, "(F6.1, a)", advance="no") c(i,j),  ' '
        end do
        print *, new_line('A')
        print *, ''
    end do

end subroutine

subroutine print_info(info)
    use data_types
    type(tInfo) :: info

    interface
        subroutine print_matrix(c)
            real(8), allocatable :: c(:,:)
        end subroutine
    end interface

    print *, info%dimen
    call print_matrix(info%cost)

end subroutine

subroutine subseq_load(s, seq, dimen, cost)
    use data_types
    integer :: dimen
    integer, dimension(dimen+1) :: s
    real(8), dimension(dimen+1, dimen+1, 3) :: seq
    real(8), dimension(dimen, dimen) :: cost

    integer :: T = 1
    integer :: C = 2
    integer :: W = 3
    real(8) :: fmax = 1073E12

    integer :: i
    integer :: j
    integer :: k
    integer :: j_prev

    do i = 1, dimen+1
        k = 1 - i - merge(0, 1, i /= 1)

        seq(i, i, T) = 0.0
        seq(i, i, C) = 0.0
        seq(i, i, W) = merge(1.0, 0.0, i /= 1)
        do j = i+1, dimen+1
            j_prev = j-1

            seq(i,j,T) = cost(s(j_prev), s(j)) + seq(i,j_prev,T)
            seq(i,j,C) = seq(i,j,T) + seq(i,j_prev,C)
            seq(i,j,W) = j+k

        end do
    end do

end subroutine

subroutine sort(arr, len, r, dimen, cost)

    implicit none 
    integer, intent(inout) :: arr(len)
    integer :: len
    integer :: r
    integer :: dimen
    real(8), dimension(dimen, dimen) :: cost(:,:)

    interface
        recursive subroutine quicksort(arr, left, right, r, dimen, cost)

            implicit none 
            integer, intent(inout) :: arr(*)
            integer, intent(in) :: left, right, r
            integer :: dimen
            real(8), dimension(dimen, dimen) :: cost(:,:)

        end subroutine
    end interface

    call quicksort(arr, 1, len, r, dimen, cost)
end subroutine sort


function partition(arr, left, right, r, dimen, cost) result(ret)

    implicit none 
    integer, intent(inout) :: arr(*)
    integer, intent(in) :: left, right, r
    integer :: dimen
    real(8), dimension(dimen, dimen) :: cost(:,:)

    integer :: pivot, i, j, temp
    integer :: ret

    pivot = arr(right)
    i = left - 1

    do j = left, right - 1
        if (cost(r, arr(j)) < cost(r, pivot)) then
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

recursive subroutine quicksort(arr, left, right, r, dimen, cost)

    implicit none 
    integer, intent(inout) :: arr(*)
    integer, intent(in) :: left, right, r
    integer :: dimen
    real(8), dimension(dimen, dimen) :: cost(:,:)

    integer :: pivot

    interface
        function partition(arr, left, right, r, dimen, cost) result(ret)
            implicit none 
            integer, intent(inout) :: arr(*)
            integer, intent(in) :: left, right, r
            integer :: dimen
            real(8), dimension(dimen, dimen) :: cost(:,:)

            integer :: ret
        end function
    end interface

    if (left < right) then
        pivot = partition(arr, left, right, r, dimen, cost)
        call quicksort(arr, left, pivot - 1, r, dimen, cost)
        call quicksort(arr, pivot + 1, right, r, dimen, cost)
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

function construction(alpha, dimen, cost, rnd) result(ret)
    use data_types
    implicit none

    real(8) :: alpha
    integer :: dimen
    real(8), dimension(dimen, dimen) :: cost
    type(tRnd) :: rnd

    integer, dimension(dimen+1) :: ret 

    integer :: i
    integer :: j
    integer, dimension(dimen+1) :: s
    integer, dimension(dimen-1) :: cL
    integer :: cL_size
    integer :: cN
    integer :: r
    integer :: rng
    integer :: index_
    real(8) :: RND_v

    interface

        subroutine sort(arr, len, r, dimen, cost)
            implicit none 
            integer, intent(inout) :: arr(len)
            integer :: len
            integer :: r
            integer :: dimen
            real(8), dimension(dimen, dimen) :: cost(:,:)
        end subroutine

    end interface

    cL_size = dimen-1

    ! init cL
    cL = (/ (I, I=2, dimen) /)

    s(1) = 1
    r = 1
    do j=2, dimen
        call sort(cL, cL_size, r, dimen, cost)

        index_ = rnd%rnd(rnd%rnd_index) + 1
        rnd%rnd_index = rnd%rnd_index + 1

        cN = cL(index_)

        ! memmove on cL
        call arr_shift(cL, index_+1, index_, cL_size-index_)
        cL_size = cL_size-1

        s(j) = cN
        r = cN
    end do

    s(dimen+1) = 1

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

subroutine search_swap(s, seq, dimen, cost, ret) 
    use data_types
    implicit none

    integer :: dimen
    integer, dimension(dimen+1) :: s
    real(8), dimension(dimen+1, dimen+1, 3) :: seq
    real(8), dimension(dimen, dimen) :: cost

    logical, intent(out) :: ret

    integer :: T = 1
    integer :: C = 2
    integer :: W = 3
    real(8) :: fmax = 1073E12

    integer :: i
    integer :: j
    integer :: i_prev
    integer :: i_next
    integer :: j_prev
    integer :: j_next

    integer :: I_best
    integer :: J_best

    real(8) :: cost_best 
    real(8) :: cost_new
    real(8) :: cost_concat_1
    real(8) :: cost_concat_2 
    real(8) :: cost_concat_3 
    real(8) :: cost_concat_4

    interface
        subroutine subseq_load(s, seq, dimen, cost)
            use data_types
            integer :: dimen
            integer, dimension(dimen+1) :: s
            real(8), dimension(dimen+1, dimen+1, 3) :: seq
            real(8), dimension(dimen, dimen) :: cost
        end subroutine
    end interface

    cost_best = fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0
    cost_concat_4 = 0.0

    do i=2, dimen-1
        i_prev = i-1
        i_next = i+1

        cost_concat_1 =                 seq(1, i_prev, T) + cost(s(i_prev), s(i_next))
        cost_concat_2 = cost_concat_1 + seq(i, i_next, T) + cost(s(i), s(i_next+1))

        cost_new = seq(1, i_prev, C)                                                    + &
                seq(i, i_next, W)               * (cost_concat_1) + cost(s(i_next), s(i))  + &
                seq(i_next+1, dimen+1, W)   * (cost_concat_2) + seq(i_next+1, dimen+1, C) 

        if (cost_new < cost_best) then
            cost_best = cost_new 
            I_best = i
            J_best = i_next
        endif

        do j=i_next+1, dimen

            j_prev = j-1
            j_next = j+1


            cost_concat_1 =                 seq(1, i_prev, T)       + cost(s(i_prev), s(j))
            cost_concat_2 = cost_concat_1                           + cost(s(j), s(i_next))
            cost_concat_3 = cost_concat_2 + seq(i_next, j_prev, T)  + cost(s(j_prev), s(i))
            cost_concat_4 = cost_concat_3                           + cost(s(i), s(j_next))


            cost_new = seq(1, i_prev, C)                                                 + &     ! 1st subseq
                    cost_concat_1 + &                                                           ! concat 2nd subseq (single node)
                    seq(i_next, j_prev, W)      * cost_concat_2 + seq(i_next, j_prev, C) + &    ! concat 3rd subseq
                    cost_concat_3 + &                                                           ! concat 4th subseq (single node)
                    seq(j_next, dimen+1, W) * cost_concat_4 + seq(j_next, dimen+1, C)   ! concat 5th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new 
                I_best = i;
                J_best = j;
            endif

        end do
        
    end do

    if (cost_best < seq(1, dimen+1, C) - EPSILON(1.0)) then
        call swap(s, I_best, J_best)
        call subseq_load(s, seq, dimen, cost)
        ret = .true.
    else
        ret = .false.
    endif

end subroutine

subroutine search_two_opt(s, seq, dimen, cost, ret) 
    use data_types
    implicit none

    integer :: dimen
    integer, dimension(dimen+1) :: s 
    real(8), dimension(dimen+1, dimen+1, 3) :: seq
    real(8), dimension(dimen, dimen) :: cost 
    logical, intent(out) :: ret

    integer :: T = 1
    integer :: C = 2
    integer :: W = 3
    real(8) :: fmax = 1073E12

    integer :: i
    integer :: j
    integer :: i_prev
    integer :: j_next

    integer :: I_best
    integer :: J_best

    real(8) :: cost_best 
    real(8) :: cost_new
    real(8) :: cost_concat_1
    real(8) :: cost_concat_2 
    real(8) :: cost_concat_3 
    real(8) :: cost_concat_4
    real(8) :: rev_seq_cost

    interface
        subroutine subseq_load(s, seq, dimen, cost)
            use data_types
            integer :: dimen
            integer, dimension(dimen+1) :: s
            real(8), dimension(dimen+1, dimen+1, 3) :: seq
            real(8), dimension(dimen, dimen) :: cost
        end subroutine
    end interface

    cost_best = fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0

    do i=2, dimen-1
        i_prev = i-1
        rev_seq_cost = seq(i, i+1, T)

        do j=i+2, dimen
            j_next = j+1

            rev_seq_cost = rev_seq_cost + cost(s(j-1), s(j)) * (seq(i, j, W)-1.0)

            cost_concat_1 =                 seq(1, i_prev, T)   + cost(s(j), s(i_prev))
            cost_concat_2 = cost_concat_1 + seq(i, j, T)        + cost(s(j_next), s(i))

            cost_new = seq(1, i_prev, C)                                                        + & !   1st subseq
                    seq(i, j, W)                * cost_concat_1 + rev_seq_cost                  + & ! concat 2nd subseq (reversed seq)
                    seq(j_next, dimen+1, W) * cost_concat_2 + seq(j_next, dimen+1, C)       ! concat 3rd subseq

            if (cost_new < cost_best) then
                cost_best = cost_new 
                I_best = i
                J_best = j
            endif

        end do

    end do

    if (cost_best < seq(1, dimen+1, C)) then
        call reverse(s, I_best, J_best)
        call subseq_load(s, seq, dimen, cost)
        ret = .true.
    else
        ret = .false.
    endif

end subroutine

subroutine search_reinsertion(s, seq, dimen, cost, opt, ret) 
    use data_types
    implicit none

    integer :: dimen
    integer, dimension(dimen+1) :: s
    real(8), dimension(dimen+1, dimen+1, 3) :: seq
    real(8), dimension(dimen, dimen) :: cost
    integer :: opt
    logical, intent(out) :: ret

    integer :: T = 1
    integer :: C = 2
    integer :: W = 3
    real(8) :: fmax = 1073E12

    integer :: i
    integer :: j
    integer :: k
    integer :: i_prev
    integer :: j_next
    integer :: k_next

    integer :: I_best
    integer :: J_best
    integer :: POS_best

    real(8) :: cost_best 
    real(8) :: cost_new
    real(8) :: cost_concat_1
    real(8) :: cost_concat_2 
    real(8) :: cost_concat_3 
    real(8) :: cost_concat_4

    interface
        subroutine subseq_load(s, seq, dimen, cost)
            use data_types
            integer :: dimen
            integer, dimension(dimen+1) :: s
            real(8), dimension(dimen+1, dimen+1, 3) :: seq
            real(8), dimension(dimen, dimen) :: cost
        end subroutine

    end interface

    cost_best = fmax
    cost_new = 0.0
    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0

    do i=2, dimen-opt+1
        j = opt+i-1
        i_prev = i-1
        j_next = j+1

        do k=1, i_prev-1
            k_next = k+1

            cost_new =  seq(1, dimen+1, C)
            cost_concat_1 =                 seq(1, k, T)            + cost(s(k), s(i))
            cost_concat_2 = cost_concat_1 + seq(i, j, T)            + cost(s(j), s(k_next))
            cost_concat_3 = cost_concat_2 + seq(k_next, i_prev, T)  + cost(s(i_prev), s(j_next))

            cost_new = seq(1, k, C)                                                             + & !       1st subseq
                    seq(i, j, W)                * cost_concat_1 + seq(i, j, C)                  + & ! concat 2nd subseq (reinserted seq)
                    seq(k_next, i_prev, W)      * cost_concat_2 + seq(k_next, i_prev, C)        + & ! concat 3rd subseq
                    seq(j_next, dimen+1, W) * cost_concat_3 + seq(j_next, dimen+1, C)       ! concat 4th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new 
                I_best = i
                J_best = j
                POS_best = k
            endif

        end do

        do k=i+opt, dimen
            k_next = k+1

            cost_concat_1 =                 seq(1, i_prev, T)   + cost(s(i_prev), s(j_next))
            cost_concat_2 = cost_concat_1 + seq(j_next, k, T)   + cost(s(k), s(i))
            cost_concat_3 = cost_concat_2 + seq(i, j, T)        + cost(s(j), s(k_next))

            cost_new = seq(1, i_prev, C)                                                        + & !       1st subseq
                    seq(j_next, k, W)           * cost_concat_1 + seq(j_next, k, C)             + & ! concat 2nd subseq
                    seq(i, j, W)                * cost_concat_2 + seq(i, j, C)                  + & ! concat 3rd subseq (reinserted seq)
                    seq(k_next, dimen+1, W) * cost_concat_3 + seq(k_next, dimen+1, C)       ! concat 4th subseq

            if (cost_new < cost_best) then
                cost_best = cost_new 
                I_best = i
                J_best = j
                POS_best = k
            endif
        end do
    end do

    if (cost_best < seq(1, dimen+1, C)) then
        call reinsert(s, I_best, J_best, POS_best+1)
        call subseq_load(s, seq, dimen, cost)
        ret = .true.
    else
        ret = .false.
    endif

end subroutine

subroutine RVND(sol, seq, dimen, cost, rnd)
    use data_types
    implicit none

    integer :: dimen
    integer, dimension(dimen+1) :: sol
    real(8), dimension(dimen+1, dimen+1, 3) :: seq
    real(8), dimension(dimen, dimen) :: cost
    type(tRnd) :: rnd

    real(8) :: RND_v


    integer, dimension(5) :: neighbd_list
    integer :: nl_size
    integer :: neighbd

    integer :: index_
    integer :: i
    logical :: improve_flag

    integer :: REINSERTION  = 1
    integer :: OR_OPT_2     = 2
    integer :: OR_OPT_3     = 3 
    integer :: SWAP         = 4
    integer :: TWO_OPT      = 5

    interface
        subroutine search_reinsertion(s, seq, dimen, cost, opt, ret) 
            use data_types
            implicit none

            integer :: dimen
            integer, dimension(dimen+1) :: s
            real(8), dimension(dimen+1, dimen+1, 3) :: seq
            real(8), dimension(dimen, dimen) :: cost
            integer :: opt
            logical, intent(out) :: ret
        end subroutine 
        subroutine search_two_opt(s, seq, dimen, cost, ret) 
            use data_types
            implicit none

            integer :: dimen
            integer, dimension(dimen+1) :: s
            real(8), dimension(dimen+1, dimen+1, 3) :: seq
            real(8), dimension(dimen, dimen) :: cost

            logical, intent(out) :: ret
        end subroutine
        subroutine search_swap(s, seq, dimen, cost, ret) 
            use data_types
            implicit none

            integer :: dimen
            integer, dimension(dimen+1) :: s
            real(8), dimension(dimen+1, dimen+1, 3) :: seq
            real(8), dimension(dimen, dimen) :: cost

            logical, intent(out) :: ret

        end subroutine


    end interface

    rnd%rnd(rnd%rnd_index) = rnd%rnd(rnd%rnd_index)

    nl_size = 5
    neighbd_list = (/ SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3 /)
    do while (nl_size > 0)

        index_ = rnd%rnd(rnd%rnd_index) + 1
        rnd%rnd_index = rnd%rnd_index + 1

        neighbd = neighbd_list(index_)
        
        improve_flag = .false.
        if (neighbd == REINSERTION) then
            call search_reinsertion(sol, seq, dimen, cost,  REINSERTION, improve_flag)
        else if (neighbd == OR_OPT_2) then
            call search_reinsertion(sol, seq, dimen, cost,  OR_OPT_2, improve_flag)
        else if (neighbd == OR_OPT_3) then
            call search_reinsertion(sol, seq, dimen, cost, OR_OPT_3, improve_flag)
        else if (neighbd == SWAP) then
            call search_swap(sol, seq, dimen, cost, improve_flag)
        else if (neighbd == TWO_OPT) then
            call search_two_opt(sol, seq, dimen, cost, improve_flag)
        endif

        if (improve_flag .eqv. .true.) then
            neighbd_list = (/ SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3 /)
            nl_size = 5
        else
            call arr_shift(neighbd_list, index_+1, index_, nl_size-index_)
            nl_size = nl_size-1
        endif

    end do

end subroutine

function perturb(sl, dimen, rnd) result(ret)
    use data_types

    implicit none
    integer :: dimen
    integer, dimension(dimen+1) :: sl
    type(tRnd) :: rnd

    integer, dimension(dimen+1) :: ret

    integer, dimension(dimen+1) :: s
    integer :: A_start, A_end
    integer :: B_start, B_end
    integer :: size_max
    integer :: size_min
    integer :: max_
    real(8) :: RND_v

    A_start = 1
    A_end = 1
    B_start = 1
    B_end = 1

    s = sl

    do while ((A_start <= B_start .and. B_start <= A_end) .or. (B_start <= A_start .and. A_start <= B_end))

       A_start = rnd%rnd(rnd%rnd_index) + 1
       rnd%rnd_index = rnd%rnd_index + 1
       A_end = A_start + rnd%rnd(rnd%rnd_index) 
       rnd%rnd_index = rnd%rnd_index + 1

       B_start = rnd%rnd(rnd%rnd_index) + 1
       rnd%rnd_index = rnd%rnd_index + 1
       B_end = B_start + rnd%rnd(rnd%rnd_index) 
       rnd%rnd_index = rnd%rnd_index + 1

    end do

    if (A_start < B_start) then
        call reinsert (s, B_start, B_end - 1, A_end)
        call reinsert (s, A_start, A_end - 1, B_end)
    else
        call reinsert (s, A_start, A_end - 1, B_end)
        call reinsert (s, B_start, B_end - 1, A_end)
    end if

    ret = s

end function

subroutine GILS_RVND(Imax, Iils, R, dimen, cost, rnd, ret)
    use data_types

    implicit none

    !! parameters
    integer :: Imax
    integer :: Iils
    real(8), dimension(26) :: R
    integer :: dimen
    real(8), allocatable :: cost (:,:)
    type(tRnd) :: rnd

    integer, allocatable, intent(out) :: ret(:)

    !! solutions
    integer, dimension(dimen+1) :: sol_best
    integer, dimension(dimen+1) :: sol_partial
    integer, dimension(dimen+1) :: sol_crnt

    real(8), dimension(dimen+1, dimen+1, 3) :: seq

    real(8) :: cost_best
    real(8) :: cost_partial
    real(8) :: cost_crnt

    !! aux variables
    integer :: i
    integer :: index_
    integer :: R_size = 26
    real(8) :: alpha
    real(8) :: RND_v
    integer :: iterILS

    integer :: C = 2

    real(8) :: fmax = 1073E12

    interface
        function construction(alpha, dimen, cost, rnd) result(ret)
            use data_types
            implicit none
            real(8) :: alpha
            type(tRnd) :: rnd
            integer :: dimen
            real(8), dimension(dimen, dimen) :: cost 
            integer, dimension(dimen+1) :: ret 
        end function
        function perturb(sl, dimen, rnd) result(ret)
            use data_types

            implicit none
            integer :: dimen
            integer, dimension(dimen+1) :: sl
            type(tRnd) :: rnd

            integer, dimension(dimen+1) :: ret
        end function
    end interface

    cost_best = fmax

    do i=1, Imax

        index_ = rnd%rnd(rnd%rnd_index) + 1
        rnd%rnd_index = rnd%rnd_index + 1

        alpha = R(index_)
        print *, "[+] Local Search ", i
        sol_crnt = construction(alpha, dimen, cost, rnd)
        sol_partial = sol_crnt

        call subseq_load(sol_crnt, seq, dimen, cost)
        cost_crnt = seq(1, dimen+1, C)
        cost_partial = cost_crnt

        print *, "        [+] Constructing Inital Solution..", cost_crnt
        print *, "        [+] Looking for the best Neighbor.."

        iterILS = 0

        do while (iterILS < Iils)
            call RVND(sol_crnt, seq, dimen, cost, rnd)

            cost_crnt = seq(1, dimen+1, C)
            if (cost_crnt < cost_partial - EPSILON(1.0)) then
                sol_partial = sol_crnt
                cost_partial = cost_crnt 
                iterILS = 0
            endif

            sol_crnt = perturb(sol_partial, dimen, rnd)
            call subseq_load(sol_crnt, seq, dimen, cost)

            iterILS = iterILS + 1
        end do

        if (cost_partial < cost_best) then
            sol_best = sol_partial
            cost_best = cost_partial
        endif

        print *, "        Current best solution cost:", cost_best
    
    end do

    print *, sol_best
    print *, "COST: ", cost_best

    ret = sol_best(:)

end subroutine

program main
    use rData
    use data_types

    implicit none

    type(tInfo) :: info
    type(tRnd) :: rnd
    real(8), dimension(26) :: r
    integer :: Iils
    integer :: Imax
    integer :: i

    integer :: dimen
    real(8), allocatable :: cost (:,:)
    integer, allocatable :: sol(:)

    INTEGER :: begin, end_, rate

    interface
        function construction(alpha, dimen, cost, rnd) result(ret)
            use data_types
            implicit none
            real(8) :: alpha
            type(tRnd) :: rnd
            integer :: dimen
            real(8), dimension(dimen, dimen) :: cost 
            integer, dimension(dimen+1) :: ret 
        end function
        subroutine GILS_RVND(Imax, Iils, R, dimen, cost, rnd, ret)
            use data_types

            implicit none

            !! parameters
            integer :: Imax
            integer :: Iils
            real(8), dimension(26) :: R
            integer :: dimen
            real(8), allocatable :: cost (:,:)
            type(tRnd) :: rnd
            integer, allocatable, intent(out) :: ret(:)
        end subroutine
        subroutine print_matrix(c)
            implicit none
            real(8), allocatable :: c(:,:)

            integer, allocatable :: nxm(:)
            integer :: dimen
            integer :: i
            integer :: j
        end subroutine

    end interface

    call load_matrix(cost, rnd%rnd, dimen)
    rnd%rnd_index = 1

    allocate(sol(dimen))

    R = (/ (i/100.0 + 0.000000001, i=1, 26) /)
    Iils = min(100, dimen)
    Imax = 10
    CALL SYSTEM_CLOCK(begin, rate)
    call GILS_RVND(Imax, Iils, R, dimen, cost, rnd, sol)
    CALL SYSTEM_CLOCK(end_)

    print *, "TIME: ", real(end_ - begin) / real(rate)

end program
