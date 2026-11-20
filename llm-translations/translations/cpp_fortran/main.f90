! Traducao de main.cpp

module mlp_mod
    use data_mod, only: dp, loadData
    implicit none

    integer, parameter :: TABLE_SZ = 26
    integer, parameter :: DBL_SZ   = 8
    integer, parameter :: INT_SZ   = 4

    real(dp), parameter :: DBL_MAX_ = huge(1.0_dp)
    real(dp), parameter :: DBL_EPSILON_ = epsilon(1.0_dp)

    type :: tData
        real(dp), allocatable :: cost(:,:)
        integer :: dimen
        integer, allocatable :: rnd(:)
        integer :: rnd_index
    end type tData

    type :: tSeqInfo
        real(dp) :: T
        real(dp) :: C
        real(dp) :: W
    end type tSeqInfo

    ! enum class Neighborhood
    integer, parameter :: REINSERTION = 1
    integer, parameter :: OR_OPT_2    = 2
    integer, parameter :: OR_OPT_3    = 3
    integer, parameter :: SWAP        = 4
    integer, parameter :: TWO_OPT     = 5

    type :: tSolution
        integer, allocatable :: s(:)
        type(tSeqInfo), allocatable :: seq(:)
        real(dp) :: cost
    end type tSolution

contains

    pure function to_1D(x, y, sz) result(a)
        integer, intent(in) :: x, y, sz
        integer :: a

        a = x*(sz+1) + y
    end function to_1D

    function Solution_init(data) result(solut)
        type(tData), intent(in) :: data
        type(tSolution) :: solut

        allocate(solut%s(0:data%dimen))
        solut%s = 0

        allocate(solut%seq(0:(data%dimen+1)*(data%dimen+1) - 1))

        solut%cost = DBL_MAX_
    end function Solution_init

    subroutine Solution_cpy(src, tgt, data)
        type(tSolution), intent(in) :: src
        type(tSolution), intent(inout) :: tgt
        type(tData), intent(in) :: data

        tgt%s = src%s
        tgt%cost = src%cost
    end subroutine Solution_cpy

    function R_table(i) result(v)
        integer, intent(in) :: i
        real(dp) :: v

        real(dp), parameter :: table(0:TABLE_SZ-1) = &
            [0.00_dp, 0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp, 0.07_dp, 0.08_dp, 0.09_dp, &
             0.10_dp, 0.11_dp, 0.12_dp, 0.13_dp, 0.14_dp, 0.15_dp, 0.16_dp, 0.17_dp, 0.18_dp, 0.19_dp, &
             0.20_dp, 0.21_dp, 0.22_dp, 0.23_dp, 0.24_dp, 0.25_dp]

        v = table(i)
    end function R_table

    subroutine print_s(s)
        integer, intent(in) :: s(0:)
        integer :: i

        do i = 0, size(s)-1
            write(*, '(I0,A)', advance='no') s(i)+1, ' '
        end do
        write(*, '(A)') ''
    end subroutine print_s

    function partition(arr, left, right, data, r) result(p)
        integer, intent(inout) :: arr(0:)
        integer, intent(in) :: left, right, r
        type(tData), intent(in) :: data
        integer :: p

        integer :: pivot, i, j, tmp

        pivot = arr(right)
        i = left - 1
        do j = left, right-1
            if (data%cost(r, arr(j)) < data%cost(r, pivot)) then
                i = i + 1
                tmp = arr(i); arr(i) = arr(j); arr(j) = tmp
            end if
        end do
        tmp = arr(i+1); arr(i+1) = arr(right); arr(right) = tmp
        p = i + 1
    end function partition

    recursive subroutine quicksort(arr, left, right, data, r)
        integer, intent(inout) :: arr(0:)
        integer, intent(in) :: left, right, r
        type(tData), intent(in) :: data

        integer :: pivot

        if (left < right) then
            pivot = partition(arr, left, right, data, r)
            call quicksort(arr, left, pivot-1, data, r)
            call quicksort(arr, pivot+1, right, data, r)
        end if
    end subroutine quicksort

    subroutine sort(arr, n, r, data)
        integer, intent(inout) :: arr(0:)
        integer, intent(in) :: n, r
        type(tData), intent(in) :: data

        call quicksort(arr, 0, n-1, data, r)
    end subroutine sort

    function construct(alpha, data) result(s)
        real(dp), intent(in) :: alpha
        type(tData), intent(inout) :: data
        integer, allocatable :: s(:)

        integer, allocatable :: cL(:)
        integer :: cl_size, s_size, i, r, index, c

        allocate(s(0:data%dimen))
        s(0) = 0
        s_size = 1

        allocate(cL(0:data%dimen-2))
        do i = 1, data%dimen-1
            cL(i-1) = i
        end do
        cl_size = data%dimen-1

        r = 0
        do while (cl_size > 0)
            call sort(cL, cl_size, r, data)

            index = data%rnd(data%rnd_index)
            data%rnd_index = data%rnd_index + 1

            c = cL(index)
            s(s_size) = c
            s_size = s_size + 1
            r = c

            ! cL.erase(cL.begin() + index)
            if (index < cl_size-1) cL(index:cl_size-2) = cL(index+1:cl_size-1)
            cl_size = cl_size - 1
        end do

        s(s_size) = 0
    end function construct

    subroutine swap_(vec, i, j)
        integer, intent(inout) :: vec(0:)
        integer, intent(in) :: i, j
        integer :: tmp

        tmp = vec(i); vec(i) = vec(j); vec(j) = tmp
    end subroutine swap_

    subroutine reverse_(vec, i, j)
        integer, intent(inout) :: vec(0:)
        integer, intent(in) :: i, j

        vec(i:j) = vec(j:i:-1)
    end subroutine reverse_

    subroutine reinsert(vec, i, j, pos)
        integer, intent(inout) :: vec(0:)
        integer, intent(in) :: i, j, pos

        integer :: seq(0:j-i)
        integer :: sz

        sz = j - i + 1
        seq = vec(i:j)

        if (pos < i) then
            vec(pos+sz : i-1+sz) = vec(pos : i-1)
            vec(pos : pos+sz-1) = seq
        else
            vec(i : pos-sz-1) = vec(j+1 : pos-1)
            vec(pos-sz : pos-1) = seq
        end if
    end subroutine reinsert

    subroutine update_subseq_info_matrix(solut, data, index)
        type(tSolution), intent(inout) :: solut
        type(tData), intent(in) :: data
        integer, intent(in), optional :: index

        integer :: i, j, j_prev, k
        integer :: from
        logical :: t

        from = 0
        if (present(index)) from = index

        do i = 0, data%dimen
            if (i == 0) then
                k = 0
            else
                k = 1 - i
            end if
            t = (i == from)

            solut%seq(to_1D(i, i, data%dimen))%T = 0.0_dp
            solut%seq(to_1D(i, i, data%dimen))%C = 0.0_dp
            if (i == 0) then
                solut%seq(to_1D(i, i, data%dimen))%W = 0.0_dp
            else
                solut%seq(to_1D(i, i, data%dimen))%W = 1.0_dp
            end if

            do j = i+1, data%dimen
                j_prev = j-1

                solut%seq(to_1D(i, j, data%dimen))%T = data%cost(solut%s(j_prev), solut%s(j)) &
                    + solut%seq(to_1D(i, j_prev, data%dimen))%T
                solut%seq(to_1D(i, j, data%dimen))%C = solut%seq(to_1D(i, j, data%dimen))%T &
                    + solut%seq(to_1D(i, j_prev, data%dimen))%C
                solut%seq(to_1D(i, j, data%dimen))%W = real(j + k, dp)
            end do

            if (t) from = from + 1
        end do

        solut%cost = solut%seq(to_1D(0, data%dimen, data%dimen))%C
    end subroutine update_subseq_info_matrix

    function search_swap(solut, data) result(improved)
        type(tSolution), intent(inout) :: solut
        type(tData), intent(in) :: data
        logical :: improved

        real(dp) :: cost_new, cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4
        real(dp) :: cost_best
        integer :: i, j, j_prev, j_next, i_prev, i_next
        integer :: I_, J_

        cost_best = DBL_MAX_
        I_ = 0
        J_ = 0

        do i = 1, data%dimen-2
            i_prev = i - 1
            i_next = i + 1

            ! nos consecutivos
            cost_concat_1 =                 solut%seq(to_1D(0, i_prev, data%dimen))%T &
                + data%cost(solut%s(i_prev), solut%s(i_next))
            cost_concat_2 = cost_concat_1 + solut%seq(to_1D(i, i_next, data%dimen))%T &
                + data%cost(solut%s(i), solut%s(i_next+1))

            cost_new = solut%seq(to_1D(0, i_prev, data%dimen))%C &
                + solut%seq(to_1D(i, i_next, data%dimen))%W * cost_concat_1 &
                    + data%cost(solut%s(i_next), solut%s(i)) &
                + solut%seq(to_1D(i_next+1, data%dimen, data%dimen))%W * cost_concat_2 &
                    + solut%seq(to_1D(i_next+1, data%dimen, data%dimen))%C

            if (cost_new < cost_best) then
                cost_best = cost_new - DBL_EPSILON_
                I_ = i
                J_ = i_next
            end if

            do j = i_next+1, data%dimen-1
                j_next = j + 1
                j_prev = j - 1

                cost_concat_1 =                 solut%seq(to_1D(0, i_prev, data%dimen))%T &
                    + data%cost(solut%s(i_prev), solut%s(j))
                cost_concat_2 = cost_concat_1 + data%cost(solut%s(j), solut%s(i_next))
                cost_concat_3 = cost_concat_2 + solut%seq(to_1D(i_next, j_prev, data%dimen))%T &
                    + data%cost(solut%s(j_prev), solut%s(i))
                cost_concat_4 = cost_concat_3 + data%cost(solut%s(i), solut%s(j_next))

                cost_new = solut%seq(to_1D(0, i_prev, data%dimen))%C &
                    + cost_concat_1 &
                    + solut%seq(to_1D(i_next, j_prev, data%dimen))%W * cost_concat_2 &
                        + solut%seq(to_1D(i_next, j_prev, data%dimen))%C &
                    + cost_concat_3 &
                    + solut%seq(to_1D(j_next, data%dimen, data%dimen))%W * cost_concat_4 &
                        + solut%seq(to_1D(j_next, data%dimen, data%dimen))%C

                if (cost_new < cost_best) then
                    cost_best = cost_new - DBL_EPSILON_
                    I_ = i
                    J_ = j
                end if
            end do
        end do

        if (cost_best < solut%cost - DBL_EPSILON_) then
            call swap_(solut%s, I_, J_)
            call update_subseq_info_matrix(solut, data, I_)
            improved = .true.
            return
        end if

        improved = .false.
    end function search_swap

    function search_two_opt(solut, data) result(improved)
        type(tSolution), intent(inout) :: solut
        type(tData), intent(in) :: data
        logical :: improved

        real(dp) :: cost_new, cost_concat_1, cost_concat_2
        real(dp) :: cost_best
        real(dp) :: rev_seq_cost
        integer :: i, j, i_prev, j_next
        integer :: I_, J_

        cost_best = DBL_MAX_
        I_ = 0
        J_ = 0

        do i = 1, data%dimen-2
            i_prev = i - 1

            rev_seq_cost = solut%seq(to_1D(i, i+1, data%dimen))%T
            do j = i+2, data%dimen-1
                j_next = j + 1

                rev_seq_cost = rev_seq_cost + data%cost(solut%s(j-1), solut%s(j)) &
                    * (solut%seq(to_1D(i, j, data%dimen))%W - 1.0_dp)

                cost_concat_1 =                 solut%seq(to_1D(0, i_prev, data%dimen))%T &
                    + data%cost(solut%s(j), solut%s(i_prev))
                cost_concat_2 = cost_concat_1 + solut%seq(to_1D(i, j, data%dimen))%T &
                    + data%cost(solut%s(j_next), solut%s(i))

                cost_new = solut%seq(to_1D(0, i_prev, data%dimen))%C &
                    + solut%seq(to_1D(i, j, data%dimen))%W * cost_concat_1 + rev_seq_cost &
                    + solut%seq(to_1D(j_next, data%dimen, data%dimen))%W * cost_concat_2 &
                        + solut%seq(to_1D(j_next, data%dimen, data%dimen))%C

                if (cost_new < cost_best) then
                    cost_best = cost_new - DBL_EPSILON_
                    I_ = i
                    J_ = j
                end if
            end do
        end do

        if (cost_best < solut%cost - DBL_EPSILON_) then
            call reverse_(solut%s, I_, J_)
            call update_subseq_info_matrix(solut, data)
            improved = .true.
            return
        end if

        improved = .false.
    end function search_two_opt

    function search_reinsertion(solut, data, opt) result(improved)
        type(tSolution), intent(inout) :: solut
        type(tData), intent(in) :: data
        integer, intent(in) :: opt
        logical :: improved

        real(dp) :: cost_new, cost_concat_1, cost_concat_2, cost_concat_3
        real(dp) :: cost_best
        integer :: i, j, k, k_next, i_prev, j_next
        integer :: I_, J_, POS

        cost_best = DBL_MAX_
        I_ = 0
        J_ = 0
        POS = 0

        do i = 1, data%dimen-opt
            j = opt + i - 1
            j_next = j + 1
            i_prev = i - 1

            ! k -> arestas
            do k = 0, i_prev-1
                k_next = k + 1

                cost_concat_1 =                 solut%seq(to_1D(0, k, data%dimen))%T &
                    + data%cost(solut%s(k), solut%s(i))
                cost_concat_2 = cost_concat_1 + solut%seq(to_1D(i, j, data%dimen))%T &
                    + data%cost(solut%s(j), solut%s(k_next))
                cost_concat_3 = cost_concat_2 + solut%seq(to_1D(k_next, i_prev, data%dimen))%T &
                    + data%cost(solut%s(i_prev), solut%s(j_next))

                cost_new = solut%seq(to_1D(0, k, data%dimen))%C &
                    + solut%seq(to_1D(i, j, data%dimen))%W * cost_concat_1 &
                        + solut%seq(to_1D(i, j, data%dimen))%C &
                    + solut%seq(to_1D(k_next, i_prev, data%dimen))%W * cost_concat_2 &
                        + solut%seq(to_1D(k_next, i_prev, data%dimen))%C &
                    + solut%seq(to_1D(j_next, data%dimen, data%dimen))%W * cost_concat_3 &
                        + solut%seq(to_1D(j_next, data%dimen, data%dimen))%C

                if (cost_new < cost_best) then
                    cost_best = cost_new - DBL_EPSILON_
                    I_ = i
                    J_ = j
                    POS = k
                end if
            end do

            do k = i+opt, data%dimen-1
                k_next = k + 1

                cost_concat_1 =                 solut%seq(to_1D(0, i_prev, data%dimen))%T &
                    + data%cost(solut%s(i_prev), solut%s(j_next))
                cost_concat_2 = cost_concat_1 + solut%seq(to_1D(j_next, k, data%dimen))%T &
                    + data%cost(solut%s(k), solut%s(i))
                cost_concat_3 = cost_concat_2 + solut%seq(to_1D(i, j, data%dimen))%T &
                    + data%cost(solut%s(j), solut%s(k_next))

                cost_new = solut%seq(to_1D(0, i_prev, data%dimen))%C &
                    + solut%seq(to_1D(j_next, k, data%dimen))%W * cost_concat_1 &
                        + solut%seq(to_1D(j_next, k, data%dimen))%C &
                    + solut%seq(to_1D(i, j, data%dimen))%W * cost_concat_2 &
                        + solut%seq(to_1D(i, j, data%dimen))%C &
                    + solut%seq(to_1D(k_next, data%dimen, data%dimen))%W * cost_concat_3 &
                        + solut%seq(to_1D(k_next, data%dimen, data%dimen))%C

                if (cost_new < cost_best) then
                    cost_best = cost_new - DBL_EPSILON_
                    I_ = i
                    J_ = j
                    POS = k
                end if
            end do
        end do

        if (cost_best < solut%cost - DBL_EPSILON_) then
            call reinsert(solut%s, I_, J_, POS+1)
            call update_subseq_info_matrix(solut, data, min(I_, POS+1))
            improved = .true.
            return
        end if

        improved = .false.
    end function search_reinsertion

    subroutine RVND(solut, data)
        type(tSolution), intent(inout) :: solut
        type(tData), intent(inout) :: data

        integer :: neighbd_list(0:4)
        integer :: list_size
        integer :: index
        integer :: neighbd
        logical :: improve_flag

        neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]
        list_size = 5

        do while (list_size > 0)

            index = data%rnd(data%rnd_index)
            data%rnd_index = data%rnd_index + 1
            neighbd = neighbd_list(index)
            improve_flag = .false.

            select case (neighbd)
            case (REINSERTION)
                improve_flag = search_reinsertion(solut, data, REINSERTION)
            case (OR_OPT_2)
                improve_flag = search_reinsertion(solut, data, OR_OPT_2)
            case (OR_OPT_3)
                improve_flag = search_reinsertion(solut, data, OR_OPT_3)
            case (SWAP)
                improve_flag = search_swap(solut, data)
            case (TWO_OPT)
                improve_flag = search_two_opt(solut, data)
            end select

            if (improve_flag) then
                neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]
                list_size = 5
            else
                if (index < list_size-1) neighbd_list(index:list_size-2) = neighbd_list(index+1:list_size-1)
                list_size = list_size - 1
            end if

        end do

    end subroutine RVND

    function perturb(solut, data) result(s)
        type(tSolution), intent(in) :: solut
        type(tData), intent(inout) :: data
        integer, allocatable :: s(:)

        integer :: A_start, A_end, B_start, B_end
        integer :: size_max, size_min

        allocate(s(0:data%dimen))
        s = solut%s

        A_start = 1
        A_end = 1
        B_start = 1
        B_end = 1

        size_max = (data%dimen+1) / 10
        if (size_max < 2) size_max = 2
        size_min = 2

        do while ((A_start <= B_start .and. B_start <= A_end) .or. (B_start <= A_start .and. A_start <= B_end))

            A_start = data%rnd(data%rnd_index)
            data%rnd_index = data%rnd_index + 1
            A_end = A_start + data%rnd(data%rnd_index)
            data%rnd_index = data%rnd_index + 1

            B_start = data%rnd(data%rnd_index)
            data%rnd_index = data%rnd_index + 1
            B_end = B_start + data%rnd(data%rnd_index)
            data%rnd_index = data%rnd_index + 1
        end do

        if (A_start < B_start) then
            call reinsert(s, B_start, B_end-1, A_end)
            call reinsert(s, A_start, A_end-1, B_end)
        else
            call reinsert(s, A_start, A_end-1, B_end)
            call reinsert(s, B_start, B_end-1, A_end)
        end if

    end function perturb

    subroutine GILS_RVND(Imax, Iils, data, verbose)
        integer, intent(in) :: Imax, Iils
        type(tData), intent(inout) :: data
        logical, intent(in) :: verbose

        type(tSolution) :: solut_partial, solut_crnt, solut_best
        integer :: i, aux, iterILS
        real(dp) :: alpha

        solut_partial = Solution_init(data)
        solut_crnt    = Solution_init(data)
        solut_best    = Solution_init(data)

        do i = 0, Imax-1
            aux = data%rnd(data%rnd_index)
            data%rnd_index = data%rnd_index + 1

            alpha = R_table(aux)

            if (verbose) then
                write(*, '(A,I0)') '[+] Search ', i+1
                write(*, '(A)')    achar(9)//'[+] Constructing..'
            end if

            solut_crnt%s = construct(alpha, data)
            call update_subseq_info_matrix(solut_crnt, data)

            call Solution_cpy(solut_crnt, solut_partial, data)
            if (verbose) then
                write(*, '(A)') achar(9)//'[+] Looking for the best Neighbor..'
                write(*, '(A,F0.3)') achar(9)//'    Construction Cost: ', solut_partial%cost
            end if

            iterILS = 0
            do while (iterILS < Iils)
                call RVND(solut_crnt, data)
                if (solut_crnt%cost < solut_partial%cost - DBL_EPSILON_) then
                    call Solution_cpy(solut_crnt, solut_partial, data)
                    iterILS = 0
                end if

                solut_crnt%s = perturb(solut_partial, data)
                call update_subseq_info_matrix(solut_crnt, data)
                iterILS = iterILS + 1
            end do

            if (solut_partial%cost < solut_best%cost - DBL_EPSILON_) then
                call Solution_cpy(solut_partial, solut_best, data)
            end if

            if (verbose) then
                write(*, '(A,A)') achar(9)//'Current best cost: ', dbl_str(solut_best%cost)

                write(*, '(A)', advance='no') 'SOLUCAO: '
                block
                    integer :: n
                    do n = 0, size(solut_best%s)-1
                        write(*, '(I0,A)', advance='no') solut_best%s(n), ' '
                    end do
                end block
                write(*, '(A)') ''
            end if

        end do

        write(*, '(A,F0.2)') 'COST: ', solut_best%cost
    end subroutine GILS_RVND

    ! Formatacao equivalente ao operator<< padrao de C++ para double
    function dbl_str(x) result(str)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: str

        character(len=64) :: buf
        character(len=16) :: fmt
        character(len=:), allocatable :: mant, expo
        integer :: ex, nd, e

        if (x == 0.0_dp) then
            str = '0'
            return
        end if

        ex = floor(log10(abs(x)))

        if (ex < -4 .or. ex >= 6) then
            write(buf, '(ES15.5E2)') x
            str = trim(adjustl(buf))
            e = index(str, 'E')
            mant = trim_zeros(str(1:e-1))
            expo = str(e+1:)
            str = mant // 'e' // expo
        else
            nd = 5 - ex
            write(fmt, '(A,I0,A)') '(F0.', nd, ')'
            write(buf, fmt) x
            str = trim_zeros(trim(adjustl(buf)))
        end if

        if (str(1:1) == '.') then
            str = '0' // str
        else if (len(str) > 1) then
            if (str(1:2) == '-.') str = '-0' // str(2:)
        end if
    end function dbl_str

    function trim_zeros(s) result(str)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: str
        integer :: p

        str = s
        if (index(str, '.') == 0) return

        p = len(str)
        do while (p > 1 .and. str(p:p) == '0')
            p = p - 1
        end do
        if (str(p:p) == '.') p = p - 1
        str = str(1:p)
    end function trim_zeros

end module mlp_mod


program main
    use data_mod, only: dp, loadData
    use mlp_mod
    implicit none

    integer :: Imax
    integer :: Iils
    logical :: verbose
    integer :: i, argc
    character(len=256) :: arg

    type(tData) :: data
    type(tSolution) :: solut

    real(dp), allocatable :: cost(:,:)
    integer, allocatable :: rnd(:)
    integer :: rnd_count

    real(dp) :: t_start, t_end, cpu_time_used, res
    integer(8) :: c1, c2, crate
    integer(8) :: duration

    Imax = 10
    verbose = .false.

    argc = command_argument_count()
    do i = 1, argc
        call get_command_argument(i, arg)
        if (trim(arg) == '-v' .or. trim(arg) == '--verbose') then
            verbose = .true.
        end if
    end do

    data%dimen = loadData(cost, rnd, rnd_count)
    call move_alloc(cost, data%cost)
    allocate(data%rnd(0:rnd_count-1))
    data%rnd = rnd(0:rnd_count-1)
    data%rnd_index = 0

    solut = Solution_init(data)

    call random_seed()

    Iils = min(data%dimen, 100)

    call system_clock(c1, crate)
    call cpu_time(t_start)
    call GILS_RVND(Imax, Iils, data, verbose)
    call cpu_time(t_end)
    call system_clock(c2)

    cpu_time_used = t_end - t_start
    duration = (c2 - c1) * 1000_8 / crate

    res = real(duration, dp) / 10e2_dp

    write(*, '(A,A)') 'TIME: ', dbl_str(cpu_time_used)
    write(*, '(A,A)') 'wall clock (s): ', dbl_str(res)

    write(*, '(A,I0)') 'Tamanho RND ', rnd_count

    stop

end program main
