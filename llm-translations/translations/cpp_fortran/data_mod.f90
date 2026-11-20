! Traducao de Data.hpp

module data_mod
    use iso_fortran_env, only: real64, iostat_end, iostat_eor
    implicit none

    integer, parameter :: dp = real64

    private
    public :: dp, loadData, tDataClass

    ! Traducao da classe Data
    type :: tDataClass
        private
        integer :: dimen
        real(dp), allocatable :: cost(:,:)
        integer, allocatable :: rnd(:)
        integer :: rnd_size = 0
    contains
        procedure :: load     => Data_load
        procedure :: getRndVec => Data_getRndVec
        procedure :: getDimen  => Data_getDimen
        procedure :: getCostPtr => Data_getCostPtr
    end type tDataClass

contains

    ! double ** matrix_allocate(size_t sz)
    subroutine matrix_allocate(matrix, sz)
        real(dp), allocatable, intent(out) :: matrix(:,:)
        integer, intent(in) :: sz

        allocate(matrix(0:sz-1, 0:sz-1))
        matrix = 0.0_dp
    end subroutine matrix_allocate

    ! std::getline(file, line)
    subroutine getline(unit, line, ok)
        integer, intent(in) :: unit
        character(len=:), allocatable, intent(out) :: line
        logical, intent(out) :: ok

        character(len=4096) :: buffer
        integer :: sz, ios

        line = ''
        ok = .false.
        do
            read(unit, '(A)', advance='no', size=sz, iostat=ios) buffer
            if (sz > 0) line = line // buffer(1:sz)

            if (ios == iostat_eor) then
                ok = .true.
                return
            else if (ios == iostat_end) then
                ok = (len(line) > 0)
                return
            else if (ios /= 0) then
                ok = .false.
                return
            end if
        end do
    end subroutine getline

    ! push_back para o vetor rnd
    subroutine rnd_push_back(rnd, rnd_size, value)
        integer, allocatable, intent(inout) :: rnd(:)
        integer, intent(inout) :: rnd_size
        integer, intent(in) :: value

        integer, allocatable :: tmp(:)
        integer :: cap

        if (.not. allocated(rnd)) allocate(rnd(0:15))

        cap = size(rnd)
        if (rnd_size >= cap) then
            allocate(tmp(0:2*cap-1))
            tmp(0:rnd_size-1) = rnd(0:rnd_size-1)
            call move_alloc(tmp, rnd)
        end if

        rnd(rnd_size) = value
        rnd_size = rnd_size + 1
    end subroutine rnd_push_back

    ! int loadData(double *** matrix, vector<int> & rnd)
    function loadData(matrix, rnd, rnd_count) result(dimension)
        real(dp), allocatable, intent(out) :: matrix(:,:)
        integer, allocatable, intent(inout) :: rnd(:)
        integer, intent(out) :: rnd_count
        integer :: dimension

        integer :: unit, ios
        character(len=:), allocatable :: line
        logical :: ok, is_open, flag
        integer :: rnd_size, i, j, pos, ivalue
        real(dp) :: value

        dimension = 0
        rnd_count = 0

        open(newunit=unit, file='../distance_matrix', status='old', action='read', iostat=ios)
        is_open = (ios == 0)

        if (is_open) then

            call getline(unit, line, ok)
            read(line, *) dimension
            call matrix_allocate(matrix, dimension)

            flag = .false.
            rnd_size = -1
            i = 0
            do
                call getline(unit, line, ok)
                if (.not. ok) exit

                j = i + 1
                do while (index(line, ' ') /= 0)
                    pos = index(line, ' ')
                    read(line(1:pos-1), *) value
                    line = line(pos+1:)
                    matrix(i, j) = value
                    matrix(j, i) = value
                    j = j + 1
                end do

                i = i + 1

                if (flag) then
                    if (rnd_size < 0) then
                        read(line, *) rnd_size
                        if (.not. allocated(rnd)) allocate(rnd(0:max(rnd_size,1)-1))
                    else
                        read(line, *) ivalue
                        call rnd_push_back(rnd, rnd_count, ivalue)
                    end if
                end if

                if (line == 'RND') flag = .true.
            end do

        end if

        close(unit)
    end function loadData

    subroutine Data_load(this)
        class(tDataClass), intent(inout) :: this

        this%dimen = loadData(this%cost, this%rnd, this%rnd_size)
    end subroutine Data_load

    function Data_getRndVec(this) result(rnd)
        class(tDataClass), intent(in) :: this
        integer, allocatable :: rnd(:)

        allocate(rnd(0:this%rnd_size-1))
        rnd = this%rnd(0:this%rnd_size-1)
    end function Data_getRndVec

    function Data_getDimen(this) result(dimen)
        class(tDataClass), intent(in) :: this
        integer :: dimen

        dimen = this%dimen
    end function Data_getDimen

    function Data_getCostPtr(this) result(cost)
        class(tDataClass), intent(in) :: this
        real(dp), allocatable :: cost(:,:)

        cost = this%cost
    end function Data_getCostPtr

end module data_mod
