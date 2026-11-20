! Traducao de readData.cpp / readData.h

module read_data_mod
    use iso_fortran_env, only: real64, iostat_end, iostat_eor
    implicit none

    integer, parameter, private :: dp = real64

    private
    public :: readData

contains

    ! equivalente a "in >> token" (leitura delimitada por espacos em branco)
    subroutine next_token(unit, tok, ok)
        integer, intent(in) :: unit
        character(len=:), allocatable, intent(out) :: tok
        logical, intent(out) :: ok

        character(len=1) :: c
        integer :: ios

        tok = ''
        ok = .false.

        do
            read(unit, '(A1)', advance='no', iostat=ios) c
            if (ios == iostat_end) return
            if (ios == iostat_eor) cycle
            if (ios /= 0) return
            if (c /= ' ' .and. c /= achar(9) .and. c /= achar(13)) exit
        end do

        tok = c
        do
            read(unit, '(A1)', advance='no', iostat=ios) c
            if (ios /= 0) exit
            if (c == ' ' .or. c == achar(9) .or. c == achar(13)) exit
            tok = tok // c
        end do

        ok = .true.
    end subroutine next_token

    ! equivalente a "in >> inteiro"
    subroutine next_int(unit, value, ok)
        integer, intent(in) :: unit
        integer, intent(out) :: value
        logical, intent(out) :: ok

        character(len=:), allocatable :: tok

        value = 0
        call next_token(unit, tok, ok)
        if (ok) read(tok, *) value
    end subroutine next_int

    ! equivalente a "in >> double"
    subroutine next_real(unit, value, ok)
        integer, intent(in) :: unit
        real(dp), intent(out) :: value
        logical, intent(out) :: ok

        character(len=:), allocatable :: tok

        value = 0.0_dp
        call next_token(unit, tok, ok)
        if (ok) read(tok, *) value
    end subroutine next_real

    ! void readData( int argc, char** argv, int* Dimension, double ***Mdist )
    subroutine readData(argc, Dimension, Mdist)
        integer, intent(in) :: argc
        integer, intent(out) :: Dimension
        real(dp), allocatable, intent(out) :: Mdist(:,:)

        integer :: N, ios, unit, i, j, tempCity
        character(len=:), allocatable :: arquivo, ewt, ewf
        character(len=1024) :: instancia
        real(dp), allocatable :: x(:), y(:), dist(:,:)
        real(dp), allocatable :: latitude(:), longitude(:)
        integer, allocatable :: tempX(:), tempY(:)
        logical :: ok

        if (argc < 2) then
            write(*, '(A)') ''
            write(*, '(A)') 'Faltando parametros'
            write(*, '(A)') ' ./exec [Instancia] '
            stop 1
        end if

        if (argc > 2) then
            write(*, '(A)') ''
            write(*, '(A)') 'Muitos parametros'
            write(*, '(A)') ' ./exec [Instancia] '
            stop 1
        end if

        call get_command_argument(1, instancia)

        open(newunit=unit, file=trim(instancia), status='old', action='read', iostat=ios)

        if (ios /= 0) then
            write(*, '(A)') 'arquivo nao pode ser aberto'
            stop 1
        end if

        arquivo = ''
        do while (arquivo /= 'DIMENSION:' .and. arquivo /= 'DIMENSION')
            call next_token(unit, arquivo, ok)
            if (.not. ok) exit
        end do

        if (arquivo == 'DIMENSION') call next_token(unit, arquivo, ok)

        call next_int(unit, N, ok)

        do while (arquivo /= 'EDGE_WEIGHT_TYPE:' .and. arquivo /= 'EDGE_WEIGHT_TYPE')
            call next_token(unit, arquivo, ok)
            if (.not. ok) exit
        end do
        if (arquivo == 'EDGE_WEIGHT_TYPE') call next_token(unit, arquivo, ok)

        call next_token(unit, ewt, ok)

        allocate(x(0:N), y(0:N))
        x = 0.0_dp
        y = 0.0_dp

        ! Alocar matriz 2D
        allocate(dist(0:N, 0:N))
        dist = 0.0_dp

        if (ewt == 'EXPLICIT') then

            do while (arquivo /= 'EDGE_WEIGHT_FORMAT:' .and. arquivo /= 'EDGE_WEIGHT_FORMAT')
                call next_token(unit, arquivo, ok)
                if (.not. ok) exit
            end do

            if (arquivo == 'EDGE_WEIGHT_FORMAT') call next_token(unit, arquivo, ok)
            call next_token(unit, ewf, ok)

            if (ewf == 'FUNCTION') then
                write(*, '(A)') 'FUNCTION - Nao suportado!'

            else if (ewf == 'FULL_MATRIX') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                do i = 1, N
                    do j = 1, N
                        call next_real(unit, dist(i,j), ok)
                    end do
                end do

            else if (ewf == 'UPPER_ROW') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                do i = 1, N-1
                    do j = i+1, N
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                    end do
                end do

                do i = 1, N
                    dist(i,i) = 0.0_dp
                end do

            else if (ewf == 'LOWER_ROW') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                do i = 2, N
                    do j = 1, i-1
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                    end do
                end do

                do i = 1, N
                    dist(i,i) = 0.0_dp
                end do

            else if (ewf == 'UPPER_DIAG_ROW') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                do i = 1, N
                    do j = i, N
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                    end do
                end do

            else if (ewf == 'LOWER_DIAG_ROW') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                do i = 1, N
                    do j = 1, i
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                    end do
                end do

            else if (ewf == 'UPPER_COL') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                do j = 2, N
                    do i = 1, j-1
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                    end do
                end do

                do i = 1, N
                    dist(i,i) = 0.0_dp
                end do

            else if (ewf == 'LOWER_COL') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                ! obs: o laco interno do codigo original incrementa j (e nao i)
                j = 1
                do while (j < N)
                    i = j + 1
                    do while (i < N+1)
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                        j = j + 1
                    end do
                    j = j + 1
                end do

                do i = 1, N
                    dist(i,i) = 0.0_dp
                end do

            else if (ewf == 'UPPER_DIAG_COL') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                do j = 1, N
                    do i = 1, j
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                    end do
                end do

            else if (ewf == 'LOWER_DIAG_COL') then

                do while (arquivo /= 'EDGE_WEIGHT_SECTION')
                    call next_token(unit, arquivo, ok)
                    if (.not. ok) exit
                end do

                ! Preencher Matriz Distancia
                ! obs: o laco interno do codigo original incrementa j (e nao i)
                j = 1
                do while (j < N+1)
                    i = j
                    do while (i < N+1)
                        call next_real(unit, dist(i,j), ok)
                        dist(j,i) = dist(i,j)
                        j = j + 1
                    end do
                    j = j + 1
                end do

            end if

        else if (ewt == 'EUC_2D') then

            do while (arquivo /= 'NODE_COORD_SECTION')
                call next_token(unit, arquivo, ok)
                if (.not. ok) exit
            end do

            ! ler coordenadas
            do i = 1, N
                call next_int(unit, tempCity, ok)
                call next_real(unit, x(i), ok)
                call next_real(unit, y(i), ok)
            end do

            ! Calcular Matriz Distancia (Euclidiana)
            do i = 1, N
                do j = 1, N
                    dist(i,j) = floor(CalcDistEuc(x, y, i, j) + 0.5_dp)
                end do
            end do

        else if (ewt == 'EUD_3D') then
            write(*, '(A)') 'EUC_3D - Nao suportado!'

        else if (ewt == 'MAX_2D') then
            write(*, '(A)') 'MAX_2D - Nao suportado!'

        else if (ewt == 'MAX_3D') then
            write(*, '(A)') 'MAX_3D - Nao Suportado!'

        else if (ewt == 'MAN_2D') then
            write(*, '(A)') 'MAN_2D - Nao suportado!'

        else if (ewt == 'MAN_3D') then
            write(*, '(A)') 'MAN_3D - Nao Suportado!'

        else if (ewt == 'CEIL_2D') then

            do while (arquivo /= 'NODE_COORD_SECTION')
                call next_token(unit, arquivo, ok)
                if (.not. ok) exit
            end do

            ! ler coordenadas
            do i = 1, N
                call next_int(unit, tempCity, ok)
                call next_real(unit, x(i), ok)
                call next_real(unit, y(i), ok)
            end do

            ! Calcular Matriz Distancia (Euclidiana)
            do i = 1, N
                do j = 1, N
                    dist(i,j) = ceiling(CalcDistEuc(x, y, i, j))
                end do
            end do

        else if (ewt == 'GEO') then

            do while (arquivo /= 'NODE_COORD_SECTION')
                call next_token(unit, arquivo, ok)
                if (.not. ok) exit
            end do

            ! ler coordenadas
            do i = 1, N
                call next_int(unit, tempCity, ok)
                call next_real(unit, x(i), ok)
                call next_real(unit, y(i), ok)
            end do

            allocate(latitude(0:N), longitude(0:N))

            call CalcLatLong(x, y, N, latitude, longitude)

            ! Calcular Matriz Distancia
            do i = 1, N
                do j = 1, N
                    dist(i,j) = CalcDistGeo(latitude, longitude, i, j)
                end do
            end do

        else if (ewt == 'ATT') then

            do while (arquivo /= 'NODE_COORD_SECTION')
                call next_token(unit, arquivo, ok)
                if (.not. ok) exit
            end do

            ! ler coordenadas
            allocate(tempX(0:N), tempY(0:N))

            do i = 1, N
                call next_int(unit, tempCity, ok)
                call next_int(unit, tempX(i), ok)
                call next_int(unit, tempY(i), ok)
                x(i) = tempX(i)
                y(i) = tempY(i)
            end do

            ! Calcular Matriz Distancia (Pseudo-Euclidiana)
            do i = 1, N
                do j = 1, N
                    dist(i,j) = CalcDistAtt(x, y, i, j)
                end do
            end do

        else if (ewt == 'XRAY1') then
            write(*, '(A)') 'XRAY1 - Nao suportado!'

        else if (ewt == 'XRAY2') then
            write(*, '(A)') 'XRAY2 - Nao suportado!'

        else if (ewt == 'SPECIAL') then
            write(*, '(A)') 'SPECIAL - Nao suportado!'

        end if

        close(unit)

        Dimension = N
        call move_alloc(dist, Mdist)
    end subroutine readData

    function CalcDistEuc(X, Y, I, J) result(d)
        real(dp), intent(in) :: X(0:), Y(0:)
        integer, intent(in) :: I, J
        real(dp) :: d

        d = sqrt((X(I) - X(J))**2 + (Y(I) - Y(J))**2)
    end function CalcDistEuc

    function CalcDistAtt(X, Y, I, J) result(dij)
        real(dp), intent(in) :: X(0:), Y(0:)
        integer, intent(in) :: I, J
        real(dp) :: dij

        ! Calcula Pseudo Distancia Euclidiana
        real(dp) :: rij, tij

        rij = sqrt(((X(I) - X(J))**2 + (Y(I) - Y(J))**2) / 10.0_dp)
        tij = floor(rij + 0.5_dp)

        if (tij < rij) then
            dij = tij + 1.0_dp
        else
            dij = tij
        end if
    end function CalcDistAtt

    subroutine CalcLatLong(X, Y, n, latit, longit)
        real(dp), intent(in) :: X(0:), Y(0:)
        integer, intent(in) :: n
        real(dp), intent(inout) :: latit(0:), longit(0:)

        real(dp), parameter :: PI = 3.141592_dp
        real(dp) :: min_
        integer :: deg, i

        do i = 1, n
            deg = int(X(i))
            min_ = X(i) - deg
            latit(i) = PI * (deg + 5.0_dp * min_ / 3.0_dp) / 180.0_dp
        end do

        do i = 1, n
            deg = int(Y(i))
            min_ = Y(i) - deg
            longit(i) = PI * (deg + 5.0_dp * min_ / 3.0_dp) / 180.0_dp
        end do
    end subroutine CalcLatLong

    function CalcDistGeo(latit, longit, I, J) result(d)
        real(dp), intent(in) :: latit(0:), longit(0:)
        integer, intent(in) :: I, J
        real(dp) :: d

        real(dp) :: q1, q2, q3
        real(dp), parameter :: RRR = 6378.388_dp

        q1 = cos(longit(I) - longit(J))
        q2 = cos(latit(I) - latit(J))
        q3 = cos(latit(I) + latit(J))

        d = real(int(RRR * acos(0.5_dp*((1.0_dp+q1)*q2 - (1.0_dp-q1)*q3)) + 1.0_dp), dp)
    end function CalcDistGeo

end module read_data_mod
