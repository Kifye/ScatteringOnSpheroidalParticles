! Contains types and functions to streamline different calculations

module communication
    use angle
    use geometry
    use logging
    use spheroidal_scatterer
    implicit none
!    private

    type, public :: ScatteringFactors
        real(knd) :: ext
        real(knd) :: sca
        !  one of values qfactors, cfactors, normalized_cfactors
        !  from enum in module constants
        integer :: factors_type
    contains
        procedure, public :: set_to_zeroes
        procedure, public :: get_c_factors_from_q
        procedure, public :: get_normalized_c_factors_from_q
        procedure, public :: get_normalized_c_factors_from_c
        procedure, public :: log => log_factors
    end type ScatteringFactors

    !  represents what is required to calculate for a mode TE or TM
    type :: ModeQuery
        !  following parameters are optional and are in order
        !  - calculate_uv_factors
        !  - calculate_pq_factors
        !  - return_uv_tmatrix 
        !  - return_pq_tmatrix
        !  - return_uv_solution
        !  - return_pq_solution
        !  if they are not present, they are assumed to be false
 
        !  should the result contain tmatrices for the
        !  respective potentials, only works if
        !  the factos for those potentials are required
        !  their sizes are
        ! PQ: (1:matrix_size, 1:matrix_size)
        ! UV: (1:2*matrix_size, 1:2*matrix_size, minm:actual_maxm)
        !  where maxtrix_size, minm are defined in SpheroidalQuery
        !  and actual_maxm is obtained during calculation
        !  from required accuracy and is <= maxm from SpheroidalQuery
        logical :: return_pq_tmatrix, return_uv_tmatrix
        !  should the result contain solution for the
        !  respective potentials
        !  their sizes are
        ! PQ: (1:matrix_size)
        ! UV: (1:2*matrix_size, minm:actual_maxm)
        logical :: return_pq_solution, return_uv_solution
        !  should the factors from PQ potentials be calculated for
        !  TE and TM modes
        logical :: calculate_pq_factors
        !  should the factors from UV potentials be calculated for
        logical :: calculate_uv_factors
    contains
        procedure, public :: set => set_mode_query
        procedure, public :: log => log_mode_query
        procedure, public :: do_not_return_arrays
    end type ModeQuery
    
    !  contains the calculation result for a specific mode
    !  the fields correspond to the logical values in ModeQuery 
    type, public :: ModeResult
        type(ScatteringFactors) :: uv_factors
        type(ScatteringFactors) :: pq_factors

        !  if matrix_size was set in the query
        !  [1:matrix_size]
        complex(knd), allocatable, dimension(:) :: pq_solution
        !  [1:2*matrix_size][minm:maxm]
        complex(knd), allocatable, dimension(:,:) :: uv_solution
        !  [1:matrix_size][1:matrix_size]
        complex(knd), allocatable, dimension(:,:) :: pq_tmatrix
        !  [1:2*matrix_size][1:2*matrix_size][minm:maxm]
        complex(knd), allocatable, dimension(:,:,:) :: uv_tmatrix

        type(ModeQuery), pointer :: query
    contains
        procedure, public :: log => log_mode_result
        procedure, private :: set_mode_result
        final :: delete_mode_result
    end type ModeResult    

    !  this type represents the query to the calculating unit
    !  the object of this type can be created with create_query function
    !  and should be passed to .... to calculate the required parameters
    type, public :: ScatteringQuery
        !  size of the matrices(Nmax) to use
        !  the T-matrix for the axisymmetric mode is matrix_size x matrix_size
        !  for the nonaxisymmetric mode 2matrix_size x 2matrix_size
        !  if matrix_size <= 0 the matrix size will be iterated over until
        !  the required accuracy is reached or it reaches the value
        !  MAXIMUM_MATRIX_SIZE from module constants
        integer :: matrix_size
        integer :: spherical_matrix_size
        !  factors for the nonaxisymmetric mode
        !  are calculated for m in [minm:maxm]
        integer :: minm, maxm
        !  required accuracy
        real(knd) :: accuracy

        !  Query body
        !  what factors for what mode should be calculated 
        !  and whether tmatrices or solution need to be returned
        type(ModeQuery) :: query_te
        type(ModeQuery) :: query_spherical_te
        type(ModeQuery) :: query_spherical_orthogonal_te
        type(ModeQuery) :: query_tm
        type(ModeQuery) :: query_spherical_tm
        type(ModeQuery) :: query_spherical_orthogonal_tm
        !  should the scattering indicatrix cut be calculated
        !  if set true then the values uv_te, uv_tm
        !  are ignored and assumed to be = .true.
        logical :: calculate_scattering_indicatrix
        !  array of directions of where the indicatrix should be calculated
        type(SphericalDirection), allocatable, dimension(:) :: indicatrix_points
    contains
        procedure, public :: log => log_query
        final :: delete_scattering_query
    end type ScatteringQuery

    public create_query
    !  this type represents the result of the calculation
    !  all fields except for accuracy_achieved
    !  are ..._value fields which correspond to the
    !  ..._query fields from the Query
    !  accuracy_achieved indicated wether the accuracy from the query
    !  was achieved during calculations
    type, public :: ScatteringResult
        type(ModeResult) :: result_te
        type(ModeResult) :: result_spherical_te
        type(ModeResult) :: result_spherical_orthogonal_te
        type(ModeResult) :: result_tm
        type(ModeResult) :: result_spherical_tm
        type(ModeResult) :: result_spherical_orthogonal_tm
        !  [1:4][1:4][1:size(query.indicatrix_points)]
        real(knd), allocatable, dimension(:,:,:) :: indicatrix
        !  is true if the internal checks passed an accuracy threshold
        !  given in query
        !  checks:
        !  - inside VB function check
        !  - m accumulation check
        logical :: accuracy_achieved
        !  if the accuracy was not achieved during calculation
        !  reason can be set to show why
        !  it can be on of values (in the order of priority):
        !  1. VB_functions
        !  2. no_m_convergeance
        integer :: reason
        integer :: actual_maxm

        type(ScatteringQuery), pointer :: query
    contains
        procedure, public :: log => log_result
        final :: delete_scattering_result
    end type ScatteringResult

    public create_result
    public :: read_input, read_long_input

contains
    !  ScatteringFactors
    subroutine set_to_zeroes(factors)
        class(ScatteringFactors) :: factors

        factors%ext = 0
        factors%sca = 0

    end subroutine set_to_zeroes

    subroutine log_factors(factors, level, fd)
        class(ScatteringFactors) :: factors
        integer :: level, fd

        if (.not. LOGS_TO_WRITE(level)) then
            return
        end if

!        if (factors%factors_type == qfactors) then
!            write(fd, *) 'Q factors values:'
!        elseif (factors%factors_type == cfactors) then
!            write(fd, *) 'Cross-sections values:'
!        elseif (factors%factors_type == normalized_cfactors) then
!            write(fd, *) 'Normalized cross-sections values:'
!        else
!            write(fd, *) 'Unrecognized factors type:'
!        end if

        write(fd, *) 'ext = ', factors%ext
        write(fd, *) 'sca = ', factors%sca
        write(fd, *) 'abs = ', abs(factors%ext) - factors%sca

    end subroutine log_factors

    real(knd) function convertQtoC(Q, f, rv, ab, alpha) result(C)
        integer, intent(in) :: f
        real(knd), intent(in) :: Q, rv, ab, alpha

        real(knd) :: a, b

        a = get_a(f, rv, ab)
        b = get_b(f, rv, ab)

        if (f == 1) then
            C = Q * PI * b * sqrt(a**2 * sin(alpha) ** 2 + b ** 2 * cos(alpha) ** 2)
        else
            C = Q * PI * a * sqrt(b**2 * sin(alpha) ** 2 + a ** 2 * cos(alpha) ** 2)
        end if
    end function convertQtoC

    type(ScatteringFactors) function get_c_factors_from_q(factors, shape)
        class(ScatteringFactors), intent(in) :: factors
        type(SpheroidalShape), intent(in) :: shape
!        if (factors%factors_type == normalized_cfactors) then
!            get_c_factors%ext = factors%ext * PI * shape%rv ** 2
!            get_c_factors%sca = factors%sca * PI * shape%rv ** 2
!        elseif (factors%factors_type == qfactors) then
            get_c_factors_from_q%ext = convertQtoC(factors%ext, shape%spheroidal_type, shape%rv, &
                    shape%ab, shape%alpha%value)
            get_c_factors_from_q%sca = convertQtoC(factors%sca, shape%spheroidal_type, shape%rv, &
                    shape%ab, shape%alpha%value)
!        end if
    end function get_c_factors_from_q

    type(ScatteringFactors) function get_normalized_c_factors_from_q(factors, shape) result(result)
        class(ScatteringFactors), intent(in) :: factors
        type(SpheroidalShape), intent(in) :: shape

        result = factors
!        if (factors%factors_type == cfactors .or. factors%factors_type == cfactors) then
!            if (factors%factors_type == qfactors) then
                result = factors%get_c_factors_from_q(shape)
!            end if
            result%ext = result%ext / (PI * shape%rv ** 2)
            result%sca = result%sca / (PI * shape%rv ** 2)
!        end if
    end function get_normalized_c_factors_from_q

    type(ScatteringFactors) function get_normalized_c_factors_from_c(factors, shape) result(result)
        class(ScatteringFactors), intent(in) :: factors
        type(SpheroidalShape), intent(in) :: shape

        result = factors
        !        if (factors%factors_type == cfactors .or. factors%factors_type == cfactors) then
        !            if (factors%factors_type == qfactors) then
        !            end if
        result%ext = result%ext / (PI * shape%rv ** 2)
        result%sca = result%sca / (PI * shape%rv ** 2)
        !        end if
    end function get_normalized_c_factors_from_c


    ! Mode Query
    subroutine set_mode_query(query, calculate_uv_factors, calculate_pq_factors, &
            return_uv_tmatrix, return_pq_tmatrix, return_uv_solution, return_pq_solution)
        
        class(ModeQuery), intent(out) :: query
        logical, optional, intent(in) :: calculate_uv_factors, calculate_pq_factors, &
                return_uv_tmatrix, return_pq_tmatrix, return_uv_solution, return_pq_solution

        !  factors
        if (present(calculate_uv_factors)) then
            query%calculate_uv_factors = calculate_uv_factors
        else
            query%calculate_uv_factors = .false.
        end if
        if (present(calculate_pq_factors)) then
            query%calculate_pq_factors = calculate_pq_factors
        else
            query%calculate_pq_factors = .false.
        end if
        !  tmatrix
        if (present(return_uv_tmatrix)) then
            query%return_uv_tmatrix = return_uv_tmatrix
        else
            query%return_uv_tmatrix = .false.
        end if
        if (present(return_pq_tmatrix)) then
            query%return_pq_tmatrix = return_pq_tmatrix
        else
            query%return_pq_tmatrix = .false.
        end if
        !  solution
        if (present(return_uv_solution)) then
            query%return_uv_solution = return_uv_solution
        else
            query%return_uv_solution = .false.
        end if
        if (present(return_pq_solution)) then
            query%return_pq_solution = return_pq_solution
        else
            query%return_pq_solution = .false.
        end if

    end subroutine set_mode_query

    subroutine log_mode_query(query, level, fd)
        class(ModeQuery) :: query
        integer :: level, fd

        if (.not. LOGS_TO_WRITE(level) .or. &
                (.not.query%calculate_uv_factors .and. &
                        .not.query%calculate_pq_factors .and. &
                                .not.query%return_pq_tmatrix .and. &
                                        .not.query%return_uv_tmatrix .and. &
                                                .not.query%return_pq_solution .and. &
                                                        .not.query%return_uv_solution)) then
            return
        end if
        if (query%calculate_uv_factors) then
            write(fd, *) 'calculate UV factors'
        end if
        if (query%calculate_pq_factors) then
            write(fd, *) 'calculate PQ factors'
        end if
        if (query%return_uv_tmatrix) then
            write(fd, *) 'return UV T-matrix'
        end if
        if (query%return_pq_tmatrix) then
            write(fd, *) 'return PQ T-matrix'
        end if
        if (query%return_uv_solution) then
            write(fd, *) 'return UV solution'
        end if
        if (query%return_pq_solution) then
            write(fd, *) 'return PQ solution'
        end if
        call log_new_line(level, fd)

    end subroutine log_mode_query

    subroutine set_mode_result(result, query, matrix_size, minm, maxm)
        class(ModeResult), intent(out) :: result
        type(ModeQuery), target, intent(in) :: query
        integer, intent(in) :: matrix_size, minm, maxm

        result%query => query
        call result%uv_factors%set_to_zeroes()
        call result%pq_factors%set_to_zeroes()

        if (query%return_uv_tmatrix) then
            allocate(result%uv_tmatrix(2 * matrix_size, 2 * matrix_size, minm:maxm))
            result%uv_tmatrix = 0
        end if
        if (query%return_pq_tmatrix) then
            allocate(result%pq_tmatrix(matrix_size, matrix_size))
            result%pq_tmatrix = 0
        end if
        if (query%return_uv_solution) then
            allocate(result%uv_solution(2 * matrix_size, minm:maxm))
            result%uv_solution = 0
        end if
        if (query%return_pq_solution) then
            allocate(result%pq_solution(matrix_size))
            result%pq_solution = 0
        end if
    end subroutine set_mode_result

    subroutine do_not_return_arrays(query)
        class(ModeQuery), intent(inout) :: query

        query%return_uv_tmatrix = .false.
        query%return_pq_tmatrix = .false.
        query%return_uv_solution = .false.
        query%return_pq_solution = .false.

    end subroutine do_not_return_arrays

    !  ModeResult
    subroutine log_mode_result(result, level, fd, border)
        class(ModeResult), intent(in) :: result
        integer, intent(in) :: level, fd
        integer, optional, intent(in) :: border

        integer :: length

        if (.not. LOGS_TO_WRITE(level)) then
            return
        end if

        if (present(border)) then
            length = border
        end if
        if (result%query%calculate_uv_factors) then
            write(fd, *) 'UV factors:'
            call result%uv_factors%log(level, fd)
        end if
        if (result%query%calculate_pq_factors) then
            write(fd, *) 'PQ factors:'
            call result%pq_factors%log(level, fd)
        end if
        if (result%query%return_uv_tmatrix .and. allocated(result%uv_tmatrix)) then
            !call log_matrix(fd, 'UV T-matrix', result%uv_tmatrix, .false., border)
        end if
        if (result%query%return_pq_tmatrix .and. allocated(result%pq_tmatrix)) then
            call log_matrix(fd, 'PQ T-matrix', result%pq_tmatrix, .false., border)
        end if
        if (result%query%return_uv_solution.and. allocated(result%uv_solution)) then
            call log_matrix(fd, 'UV solution', result%uv_solution,.false., border)
        end if
        if (result%query%return_pq_solution.and. allocated(result%pq_solution)) then
            call log_array(fd, 'PQ solution', result%pq_solution, border)
        end if
    end subroutine log_mode_result

    subroutine delete_mode_result(result)
        type(ModeResult), intent(inout) :: result

        if (allocated(result%uv_tmatrix)) then
            deallocate(result%uv_tmatrix)
        end if
        if (allocated(result%pq_tmatrix)) then
            deallocate(result%pq_tmatrix)
        end if
        if (allocated(result%uv_solution)) then
            deallocate(result%uv_solution)
        end if
        if (allocated(result%pq_solution)) then
            deallocate(result%pq_solution)
        end if
    end subroutine delete_mode_result

    !  ScatteringQuery
    subroutine log_query(query, level, fd)
        class(ScatteringQuery) :: query
        integer :: level, fd

        if (.not. LOGS_TO_WRITE(level)) then
            return
        end if

        write(fd, *) 'Query:'
        write(fd, *) 'matrix_size = ', query%matrix_size
        write(fd, *) 'minm = ', query%minm, 'maxm = ', query%maxm
        write(fd, *) 'expected accuracy ', query%accuracy
        write(fd, *) 'Required:'
        write(fd, *) ' TE:'
        call query%query_te%log(level, fd)
        write(fd, *) ' TM:'
        call query%query_tm%log(level, fd)
        if (query%calculate_scattering_indicatrix) then
            write(fd, *) 'Indicatrix for ', size(query%indicatrix_points), ' points'
        end if
        call log_new_line(level, fd)

    end subroutine log_query

    subroutine delete_scattering_query(query)
        type(ScatteringQuery), intent(inout) :: query

        if (allocated(query%indicatrix_points)) then
            deallocate(query%indicatrix_points)
        end if

    end subroutine delete_scattering_query

    function create_query(matrix_size, minm, maxm, accuracy, &
            uv_te, uv_tm, pq_te, pq_tm, need_indicatrix, indicatrix_points, &
            return_tmatrix, return_solution, double_for_spherical, spherical_matrix_size) result(this)

        type(ScatteringQuery) :: this
        logical, intent(in), optional :: uv_te, uv_tm, pq_te, pq_tm, return_tmatrix, return_solution, need_indicatrix, &
                double_for_spherical
        integer, intent(in)  :: matrix_size, minm, maxm
        integer, optional :: spherical_matrix_size
        type(SphericalDirection), dimension(:), intent(in), optional :: indicatrix_points
        real(knd), intent(in)  :: accuracy

        this%matrix_size = matrix_size
        this%minm = minm
        this%maxm = maxm
        this%accuracy = accuracy
        
        call this%query_te%set(uv_te, pq_te, return_tmatrix, return_tmatrix, return_solution, return_solution)
        call this%query_tm%set(uv_tm, pq_tm, return_tmatrix, return_tmatrix, return_solution, return_solution)
        call this%query_spherical_te%set()
        call this%query_spherical_tm%set()
        if (present(double_for_spherical)) then
            if (double_for_spherical) then
                call this%query_spherical_te%set(uv_te, pq_te, return_tmatrix, return_tmatrix, return_solution, return_solution)
                call this%query_spherical_orthogonal_te%set(uv_te, pq_te, return_tmatrix, return_tmatrix, &
                        return_solution, return_solution)
                call this%query_spherical_tm%set(uv_tm, pq_tm, return_tmatrix, return_tmatrix, &
                        return_solution, return_solution)
                call this%query_spherical_orthogonal_tm%set(uv_tm, pq_tm, return_tmatrix, return_tmatrix, &
                        return_solution, return_solution)
            end if
        end if
        this%spherical_matrix_size = this%matrix_size
        if (present(spherical_matrix_size)) then
            this%spherical_matrix_size = spherical_matrix_size
        end if
        this%calculate_scattering_indicatrix = .false.
        if (present(need_indicatrix)) then
            if (need_indicatrix) then
                if (allocated(this%indicatrix_points)) then
                    deallocate(this%indicatrix_points)
                end if
                allocate(this%indicatrix_points(size(indicatrix_points)))
                this%indicatrix_points = indicatrix_points
            end if
        end if
        
        call log_message('Created query', DETAIL, FILE_DESCRIPTOR(DETAIL))
        call log_query(this, DETAIL, FILE_DESCRIPTOR(DETAIL))
    end function create_query

    !  ScatteringResult
    subroutine log_result(result, level, fd)
        class(ScatteringResult) :: result
        integer :: level, fd, i

        if (.not. LOGS_TO_WRITE(level)) then
            return
        end if
        write(fd, *) 'Result:'
        write(fd, *) 'actual maxm = ', result%actual_maxm
        write(fd, *) 'TE:'
        call result%result_te%log(level, fd)
        write(fd,*) 'spherical'
        call result%result_spherical_te%log(level, fd)
        write(fd,*) 'orthogonal'
        call result%result_spherical_orthogonal_te%log(level, fd)
        write(fd, *) 'TM:'
        call result%result_tm%log(level, fd)
        write(fd,*) 'spherical'
        call result%result_spherical_tm%log(level, fd)
        write(fd,*) 'orthogonal'
        call result%result_spherical_orthogonal_tm%log(level, fd)
        if (level == DETAIL .and. result%query%calculate_scattering_indicatrix .and. allocated(result%indicatrix)) then
            do i = 1, size(result%indicatrix)
                write(fd, *) result%query%indicatrix_points(i)%theta%in_degrees(), &
                        result%query%indicatrix_points(i)%phi%in_degrees(), &
                        result%indicatrix(:,:,i)
            end do
        end if
        if (result%accuracy_achieved) then
            write(fd, *) 'accuracy achieved'
        else
            write(fd, *) 'accuracy not achieved, reason =', result%reason
        end if
        call log_new_line(level, fd)

    end subroutine log_result

    subroutine delete_scattering_result(result)
        type(ScatteringResult), intent(inout) :: result

        if (allocated(result%indicatrix)) then
            deallocate(result%indicatrix)
        end if
    end subroutine delete_scattering_result

    type(ScatteringResult) function create_result(query) result(result)
        type(ScatteringQuery), target, intent(in) :: query

        result%query => query
        call result%result_te%set_mode_result(query%query_te, query%matrix_size, query%minm, query%maxm)
        call result%result_spherical_te%set_mode_result(query%query_spherical_te, query%spherical_matrix_size, &
                query%minm, query%maxm)
        call result%result_spherical_orthogonal_te%set_mode_result(query%query_spherical_orthogonal_te, &
                query%spherical_matrix_size, &
                query%minm, query%maxm)
        call result%result_tm%set_mode_result(query%query_tm, query%matrix_size, query%minm, query%maxm)
        call result%result_spherical_tm%set_mode_result(query%query_spherical_tm, query%spherical_matrix_size, &
                query%minm, query%maxm)
        call result%result_spherical_orthogonal_tm%set_mode_result(query%query_spherical_orthogonal_tm, &
                query%spherical_matrix_size, &
                query%minm, query%maxm)

        if (query%calculate_scattering_indicatrix) then
            if (allocated(result%indicatrix)) then
                deallocate(result%indicatrix)
            end if
            allocate(result%indicatrix(4,4,size(query%indicatrix_points)))
            result%indicatrix = 0
        end if
    end function create_result

    !  Reads input values from the file filename
    !  all values should be written in the form "key = value"
    !  in different lines
    !  out of rv, xv only one can be present, the other shall be
    !  calculated from it
    !  if both are present, xv is calculated from rv and the
    !  given xv is ignored
    !  if values are not present they are set to default:
    !    f = 1 (prolate)
    !    rv = 1
    !    xv = 1
    !    ab = 2
    !    alpha = 0
    !    lambda = 2Pi
    !    ri = 1.5 + 0.0i
    subroutine read_input(filename, f, nol, rv, xv, ab, alpha, lambda, ri, matrix_size, minm, maxm)
        character(*), intent(in) :: filename
        integer, intent(inout) :: f, nol, matrix_size, minm, maxm
        real(knd), allocatable, intent(inout) :: rv(:), xv(:), ab(:)
        real(knd), intent(inout) :: alpha, lambda
        complex(knd), allocatable, intent(inout) :: ri(:)

        character(1024) :: line, name, eq
        integer :: ios, i
        logical :: found_rv, found_xv

        open(FD_GENERAL, file=filename)
        read(FD_GENERAL, *) f
        read(FD_GENERAL, *) nol
        allocate(rv(nol), xv(nol), ab(nol), ri(0:nol))
        read(FD_GENERAL, *) xv
        read(FD_GENERAL, *) ab
        read(FD_GENERAL, *) ri
        read(FD_GENERAL, *) lambda
        read(FD_GENERAL, *) alpha
        read(FD_GENERAL, *) matrix_size
        read(FD_GENERAL, *) minm
        read(FD_GENERAL, *) maxm

        close(FD_GENERAL)

        if (matrix_size == 0) then
            matrix_size = size_of_matrices(f, xv(1), ab(1), ri(1))
        end if
        alpha = alpha / 180_knd * PI
        alpha = min(alpha, PI / 2q0 - 1q-24)
        rv = xv * lambda / (2.0_knd * PI)

        write(FILE_DESCRIPTOR(INFO), *) 'Read input:'
        write(FILE_DESCRIPTOR(INFO), *) 'f = ', f
        write(FILE_DESCRIPTOR(INFO), *) 'rv = ', rv
        write(FILE_DESCRIPTOR(INFO), *) 'xv = ', xv
        write(FILE_DESCRIPTOR(INFO), *) 'ab = ', ab
        write(FILE_DESCRIPTOR(INFO), *) 'alpha = ', alpha, 'radian = ', alpha / PI * 180_knd, 'degrees'
        write(FILE_DESCRIPTOR(INFO), *) 'lambda = ', lambda
        write(FILE_DESCRIPTOR(INFO), *) 'ri = ', ri
        write(FILE_DESCRIPTOR(INFO), *) 'lnum = ', matrix_size
        write(FILE_DESCRIPTOR(INFO), *) 'm = ', minm, ':', maxm
    end subroutine read_input

    subroutine read_long_input(filename, f, rv_s, xv_s, ab_s, alpha, lambda, ri_s)
        character(*), intent(in) :: filename
        integer, intent(inout) :: f
        real(knd), intent(inout) :: alpha, lambda
        real(knd), allocatable, dimension(:), intent(inout) :: rv_s, xv_s, ab_s
        complex(knd), allocatable, dimension(:), intent(inout) :: ri_s

        character(1024) :: line, name, eq
        integer :: ios, i, num
        logical :: found_rv, found_xv

        ios = 0
        f = 1
        rv_s = (/1q0/)
        xv_s = (/1q0/)
        ab_s = (/2q0/)
        alpha = 0q0
        lambda = 2q0 * PI
        ri_s = (/cmplx(1.5q0, 0q0, knd)/)
        found_rv = .false.
        found_xv = .false.

        open(FD_GENERAL, file=filename)
        do while (ios == 0)
            read(FD_GENERAL, '(A)', iostat=ios) line
            read(line, *) name, eq

            if (eq == '=') then
                if (name == "f") then
                    read(line, *) name, eq, f
                end if
                if (name == "rv") then
                    read(line, *) name, eq, num
                    deallocate(rv_s)
                    allocate(rv_s(num))
                    read(line, *) name, eq, num, rv_s
                    found_rv = .true.
                end if
                if (name == "xv") then
                    read(line, *) name, eq, num
                    deallocate(xv_s)
                    allocate(xv_s(num))
                    read(line, *) name, eq, num, xv_s
                    found_xv = .true.
                end if
                if (name == "ab") then
                    read(line, *) name, eq, num
                    deallocate(ab_s)
                    allocate(ab_s(num))
                    read(line, *) name, eq, num, ab_s

                end if
                if (name == "alpha") then
                    read(line, *) name, eq, alpha
                end if
                if (name == "lambda") then
                    read(line, *) name, eq, lambda
                end if
                if (name == "ri") then
                    read(line, *) name, eq, num
                    deallocate(ri_s)
                    allocate(ri_s(num))
                    read(line, *) name, eq, num, ri_s
                end if
            end if
        end do
        close(FD_GENERAL)

        alpha = alpha / 180_knd * PI
        alpha = min(alpha, PI / 2q0 - 1q-24)
        if (found_rv) then
            xv_s = 2.0_knd * rv_s * PI / lambda
        else
            rv_s = xv_s * lambda / (2.0_knd * PI)
        end if

        write(FILE_DESCRIPTOR(INFO), *) 'Read input:'
        write(FILE_DESCRIPTOR(INFO), *) 'f = ', f
        write(FILE_DESCRIPTOR(INFO), *) 'rv = ', rv_s
        write(FILE_DESCRIPTOR(INFO), *) 'xv = ', xv_s
        write(FILE_DESCRIPTOR(INFO), *) 'ab = ', ab_s
        write(FILE_DESCRIPTOR(INFO), *) 'alpha = ', alpha
        write(FILE_DESCRIPTOR(INFO), *) 'lambda = ', lambda
        write(FILE_DESCRIPTOR(INFO), *) 'ri = ', ri_s
    end subroutine read_long_input

end module communication