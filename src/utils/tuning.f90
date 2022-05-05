! Created by drakosha on 04.09.2021.

module tuning
    use constants
    use communication
    use spheroidal_scattering
    implicit none
    integer, parameter :: BIG_STEP = 8
    integer, parameter :: SMALL_STEP = 4
    integer, parameter :: BASE_STEP = 2
    real(knd), parameter :: BIG_MULTIPLIER = 1e6_knd
    real(knd), parameter :: SMALL_MULTIPLIER = 1e3_knd

contains
    real(knd) function relative_difference(x, y)
        real(knd) :: x, y

        relative_difference = abs(abs(x) - abs(y)) / abs(abs(x) + abs(y))
    end function relative_difference

    real(knd) function get_difference(tuning_type, current, previous)
        integer, intent(in) :: tuning_type
        type(ScatteringFactors), intent(in) :: current, previous

        if (tuning_type == consequential) then
            get_difference = max(relative_difference(current%ext, previous%ext), &
                    relative_difference(current%sca, previous%sca))
        else
            get_difference = relative_difference(current%ext, current%sca)
        end if
    end function get_difference

    integer function find_optimal_nonsymmetric_matrix_size(tuning_type, f, xv, ab, alpha, lambda, ri, query, fd, guess) result(ms)
        integer, intent(in) :: f, fd, tuning_type
        real(knd), intent(in) :: xv, ab, alpha, lambda
        complex(knd), intent(in) :: ri
        type(ScatteringQuery), intent(inout) :: query
        integer, optional, intent(in) :: guess

        type(SpheroidalScatterer) :: scatterer
        type(WavelengthPoint) :: point
        type(ScatteringResult) :: current_result, previous_result
        integer :: matrix_size
        real(knd) :: rel_diff, current_difference

        call scatterer%set(f, (/xv/), (/ab/), alpha, 1)
        call point%initialize(lambda, 1, (/ cmplx(1.0q0, 0.0q0, knd), ri /))
        
        call query%query_te%do_not_return_arrays()
        call query%query_tm%do_not_return_arrays()
        query%calculate_scattering_indicatrix = .false.

        matrix_size = 8
        if (present(guess)) then
            matrix_size = guess
        end if
        if (tuning_type == consequential) then
            write(fd,*) 'consequential tuning'
        else
            write(fd,*) 'no absorbtion tuning'
        end if
        write(fd,*) 'f = ', f, 'ri = ', ri, 'ab = ', ab, 'xv = ', xv, 'starting from', matrix_size
        do while (matrix_size < MAXIMUM_MATRIX_SIZE)
            query%matrix_size = matrix_size - 4
            previous_result = calculate_scattering(scatterer, point, query)
            query%matrix_size = matrix_size
            current_result = calculate_scattering(scatterer, point, query)
            current_difference = 0
            if (query%query_te%calculate_uv_factors) then
                rel_diff = get_difference(tuning_type, current_result%result_te%uv_factors, previous_result%result_te%uv_factors)
                if (rel_diff < query%accuracy) then
                    query%query_te%calculate_uv_factors = .false.
                    write(fd,*) 'achieved accuracy ', query%accuracy, 'for TE UV for matrix_size = ', matrix_size
                    call current_result%result_te%uv_factors%log(INFO, FILE_DESCRIPTOR(INFO))
                end if
                current_difference = max(current_difference, rel_diff)
            end if
            if (query%query_te%calculate_pq_factors) then
                rel_diff = get_difference(tuning_type, current_result%result_te%pq_factors, previous_result%result_te%pq_factors)
                if (rel_diff < query%accuracy) then
                    query%query_te%calculate_uv_factors = .false.
                    write(fd,*) 'achieved accuracy ', query%accuracy, 'for TE PQ for matrix_size = ', matrix_size
                    call current_result%result_te%pq_factors%log(INFO, FILE_DESCRIPTOR(INFO))
                end if
                current_difference = max(current_difference, rel_diff)
            end if
            if (query%query_tm%calculate_uv_factors) then
                rel_diff = get_difference(tuning_type, current_result%result_tm%uv_factors, previous_result%result_tm%uv_factors)
                if (rel_diff < query%accuracy) then
                    query%query_tm%calculate_uv_factors = .false.
                    write(fd,*) 'achieved accuracy ', query%accuracy, 'for TE UV for matrix_size = ', matrix_size
                    call current_result%result_tm%uv_factors%log(INFO, FILE_DESCRIPTOR(INFO))
                end if
                current_difference = max(current_difference, rel_diff)
            end if
            if (query%query_tm%calculate_pq_factors) then
                rel_diff = get_difference(tuning_type, current_result%result_tm%pq_factors, previous_result%result_tm%pq_factors)
                if (rel_diff < query%accuracy) then
                    query%query_tm%calculate_uv_factors = .false.
                    write(fd,*) 'achieved accuracy ', query%accuracy, 'for TE PQ for matrix_size = ', matrix_size
                    call current_result%result_tm%pq_factors%log(INFO, FILE_DESCRIPTOR(INFO))
                end if
                current_difference = max(current_difference, rel_diff)
            end if
            write(fd,*) 'ms = ', matrix_size, 'rel_diff = ', current_difference
            if (current_difference > BIG_MULTIPLIER * query%accuracy) then
                matrix_size = matrix_size + BIG_STEP
            elseif (current_difference > SMALL_MULTIPLIER * query%accuracy) then
                matrix_size = matrix_size + SMALL_STEP
            elseif (current_difference > query%accuracy) then
                matrix_size = matrix_size + BASE_STEP
            else
                exit
            end if
            previous_result = current_result
        end do
        ms = matrix_size
        write(fd,*)
    end function find_optimal_nonsymmetric_matrix_size
end module tuning