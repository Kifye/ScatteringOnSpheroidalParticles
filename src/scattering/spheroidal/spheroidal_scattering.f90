! Created by drakosha on 15.07.2021.

module spheroidal_scattering
    use spheroidal_tmatrix
    use spheroidal_initial
    use spheroidal
    use spheroidal_scatterer
    use wavelength_point
    use angle
    use legendre_functions
    use spherical
    implicit none
    integer, private, parameter :: general_log_fd = 111
contains
    !  Calculte solution
    subroutine set_nonaxisymmetric_solution(scatterer, calculation_point, outside_layer, inside_layer, &
            matrix_size, accuracy, is_te_mode, initial, tmatrix, solution, argnum)

        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalCalculation), intent(in) :: outside_layer, inside_layer
        integer, intent(in) :: matrix_size
        real(knd), intent(in) :: accuracy
        logical, intent(in) :: is_te_mode
        complex(knd), intent(out) :: initial(2 * matrix_size), tmatrix(2 * matrix_size, 2 * matrix_size), &
                solution(2 * matrix_size)
        integer, intent(in) :: argnum

        if (is_te_mode) then
            call set_nonaxisymmetric_initial_te(outside_layer, calculation_point, scatterer%alpha, matrix_size, initial, argnum)
        else
            call set_nonaxisymmetric_initial_tm(outside_layer, calculation_point, scatterer%alpha, matrix_size, initial, argnum)
        end if

        solution = matmul(tmatrix, initial)

    end subroutine set_nonaxisymmetric_solution

    !  Calculate factors
    subroutine set_axisymmetric_solution(scatterer, calculation_point, outside_layer, inside_layer, &
            matrix_size, accuracy, is_te_mode, initial, tmatrix, solution, argnum)

        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalCalculation), intent(in) :: outside_layer, inside_layer
        integer, intent(in) :: matrix_size
        real(knd), intent(in) :: accuracy
        logical, intent(in) :: is_te_mode
        complex(knd), intent(out) :: initial(matrix_size), tmatrix(matrix_size, matrix_size), &
                solution(matrix_size)
        integer, intent(in) :: argnum

        if (is_te_mode) then
            call set_axisymmetric_initial_te(outside_layer, matrix_size, initial, argnum)
        else
            call set_axisymmetric_initial_tm(outside_layer, calculation_point, matrix_size, initial, argnum)
        end if

        !call log_array('initial', matrix_size, initial, matrix_size, FILE_DESCRIPTOR(INFO))
        !call log_matrix('tmatrix', matrix_size, matrix_size, tmatrix, .false., matrix_size, matrix_size, FILE_DESCRIPTOR(INFO))
        solution = matmul(tmatrix, initial)
        !call log_array('solution', matrix_size, solution, matrix_size, FILE_DESCRIPTOR(INFO))

    end subroutine set_axisymmetric_solution

    ! Extinction factor
    function get_symmetric_extinction_factor(scatterer, outside_layer, matrix_size, solution, argnum) result(ext)
        type(SpheroidalScatterer), intent(in) :: scatterer
        type(SpheroidalCalculation), intent(in) :: outside_layer
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(matrix_size)

        integer :: i
        real(knd) :: ext

        ext = 0

        do i = 1, matrix_size
            ext = ext + real(solution(i) * IDEG(mod(3 * i, 4)) * outside_layer%s1(i, argnum), knd)
        enddo

        !ext = -ext * 4.0q0 * scatterer%common_factor

    end function get_symmetric_extinction_factor

    function get_non_symmetric_extinction_factor(scatterer, calculation_point, outside_layer, &
            matrix_size, solution, argnum) result(ext)
        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalCalculation), intent(in) :: outside_layer
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)

        real(knd) :: ext, k1
        integer :: i

        ext = 0
        k1 = calculation_point%k

        do i = 1, matrix_size
            ext = ext + real(NEGIDEG(mod(i + outside_layer%m + 2, 4)) * (k1 * solution(i) * outside_layer%s1(i, argnum) - &
                    solution(i + matrix_size) * NEGIDEG(1) * outside_layer%s1d(i, argnum)), knd)
        enddo

        if (outside_layer%m == 0) then
            ext = ext * 0.5_knd
        end if

    end function get_non_symmetric_extinction_factor

    function get_symmetric_scattering_factor(scatterer, matrix_size, solution) result(sca)
        type(SpheroidalScatterer), intent(in) :: scatterer
        integer, intent(in) :: matrix_size
        complex(knd), intent(in) :: solution(matrix_size)

        real(knd) :: sca

        sca = sum(abs(solution)**2) * 2.0_knd !* scatterer%common_factor

    end function get_symmetric_scattering_factor

    function get_non_symmetric_scattering_factor(scatterer, calculation_point, outside_layer, matrix_size, solution) result(sca)
        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalCalculation), intent(in) :: outside_layer
        integer, intent(in) :: matrix_size
        complex(knd), intent(in) :: solution(2 * matrix_size)

        integer :: i, j, m
        complex(knd) :: Omega(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Tau(matrix_size, matrix_size)
        real(knd) :: sca, k1

        k1 = calculation_point%k
        sca = 0

        call calculate_omega(outside_layer, outside_layer, Omega, matrix_size)
        call calculate_kappa(outside_layer, outside_layer, Kappa, matrix_size)
        call calculate_tau(outside_layer, outside_layer, Tau, matrix_size)
        do i = 1, matrix_size
            do j = 1, matrix_size
                sca = sca + real(IDEG(mod(j + i * 3, 4)) * (&
                        k1**2 * solution(i) * conjg(solution(j)) * Omega(i, j) + &
                                IDEG(1) * k1 * (&
                                        solution(i + matrix_size) * conjg(solution(j)) * Kappa(j, i) - &
                                                conjg(solution(j + matrix_size)) * solution(i) * Kappa(i, j)) + &
                                solution(i + matrix_size) * conjg(solution(j + matrix_size)) * &
                                        Tau(i, j)), knd)
            end do
        end do
        if (outside_layer%m == 0) then
            sca = sca * 0.5_knd
        end if

    end function get_non_symmetric_scattering_factor

    subroutine calculate_base_matrices(layer0, layer1, matrix_size, accuracy, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon)
        type(SpheroidalCalculation), intent(in) :: layer0, layer1
        integer, intent(in) :: matrix_size
        real(knd), intent(in) :: accuracy
        complex(knd) :: Delta(matrix_size, matrix_size), &
                Q01(matrix_size, matrix_size), Q11(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size), &
                Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size)

        complex(knd), allocatable, dimension(:, :) :: tmp, result, identity
        real(knd) :: ksi
        integer :: m, full_size

        if (layer0%m /= layer1%m) then
            call log_message('different m in calculate_base_matrices', ERROR, FILE_DESCRIPTOR(ERROR))
            return
        endif

        m = layer1%m
        ksi = layer0%ksi

        call calculate_kappa(layer1, layer1, Kappa, matrix_size)
        call calculate_gamma(layer1, layer1, Gamma11, matrix_size)
        call calculate_epsilon(layer1, layer1, Epsilon, matrix_size)

        full_size = min(get_full_matrix_size(matrix_size), min(layer0%lnum, layer1%lnum))

        allocate(tmp(full_size, full_size), identity(full_size, full_size), result(full_size, full_size))

        call get_identity_matrix(identity, full_size)
        call calculate_omega(layer1, layer1, tmp, full_size)
        tmp = (ksi**2 - layer0%spheroidal_type) * identity + layer0%spheroidal_type * tmp
        call quick_inverse_matrix(tmp, full_size, result)
        Q11 = result(1:matrix_size, 1:matrix_size)

        identity = 0
        call calculate_delta(layer0, layer1, identity, full_size)
        Delta = 0
        Delta = identity(1:matrix_size, 1:matrix_size)
        tmp = matmul(identity, result)
        Q01 = tmp(1:matrix_size, 1:matrix_size)
        result = matmul(tmp, result)
        Q01Q11 = result(1:matrix_size, 1:matrix_size)

        deallocate(tmp, identity, result)

    end subroutine calculate_base_matrices

    subroutine calculate_symmetric_qfactors(scatterer, calculation_point, query, result, layers, Delta)
        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalCalculation), intent(inout) :: layers(0:1,scatterer%number_of_layers)
        type(ScatteringQuery), intent(in) :: query
        type(ScatteringResult), intent(inout) :: result
        complex(knd) :: Delta(query%matrix_size, query%matrix_size, scatterer%number_of_layers)
        integer :: matrix_size, spherical_matrix_size
        complex(knd), allocatable, dimension(:) :: initial, solution
        complex(knd), allocatable, dimension(:, :) :: tmatrix

        type(LegendreCalculation) :: legendre, legendre0
        complex(knd) :: spherical_initial(query%spherical_matrix_size), init(query%spherical_matrix_size)
        complex(knd) :: spherical_solution(query%spherical_matrix_size)
        complex(knd) :: spherical_tmatrix(query%spherical_matrix_size, query%spherical_matrix_size)
        complex(knd) :: tmp(query%spherical_matrix_size, query%spherical_matrix_size)
        complex(knd) :: connecting_matrix(query%matrix_size, query%spherical_matrix_size)

        matrix_size = query%matrix_size
        spherical_matrix_size = query%spherical_matrix_size
        call legendre%set(layers(0,1)%m, spherical_matrix_size + layers(0,1)%m - 1, scatterer%alpha%angle_cos)
        call legendre%calculate()
        call legendre0%set(0, spherical_matrix_size + 1, scatterer%alpha%angle_cos)
        call legendre0%calculate()
        allocate(initial(matrix_size), solution(matrix_size), tmatrix(matrix_size, matrix_size))
        if (query%query_te%calculate_pq_factors) then
            tmatrix = 0
            initial = 0
            solution = 0
            call calculate_axisymmetric_tmatrix(scatterer, calculation_point, layers, &
                    matrix_size, query%accuracy, .true., tmatrix, Delta)
!            call log_matrix(FILE_DESCRIPTOR(WARNING), 'tmatrix', tmatrix, .false., matrix_size)
            call set_axisymmetric_solution(scatterer, calculation_point, layers(0,1), layers(0,1), &
                    matrix_size, query%accuracy, .true., initial, tmatrix, solution, 1)

                result%result_te%pq_factors%ext = get_symmetric_extinction_factor(&
                        scatterer, layers(0,1), matrix_size, solution, 1) &
                        * 4.0_knd * scatterer%common_factor(1)
                write(general_log_fd, *) 'symmetric Qext_te = ', result%result_te%pq_factors%ext
                result%result_te%pq_factors%sca = get_symmetric_scattering_factor(&
                        scatterer, matrix_size, solution) * scatterer%common_factor(1)
                write(general_log_fd, *) 'symmetric Qsca_te = ', result%result_te%pq_factors%sca
            if (query%query_spherical_te%calculate_pq_factors) then
                call set_connecting_matrix(layers(0,1), matrix_size, spherical_matrix_size, connecting_matrix)
                call set_spherical_tmatrix(matrix_size, tmatrix, &
                        spherical_matrix_size, spherical_tmatrix, connecting_matrix, legendre)
                call calculate_axisymmetric_spherical_initial_te(legendre, &
                        spherical_matrix_size, spherical_initial)
!                write(*,*) 'sph_init = ', spherical_initial
                init = spherical_initial
                spherical_solution = matmul(spherical_tmatrix, spherical_initial)
                result%result_spherical_te%pq_factors%ext = result%result_spherical_te%pq_factors%ext + &
                        get_axisymmetric_spherical_extinction(legendre, calculation_point, scatterer, &
                                spherical_matrix_size, spherical_solution)
                result%result_spherical_te%pq_factors%sca = result%result_spherical_te%pq_factors%sca + &
                        get_axisymmetric_spherical_scattering(legendre, calculation_point, &
                                spherical_matrix_size, spherical_solution)
            end if
        end if
        if (query%query_tm%calculate_pq_factors) then
            tmatrix = 0
            initial = 0
            solution = 0
            call calculate_axisymmetric_tmatrix(scatterer, calculation_point, layers, &
                    matrix_size, query%accuracy, .false., tmatrix, Delta)
!            call log_matrix(FILE_DESCRIPTOR(WARNING), 'tmatrix', tmatrix, .false., matrix_size)
            call set_axisymmetric_solution(scatterer, calculation_point, layers(0,1), layers(1,1), &
                    matrix_size, query%accuracy, .false., initial, tmatrix, solution, 1)

                result%result_tm%pq_factors%ext = get_symmetric_extinction_factor(&
                        scatterer, layers(0,1), matrix_size, solution, 1) &
                        * 4.0_knd * scatterer%common_factor(1)
                write(general_log_fd, *) 'symmetric Qext_tm = ', result%result_tm%pq_factors%ext
                result%result_tm%pq_factors%sca = get_symmetric_scattering_factor(&
                        scatterer, matrix_size, solution) * scatterer%common_factor(1)
                write(general_log_fd, *) 'symmetric Qsca_tm = ', result%result_tm%pq_factors%sca
            if (query%query_spherical_tm%calculate_pq_factors) then
                call set_connecting_matrix(layers(0,1), matrix_size, spherical_matrix_size, connecting_matrix)
                call set_spherical_tmatrix(matrix_size, tmatrix, &
                        spherical_matrix_size, spherical_tmatrix, connecting_matrix, legendre)
                call calculate_axisymmetric_spherical_initial_tm(legendre, &
                        spherical_matrix_size, spherical_initial)
                spherical_solution = matmul(spherical_tmatrix, spherical_initial)
                result%result_spherical_tm%pq_factors%ext = result%result_spherical_tm%pq_factors%ext + &
                        get_axisymmetric_spherical_extinction(legendre, calculation_point, scatterer, &
                                spherical_matrix_size, spherical_solution)
                result%result_spherical_tm%pq_factors%sca = result%result_spherical_tm%pq_factors%sca + &
                        get_axisymmetric_spherical_scattering(legendre, calculation_point, &
                                spherical_matrix_size, spherical_solution)
            end if
        end if
    end subroutine calculate_symmetric_qfactors

    logical function calculate_nonsymmetric_qfactors(scatterer, calculation_point, query, result, layers, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon, tmatrix, initial, solution_te, solution_tm) result(reached_accuracy)
        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalCalculation), intent(inout) :: layers(0:1, scatterer%number_of_layers)
        type(ScatteringQuery), intent(in) :: query
        type(ScatteringResult), intent(inout) :: result
        complex(knd), intent(in) :: Delta(query%matrix_size, query%matrix_size, scatterer%number_of_layers), &
                Q01(query%matrix_size, query%matrix_size, scatterer%number_of_layers), &
                Q11(query%matrix_size, query%matrix_size, scatterer%number_of_layers), &
                Q01Q11(query%matrix_size, query%matrix_size, scatterer%number_of_layers), &
                Kappa(query%matrix_size, query%matrix_size, scatterer%number_of_layers), &
                Gamma11(query%matrix_size, query%matrix_size, scatterer%number_of_layers), &
                Epsilon(query%matrix_size, query%matrix_size, scatterer%number_of_layers)
        complex(knd), intent(inout) :: tmatrix(2 * query%matrix_size, 2 * query%matrix_size), &
                initial(2 * query%matrix_size), solution_te(2 * query%matrix_size), solution_tm(2 * query%matrix_size)
        real(knd) :: part
        integer :: matrix_size, m, spherical_matrix_size
        type(LegendreCalculation) :: legendre
        complex(knd) :: spherical_initial(2 * query%spherical_matrix_size)
        complex(knd) :: spherical_solution(2 * query%spherical_matrix_size)
        complex(knd) :: spherical_tmatrix(2 * query%spherical_matrix_size, 2 * query%spherical_matrix_size)
        complex(knd) :: orth_tmatrix(2 * query%spherical_matrix_size, 2 * query%spherical_matrix_size)
        complex(knd) :: connecting_matrix(query%matrix_size, query%spherical_matrix_size)

        matrix_size = query%matrix_size
        m = layers(0,1)%m
        reached_accuracy = .true.
        spherical_matrix_size = query%spherical_matrix_size
        call legendre%set(m, spherical_matrix_size + m - 1, scatterer%alpha%angle_cos)
        call legendre%calculate()
        if (query%query_te%calculate_uv_factors) then
            tmatrix = 0
            initial = 0
            solution_te = 0
            call calculate_nonaxisymmetric_tmatrix(scatterer, calculation_point, layers, &
                    matrix_size, query%accuracy, .true., tmatrix, Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon)
            call set_nonaxisymmetric_solution(scatterer, calculation_point, layers(0, 1), layers(1, 1), &
                    matrix_size, query%accuracy, .true., initial, tmatrix, solution_te, 1)

                part = get_non_symmetric_extinction_factor(scatterer, calculation_point, layers(0, 1), &
                        matrix_size, solution_te, 1)
!                write(general_log_fd,*) 'm = ', m, 'Qext_te part = ', part * 4.0q0 * &
!                        scatterer%common_factor * scatterer%alpha%angle_sin
                result%result_te%uv_factors%ext = result%result_te%uv_factors%ext + part
                reached_accuracy = reached_accuracy .and. &
                        (abs(part / result%result_te%uv_factors%ext) < query%accuracy)
                part = get_non_symmetric_scattering_factor(scatterer, calculation_point, layers(0, 1), &
                        matrix_size, solution_te)
!                write(general_log_fd,*) 'm = ', m, 'Qsca_te part = ', part * scatterer%common_factor
                result%result_te%uv_factors%sca = result%result_te%uv_factors%sca + part
                reached_accuracy = reached_accuracy .and. &
                        (abs(part / result%result_te%uv_factors%sca) < query%accuracy)

            if (query%query_spherical_te%calculate_uv_factors) then
                call set_connecting_matrix(layers(0, 1), matrix_size, spherical_matrix_size, connecting_matrix)
                call set_nonaxisymmetric_spherical_tmatrix(matrix_size, tmatrix, &
                        spherical_matrix_size, spherical_tmatrix, connecting_matrix, legendre)
!                call log_matrix(FILE_DESCRIPTOR(WARNING), 'T_UV_spherical', tmatrix, .false., 2 * matrix_size)
                call calculate_nonaxisymmetric_spherical_initial_te(legendre, calculation_point, scatterer, &
                        spherical_matrix_size, spherical_initial)
                spherical_solution = matmul(spherical_tmatrix, spherical_initial)
!                call log_array(FILE_DESCRIPTOR(WARNING), 'spherical_a', &
!                        spherical_solution(1:spherical_matrix_size), spherical_matrix_size)
!                call log_array(FILE_DESCRIPTOR(WARNING), 'spherical_b', &
!                        spherical_solution(spherical_matrix_size + 1 : (2 * spherical_matrix_size)), spherical_matrix_size)
                result%result_spherical_te%uv_factors%ext = result%result_spherical_te%uv_factors%ext + &
                        get_nonaxisymmetric_spherical_extinction(legendre, calculation_point, scatterer, &
                                spherical_matrix_size, spherical_solution)
                result%result_spherical_te%uv_factors%sca = result%result_spherical_te%uv_factors%sca + &
                        get_nonaxisymmetric_spherical_scattering(legendre, calculation_point, scatterer, &
                                spherical_matrix_size, spherical_solution)

                if (query%query_spherical_te%return_uv_tmatrix) then
                    call set_spherical_orth_tmatrix(spherical_tmatrix, matrix_size, orth_tmatrix, legendre, calculation_point%k)
                    result%result_spherical_te%uv_tmatrix(:,:,m) = orth_tmatrix
                end if

            end if
        end if

        !  TM mode
        if (query%query_tm%calculate_uv_factors) then
            tmatrix = 0
            initial = 0
            solution_tm = 0
            call calculate_nonaxisymmetric_tmatrix(scatterer, calculation_point, layers, &
                    matrix_size, query%accuracy, .false., tmatrix, Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon)
!            call log_matrix(FILE_DESCRIPTOR(WARNING), 'T_UV_spheroidal', tmatrix, .false., 2 * matrix_size)
            call set_nonaxisymmetric_solution(scatterer, calculation_point, layers(0, 1), layers(1, 1), &
                    matrix_size, query%accuracy, .false., initial, tmatrix, solution_tm, 1)

                part = get_non_symmetric_extinction_factor(scatterer, calculation_point, layers(0, 1), &
                        matrix_size, solution_tm, 1)
!                write(general_log_fd,*) 'm = ', m, 'Qext_tm part = ', part * 4.0_knd * &
!                        scatterer%common_factor * scatterer%alpha%angle_sin
                result%result_tm%uv_factors%ext = result%result_tm%uv_factors%ext + part
                reached_accuracy = reached_accuracy .and. &
                        (abs(part / result%result_tm%uv_factors%ext) < query%accuracy)

                part = get_non_symmetric_scattering_factor(scatterer, calculation_point, layers(0, 1), &
                        matrix_size, solution_tm)
!                write(general_log_fd,*) 'm = ', m, 'Qsca_tm part = ', part * scatterer%common_factor
                result%result_tm%uv_factors%sca = result%result_tm%uv_factors%sca + part
                reached_accuracy = reached_accuracy .and. &
                        (abs(part / result%result_tm%uv_factors%sca) < query%accuracy)
            if (query%query_spherical_tm%calculate_uv_factors) then
                call set_connecting_matrix(layers(0, 1), matrix_size, spherical_matrix_size, connecting_matrix)
                call set_nonaxisymmetric_spherical_tmatrix(matrix_size, tmatrix, &
                        spherical_matrix_size, spherical_tmatrix, connecting_matrix, legendre)
!                call log_matrix(FILE_DESCRIPTOR(WARNING), 'T_UV_spherical', spherical_tmatrix, .false., 2 * matrix_size)
                call calculate_nonaxisymmetric_spherical_initial_tm(legendre, calculation_point, scatterer, &
                        spherical_matrix_size, spherical_initial)
                spherical_solution = matmul(spherical_tmatrix, spherical_initial)
                result%result_spherical_tm%uv_factors%ext = result%result_spherical_tm%uv_factors%ext + &
                        get_nonaxisymmetric_spherical_extinction(legendre, calculation_point, scatterer, &
                                spherical_matrix_size, spherical_solution)
                result%result_spherical_tm%uv_factors%sca = result%result_spherical_tm%uv_factors%sca + &
                        get_nonaxisymmetric_spherical_scattering(legendre, calculation_point, scatterer, &
                                spherical_matrix_size, spherical_solution)

                if (query%query_spherical_tm%return_uv_tmatrix) then
                    call set_spherical_orth_tmatrix(spherical_tmatrix, matrix_size, orth_tmatrix, legendre, calculation_point%k)
                    result%result_spherical_tm%uv_tmatrix(:,:,m) = orth_tmatrix
                end if
            end if
        end if
    end function calculate_nonsymmetric_qfactors

    !  General function iterating over m
    function calculate_scattering(scatterer, calculation_point, query) result(result)
        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(ScatteringQuery), target, intent(in) :: query

        type(ScatteringResult) :: result

        type(SpheroidalCalculation) :: layers(0:1, scatterer%number_of_layers)
        integer :: m, matrix_size, j, nol
!        complex(knd) :: c
        complex(knd), allocatable, dimension(:) :: initial, solution_te, solution_tm
        complex(knd), allocatable, dimension(:, :, :) :: Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon
        complex(knd), allocatable, dimension(:,:) :: tmatrix
        real(knd) :: part
        logical :: calculated_symmetric, reached_accuracy

        matrix_size = query%matrix_size
!        c = scatterer%c0 * calculation_point%get_refractive_index(1)
        nol = scatterer%number_of_layers
        result = create_result(query)

        calculated_symmetric = .false.

        if (query%query_te%calculate_uv_factors .or. query%query_tm%calculate_uv_factors) then
            allocate(Delta(matrix_size, matrix_size, nol), Q01(matrix_size, matrix_size, nol), &
                    Q11(matrix_size, matrix_size, nol), &
                    Q01Q11(matrix_size, matrix_size, nol), Kappa(matrix_size, matrix_size, nol), &
                    Gamma11(matrix_size, matrix_size, nol), &
                    Epsilon(matrix_size, matrix_size, nol))
            allocate(tmatrix(2 * matrix_size, 2 * matrix_size), initial(2 * matrix_size), &
                    solution_te(2 * matrix_size), solution_tm(2 * matrix_size))

            do m = query%minm, query%maxm
                ! layer(0, j) outside layer of the j layer of the particle
                ! layer(1, j) - inside layer of j layer
                ! j = 1 - mantle
                ! j = 2 - core
                do j = 1, scatterer%number_of_layers
                call layers(0,j)%calculate(m, get_full_function_size(matrix_size), &
                        scatterer%c0(j) * calculation_point%get_refractive_index(j - 1), scatterer%ksi(j), &
                        1, (/ scatterer%alpha%angle_cos /), scatterer%spheroidal_type)
                call layers(1,j)%calculate(m, get_full_function_size(matrix_size), &
                        scatterer%c0(j) * calculation_point%get_refractive_index(j), scatterer%ksi(j), &
                        1, (/ scatterer%alpha%angle_cos /), scatterer%spheroidal_type)
                end do
                Delta = 0
                Q01 = 0
                Q11 = 0
                Q01Q11 = 0
                Kappa = 0
                Gamma11 = 0
                Epsilon = 0
                do j = 1, scatterer%number_of_layers
                call calculate_base_matrices(layers(0, j), layers(1, j), query%matrix_size, query%accuracy, &
                        Delta(:,:,j), Q01(:,:,j), Q11(:,:,j), Q01Q11(:,:,j), Kappa(:,:,j), Gamma11(:,:,j), Epsilon(:,:,j))
                end do

                reached_accuracy = calculate_nonsymmetric_qfactors(scatterer, calculation_point, query, result, &
                        layers, Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon, tmatrix, initial, &
                        solution_te, solution_tm)
                if (query%query_te%return_uv_solution) then
                    result%result_te%uv_solution(:,m) = solution_te
                end if
                if (query%query_tm%return_uv_solution) then
                    result%result_tm%uv_solution(:,m) = solution_tm
                end if
                result%actual_maxm = m
                ! TODO return when layered PQ potentials
                if (m == 1 .and. (query%query_te%calculate_pq_factors .or. query%query_tm%calculate_pq_factors)) then
                    call calculate_symmetric_qfactors(scatterer, calculation_point, query, result, &
                            layers, Delta)
                    calculated_symmetric = .true.
                end if

                if (reached_accuracy) then
                    exit
                end if

            end do

            result%result_te%uv_factors%ext = result%result_te%uv_factors%ext * &
                    4.0_knd * scatterer%common_factor(1) * scatterer%alpha%angle_sin
            result%result_tm%uv_factors%ext = result%result_tm%uv_factors%ext * &
                    4.0_knd * scatterer%common_factor(1) * scatterer%alpha%angle_sin
            result%result_te%uv_factors%sca = result%result_te%uv_factors%sca * scatterer%common_factor(1)
            result%result_tm%uv_factors%sca = result%result_tm%uv_factors%sca * scatterer%common_factor(1)

            deallocate(Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon)
            deallocate(tmatrix, initial, solution_te, solution_tm)
        end if

        if ((query%query_te%calculate_pq_factors .or. query%query_tm%calculate_pq_factors) .and. .not. calculated_symmetric) then
            m = 1
            allocate(Delta(matrix_size, matrix_size, scatterer%number_of_layers))
            do j = 1, scatterer%number_of_layers
                call layers(0,j)%calculate(m, get_full_function_size(matrix_size), &
                        scatterer%c0(j) * calculation_point%get_refractive_index(j - 1), scatterer%ksi(j), &
                        1, (/ scatterer%alpha%angle_cos /), scatterer%spheroidal_type)
                call layers(1,j)%calculate(m, get_full_function_size(matrix_size), &
                        scatterer%c0(j) * calculation_point%get_refractive_index(j), scatterer%ksi(j), &
                        1, (/ scatterer%alpha%angle_cos /), scatterer%spheroidal_type)
                call calculate_delta(layers(0,j), layers(1,j), Delta(:,:,j), matrix_size)
            end do
            call calculate_symmetric_qfactors(scatterer, calculation_point, query, result, &
                    layers, Delta)
            deallocate(Delta)
        end if

        call log_message('eventual result', INFO, FILE_DESCRIPTOR(INFO))
        call result%log(INFO, FILE_DESCRIPTOR(INFO))
    end function calculate_scattering
end module spheroidal_scattering