module scattering_functions
    use communication
    use constants
    use geometry
    use spheroidal_scatterer
    use spheroidal_scattering
    use spheroidal_indicatrix
    use wavelength_point
    implicit none
contains
!    integer(4) function FindOptimalMatrixSizeExtSca(f, xv, ab, alpha, lambda, ri, relative_accuracy)
!        integer, intent(in) :: f
!        real(knd), intent(in) :: xv, ab, alpha, lambda, relative_accuracy
!        complex(knd), intent(in) :: ri
!
!        type(SpheroidalScatterer) :: scatterer
!        type(WavelengthPoint) :: point
!        type(ScatteringQuery) :: query
!        type(ScatteringResult) :: result
!
!        call scatterer%set(f, xv, ab, alpha, 1)
!        call point%initialize(lambda, 1, (/ cmplx(1.0q0, 0.0q0, knd), ri /))
!
!        FindOptimalMatrixSizeExtSca = 0
!        do matrix_size = 4, MAXIMUM_MATRIX_SIZE, 2
!            query%matrix_size = reate_query(size_of_matrices(f, xv, ab, ri), &
!                    0, get_max_m(f, xv, ab), &
!                      .true., .true., .true., .true., &
!                    .true., .true., .true., .true., &
!                    !.false., .false., .false., .false., &
!                    BASE_ACCURACY)
!            result = calculate_scattering(scatterer, point, query)
!
!            if (FindOptimalMatrixSizeExtSca(0) == 0 .and. &
!                    relative_difference(&
!                            result%symmetric_factors_value%Qext_te, &
!                            result%symmetric_factors_value%Qsca_te) < relative_accuracy) then
!                FindOptimalMatrixSizeExtSca(0) = matrix_size
!            end if
!            if (FindOptimalMatrixSizeExtSca(1) == 0 .and. &
!                    relative_difference(&
!                            result%symmetric_factors_value%Qext_tm, &
!                            result%symmetric_factors_value%Qsca_tm) < relative_accuracy) then
!                FindOptimalMatrixSizeExtSca(1) = matrix_size
!            end if
!            if (FindOptimalMatrixSizeExtSca(2) == 0 .and. &
!                    relative_difference(&
!                            result%nonsymmetric_factors_value%Qext_te, &
!                            result%nonsymmetric_factors_value%Qsca_te) < relative_accuracy) then
!                FindOptimalMatrixSizeExtSca(2) = matrix_size
!            end if
!            if (FindOptimalMatrixSizeExtSca(3) == 0 .and. &
!                    relative_difference(&
!                            result%nonsymmetric_factors_value%Qext_tm, &
!                            result%nonsymmetric_factors_value%Qsca_tm) < relative_accuracy) then
!                FindOptimalMatrixSizeExtSca(3) = matrix_size
!            end if
!
!            if (FindOptimalMatrixSizeExtSca(0) > 0 .and. FindOptimalMatrixSizeExtSca(1) > 0 .and. &
!                    FindOptimalMatrixSizeExtSca(2) > 0 .and. FindOptimalMatrixSizeExtSca(3) > 0) then
!                exit
!            end if
!        end do
!    end function FindOptimalMatrixSizeExtSca

    type(ScatteringResult) function GetSpheroidalQFactors(f, number_of_layers, xv, ab, alpha, lambda, ri, matrix_size, minm, maxm)
        integer, intent(in) :: f, number_of_layers
        real(knd), intent(in) :: xv(number_of_layers), ab(number_of_layers), alpha, lambda
        complex(knd), intent(in) :: ri(0:number_of_layers)
        integer, intent(in) :: matrix_size
        integer, optional, intent(in) :: minm, maxm

        type(SpheroidalScatterer) :: scatterer
        type(WavelengthPoint) :: point
        type(ScatteringQuery) :: query
        type(ScatteringResult) :: result
        integer :: max_border, min_border

        min_border = 0
        max_border = get_max_m(f, xv(1), ab(1))
        if (present(minm)) then
            min_border = minm
        end if
        if (present(maxm)) then
            max_border = maxm
        end if

        call scatterer%set(f, xv, ab, alpha, number_of_layers)
        query = create_query(matrix_size, & !size_of_matrices(f, xv(1), ab(1), ri(1)), &
                min_border, max_border, &
                BASE_ACCURACY, &
                uv_te = .true., uv_tm = .true., pq_te = .true., pq_tm = .true., double_for_spherical = .false.)
!                uv_te = .true., uv_tm = .true., pq_te = .true., pq_tm = .true., double_for_spherical = .true.)
!                .true., .true., .true., .true.)
!                double_for_spherical=.true., spherical_matrix_size=2 * size_of_matrices(f, xv, ab, ri))
        call point%initialize(lambda, number_of_layers, ri)
        GetSpheroidalQFactors = calculate_scattering(scatterer, point, query)

    end function GetSpheroidalQFactors

    type(ScatteringResult) function GetSphericalTmatrix(f, number_of_layers, xv, ab, alpha, lambda, ri, matrix_size, minm, maxm)
        integer, intent(in) :: f, number_of_layers
        real(knd), intent(in) :: xv(number_of_layers), ab(number_of_layers), alpha, lambda
        complex(knd), intent(in) :: ri(0:number_of_layers)
        integer, intent(in) :: matrix_size
        integer, optional, intent(in) :: minm, maxm

        type(SpheroidalScatterer) :: scatterer
        type(WavelengthPoint) :: point
        type(ScatteringQuery) :: query
        type(ScatteringResult) :: result
        integer :: max_border, min_border

        min_border = 0
        max_border = get_max_m(f, xv(1), ab(1))
        if (present(minm)) then
            min_border = minm
        end if
        if (present(maxm)) then
            max_border = maxm
        end if

        call scatterer%set(f, xv, ab, alpha, number_of_layers)
        query = create_query(matrix_size, &
                min_border, max_border, &
                BASE_ACCURACY, &
                uv_te = .true., uv_tm = .true., return_tmatrix = .true., double_for_spherical = .true.)
        call point%initialize(lambda, number_of_layers, ri)
        GetSphericalTmatrix = calculate_scattering(scatterer, point, query)

    end function GetSphericalTmatrix

    subroutine print_indicatrix(f, xv, ab, alpha, lambda, ri, minm, maxm, &
            ntheta, theta_begin, theta_end, nphi, phi_begin, phi_end, fd)
        integer, intent(in) :: f, minm, maxm, ntheta, nphi, fd
        real(knd), intent(in) :: xv, ab, alpha, lambda, theta_begin, theta_end, phi_begin, phi_end
        complex(knd), intent(in) :: ri

        type(SpheroidalScatterer) :: scatterer
        type(WavelengthPoint) :: point
        type(ScatteringQuery) :: query
        type(ScatteringResult) :: result
        integer :: matrix_size, i, m
        type(AngleType) :: thetas(ntheta + 1), phis(nphi + 1, minm:maxm)
        real(knd) :: dtheta, dphi, theta_value

        call scatterer%set(f, (/xv/), (/ab/), alpha, 1)
        matrix_size = size_of_matrices(f, xv, ab, ri)
        query = create_query(matrix_size, &
                minm, maxm, BASE_ACCURACY, &
                .true., return_solution = .true.)
        call point%initialize(lambda, 1, (/ cmplx(1.0q0, 0.0q0, knd), ri /))
        result = calculate_scattering(scatterer, point, query)
        if (ntheta == 0) then
            dtheta = 0
        else
            dtheta = (theta_end - theta_begin) / ntheta
        end if
        if (nphi == 0) then
            dphi = 0
        else
            dphi = (phi_end - phi_begin) / nphi
        end if
        do i = 1, ntheta + 1
            theta_value = theta_begin + dtheta * (i - 1)
            if (abs(theta_value) < 1q-16) then
                theta_value = theta_value + 1q-12
            end if
            call thetas(i)%set(theta_begin + dtheta * (i - 1))
        end do
        do m = query%minm, result%actual_maxm
            do i = 1, nphi + 1
                call phis(i, m)%set(m * (phi_begin + dphi * (i - 1)))
            end do
        end do
        call print_scattering_indicatrix(scatterer, point, query%minm, result%actual_maxm, &
                matrix_size, result%result_te%uv_solution, result%result_tm%uv_solution, &
                ntheta, thetas, nphi, phis(:, query%minm:result%actual_maxm), fd)
    end subroutine print_indicatrix

!    subroutine print_indicatrix_in_xz_plane(f, xv, ab, alpha, lambda, ri, minm, maxm, &
!            nxi, xi_begin, xi_end, fd)
!        integer, intent(in) :: f, minm, maxm, nxi, fd
!        real(knd), intent(in) :: xv, ab, alpha, lambda, xi_begin, xi_end
!        complex(knd), intent(in) :: ri
!
!        type(SpheroidalScatterer) :: scatterer
!        type(WavelengthPoint) :: point
!        type(ScatteringQuery) :: query
!        type(ScatteringResult) :: result
!        integer :: matrix_size, i, m
!        type(AngleType) :: thetas(ntheta + 1), phis(nphi + 1, minm:maxm)
!        real(knd) :: dtheta, theta_value
!
!        call scatterer%set(f, xv, ab, alpha, 1)
!        matrix_size = size_of_matrices(f, xv, ab, ri)
!        query = create_query(matrix_size, &
!                minm, maxm, &
!                .true., .true., .true., .true., &
!                .true., .true., .true., .true., &
!                !                .false., .false., .false., .false., &
!                BASE_ACCURACY)
!        call point%initialize(lambda, 1, (/ cmplx(1.0q0, 0.0q0, knd), ri /))
!        result = calculate_scattering(scatterer, point, query)
!        if (nxi == 0) then
!            dxi = 0
!        else
!            dxi = (xi_end - xi_begin) / nxi
!        end if
!        do i = 1, nxi + 1
!            theta_value = xi_begin + dxi * (i - 1)
!            if (abs(theta_value - alpha) < 1q-16) then
!                theta_value = theta_value + 1q-12
!            end if
!            if (theta_value < scatterer%alpha) then
!            call thetas(i)%set(abs(alpha - theta_value))
!
!            !end if
!        end do
!        do m = query%minm, result%actual_maxm
!            do i = 1, nphi + 1
!                call phis(i, m)%set(m * (phi_begin + dphi * (i - 1)))
!            end do
!        end do
!        call print_scattering_indicatrix(scatterer, point, query%minm, result%actual_maxm, &
!                matrix_size, result%solution_te, result%solution_tm, &
!                ntheta, thetas, nphi, phis(:, query%minm:result%actual_maxm), fd)
!    end subroutine print_indicatrix_in_xz_plane


end module scattering_functions