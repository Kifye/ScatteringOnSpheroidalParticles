! Created by drakosha on 13.07.2021.

module spheroidal_tmatrix
    use communication
    use wavelength_point
    use spheroidal_scatterer
    use spheroidal_integrals
    use spheroidal
    use logging
    use matrix
    use spherical
    implicit none

contains
    ! formulas (23)-(26) AppliedOptics 2012
    ! A_ik = W(left_R Delta - mu_j/mu_{j+1}Delta right_R - (mu_j/mu_{j+1}-1)ksi/(ksi^20f)Delta)P
    function get_part(f, ksi, mu, left_R, right_R, W, P, &
            Delta, matrix_size) result(res)
        complex(knd) :: mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i
        complex(knd) :: left_R(matrix_size), right_R(matrix_size), W(matrix_size), P(matrix_size, matrix_size), &
                Delta(matrix_size, matrix_size), tmp(matrix_size, matrix_size), res(matrix_size, matrix_size), val

        tmp = Delta
        call multiply_by_diag_left(tmp, matrix_size, left_R)
        res = tmp

        tmp = Delta
        call multiply_by_diag_right(tmp, matrix_size, right_R)
        res = res - mu(0) / mu(1) * tmp - (mu(0) / mu(1) - 1.0_knd)*ksi / (ksi**2 - f) * Delta
        call multiply_by_diag_left(res, matrix_size, W)

        res = matmul(res, P)
    end function get_part

    subroutine calculate_axisymmetric_tmatrix(scatterer, calculation_point, layers, &
            matrix_size, accuracy, is_te_mode, tmatrix, Delta)

        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalCalculation), intent(in) :: layers(0:1, scatterer%number_of_layers)
        integer :: matrix_size, j, nol
        real(knd), intent(in) :: accuracy
        logical, intent(in) :: is_te_mode
        complex(knd), intent(out) :: tmatrix(matrix_size, matrix_size)
        complex(knd), intent(in) :: Delta(matrix_size, matrix_size, scatterer%number_of_layers)

        integer :: n, i, l
        complex(knd), allocatable, dimension(:, :) :: adder, A11, A31, A31inv
        complex(knd) :: R11(matrix_size), R31(matrix_size), R12(matrix_size), W1(matrix_size), &
                big_matr(2*matrix_size, matrix_size),&
                tmp(2*matrix_size, 2*matrix_size), P(matrix_size, matrix_size), R32(matrix_size)
        complex(knd) :: mu(0:scatterer%number_of_layers)

!        if (.not. outside_layer%functions_calculated) then
!            call log_message('Outside layer relating functions arent calculated!', ERROR, FILE_DESCRIPTOR(ERROR))
!            return
!        end if
!
!        if (.not. inside_layer%functions_calculated) then
!            call log_message('Inside layer relating functions arent calculated!', ERROR, FILE_DESCRIPTOR(ERROR))
!            return
!        end if

        if (scatterer%number_of_layers + 1 < size(calculation_point%mu)) then
            call log_message('Not enough mu data for layers', ERROR, FILE_DESCRIPTOR(ERROR))
            return
        end if

        allocate(adder(matrix_size, matrix_size), &
                A11(matrix_size, matrix_size), &
                A31(matrix_size, matrix_size), A31inv(matrix_size, matrix_size))

        mu = calculation_point%mu
        if (.not. is_te_mode) then
            mu = calculation_point%eps
        endif

        nol = scatterer%number_of_layers
        ! setting up core
        R11 = layers(0, nol)%r1d(1:matrix_size) / layers(0, nol)%r1(1:matrix_size)
        R31 = layers(0, nol)%r3d(1:matrix_size) / layers(0, nol)%r3(1:matrix_size)
        R12 = layers(1, nol)%r1d(1:matrix_size) / layers(1, nol)%r1(1:matrix_size)
        W1 = -1.0_knd / (R31 - R11)
        call get_identity_matrix(P, matrix_size)

        big_matr = 0
        big_matr(matrix_size+1:,:) = get_part(scatterer%spheroidal_type, scatterer%ksi(nol), &
                mu(nol - 1:nol), R11, R12, W1, P, &
                Delta(:,:,nol), matrix_size)
        big_matr(1 : matrix_size,:) = -get_part(scatterer%spheroidal_type, scatterer%ksi(nol), &
                mu(nol - 1:nol), R31, R12, W1, P, &
                Delta(:,:,nol), matrix_size)
        A31 = big_matr(1:matrix_size, :)
        A11 = big_matr(matrix_size + 1:, :)
        do j = nol - 1, 1, -1
            R11 = layers(0, j)%r1d(1:matrix_size) / layers(0, j)%r1(1:matrix_size)
            R31 = layers(0, j)%r3d(1:matrix_size) / layers(0, j)%r3(1:matrix_size)
            R12 = layers(1, j)%r1d(1:matrix_size) / layers(1, j)%r1(1:matrix_size)
            R32 = layers(1, j)%r3d(1:matrix_size) / layers(1, j)%r3(1:matrix_size)
            W1 = -1.0_knd / (R31 - R11)
            tmp = 0
            call set_pi_matrix(layers(1,j), layers(0,j + 1), matrix_size, P)
            do n = 1, matrix_size
                do l = 1, matrix_size
                    if (abs(P(n,l)) < 1q-60) then
                        P(n,l) = 0q0
                    end if
                    P(n,l) = P(n,l) * layers(1,j)%r1(n) / layers(0,j + 1)%r1(l)
                end do
            end do

            tmp(matrix_size + 1: 2*matrix_size, 1:matrix_size) = get_part(scatterer%spheroidal_type, scatterer%ksi(j), &
                    mu(j - 1:j), R11, R12, W1, P, &
                    Delta(:,:,j), matrix_size)
            tmp(1:matrix_size, 1:matrix_size) = -get_part(scatterer%spheroidal_type, scatterer%ksi(j), &
                    mu(j - 1:j), R31, R12, W1, P, &
                    Delta(:,:,j), matrix_size)

            call set_pi_matrix(layers(1,j), layers(0,j + 1), matrix_size, P)
            do n = 1, matrix_size
                do l = 1, matrix_size
                    P(n,l) = P(n,l) * layers(1,j)%r3(n) / layers(0,j + 1)%r3(l)
                end do
            end do

            tmp(matrix_size + 1: 2*matrix_size, matrix_size+1:2*matrix_size) = get_part(scatterer%spheroidal_type, &
                    scatterer%ksi(j), &
                    mu(j - 1:j), R11, R32, W1, P, &
                    Delta(:,:,j), matrix_size)
            tmp(1:matrix_size, matrix_size+1:) = -get_part(scatterer%spheroidal_type, scatterer%ksi(j), &
                    mu(j - 1:j), R31, R32, W1, P, &
                    Delta(:,:,j), matrix_size)
            big_matr = matmul(tmp, big_matr)
        enddo

        A11 = 0
        A31 = 0
        A31inv = 0
        A31 = big_matr(1:matrix_size, :)
        A11 = big_matr(matrix_size + 1:2*matrix_size, :)

        call multiply_by_diag_left(A11, matrix_size, 1.0_knd / layers(0,1)%r3(1:matrix_size))
        call quick_inverse_matrix(A31, matrix_size, A31inv)
        call multiply_by_diag_right(A31inv, matrix_size, layers(0,1)%r1(1:matrix_size))


        tmatrix = -matmul(A11, A31inv)

        deallocate(adder, A11, A31, A31inv)

    end subroutine calculate_axisymmetric_tmatrix

    !  Nonaxisymmetric part

    function get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Epsilon, matrix_size) result(res)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size), identity(matrix_size, matrix_size), &
                res(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        res = -(eps(1) / eps(0) - mu(0) / mu(1)) * f * ksi / (ksi ** 2 - f) * matmul(Q01, Epsilon) - &
                (eps(1) / eps(0) - 1q0) * ksi * (Q01 - 2q0 * ksi**2 * Q01Q11)

        adder = (mu(0) / mu(1) - 1q0) * ksi**2 * Q01 - mu(0) / mu(1) * Delta
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        res = res + adder

        adder = Delta + (eps(1) / eps(0) - 1q0) * ksi**2 * Q01
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        res = res + adder

    end function get_part_11

    function get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1), identity(matrix_size, matrix_size)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                result(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        result = -(eps(1) / eps(0) - mu(0) / mu(1)) * f / (ksi ** 2 - f) * &
                (matmul((ksi**2 * Q01 - Delta), Kappa) + matmul(Delta, Gamma11)) + &
                (eps(1) / eps(0) - 1q0) * f * ksi**2 * 2q0 * matmul(Q01Q11, Gamma11)
        adder = (mu(0) / mu(1) - 1q0) * f * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = (eps(1) / eps(0) - 1q0) * f * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_12

    function get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                result(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        result = (eps(1) / eps(0) - mu(0) / mu(1)) * ksi**2 / (ksi ** 2 - f) * matmul(Q01, Kappa) - &
                (eps(1) / eps(0) - 1q0) * ksi**2 * 2q0 * matmul(Q01Q11, Gamma11)
        adder = -(mu(0) / mu(1) - 1q0) * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = -(eps(1) / eps(0) - 1q0) * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_21

    function get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Epsilon, matrix_size) result(result)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size), identity(matrix_size, matrix_size), &
                result(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        result = (eps(1) / eps(0) - mu(0) / mu(1)) * ksi / (ksi ** 2 - f) * (f * matmul(Q01, Epsilon) + Delta) + &
                (eps(1) / eps(0) - 1q0) * ksi * (Q01 - 2q0 * ksi**2 * Q01Q11)

        adder = -(mu(0) / mu(1) - 1q0) * ksi**2 * Q01 - Delta
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = eps(1) / eps(0) * Delta - (eps(1) / eps(0) - 1q0) * ksi**2 * Q01
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_22

    subroutine set_full_matrix(f, ksi, eps, mu, first_multiplier, left_multiplier, right_multiplier, last_multiplier, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon, matrix_size, result, ddd)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i, j
        real(knd), optional :: ddd
        complex(knd) :: first_multiplier(matrix_size), left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                Epsilon(matrix_size, matrix_size), result(2 * matrix_size, 2 * matrix_size), &
                Q01Q11(matrix_size, matrix_size), last_multiplier(matrix_size, matrix_size)

        result = 0
        result(1:matrix_size, 1:matrix_size) = &
                get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Epsilon, matrix_size)
!        call log_matrix(FILE_DESCRIPTOR(WARNING), 'result', result(1:matrix_size, 1:matrix_size), .false., matrix_size)

        result(1:matrix_size, (matrix_size + 1):(2 * matrix_size)) = &
                get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Kappa, Gamma11, matrix_size)
        result((matrix_size + 1):(2 * matrix_size), 1:matrix_size) = &
                get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Kappa, Gamma11, matrix_size)
        result((matrix_size + 1):(2 * matrix_size), (matrix_size + 1):(2 * matrix_size)) = &
                get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Epsilon, matrix_size)

        do i = 0, 1
            do j = 0, 1
                call multiply_by_diag_left(&
                        result((i * matrix_size + 1):((i + 1) * matrix_size), (j * matrix_size + 1):((j + 1) * matrix_size)), &
                        matrix_size, first_multiplier)
                result((i * matrix_size + 1):((i + 1) * matrix_size), (j * matrix_size + 1):((j + 1) * matrix_size)) =&
                        matmul(result((i * matrix_size + 1):((i + 1) * matrix_size), &
                                (j * matrix_size + 1):((j + 1) * matrix_size)), last_multiplier)

            end do
        end do
        if (present(ddd)) then
            result(:, matrix_size + 1 : 2 * matrix_size) = ddd * result(:, matrix_size + 1 : 2 * matrix_size)
        end if

    end subroutine set_full_matrix

    subroutine calculate_nonaxisymmetric_tmatrix(scatterer, calculation_point, layers, &
            matrix_size, accuracy, is_te_mode, tmatrix, Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon)

        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        ! layers(0, j) = otside_layer for layer j (with ref_ind related to the outside and c, ksi to this layer's size
        ! layers(1,j) - inside_layer
        type(SpheroidalCalculation), intent(in) :: layers(0:1, scatterer%number_of_layers)
        integer, intent(in) :: matrix_size
        real(knd), intent(in) :: accuracy
        logical, intent(in) :: is_te_mode
        complex(knd), intent(out) :: tmatrix(2 * matrix_size, 2 * matrix_size)

        integer :: n, i, l
        complex(knd), allocatable, dimension(:, :) :: A11, A31, A31inv
        complex(knd) :: R11(matrix_size), R31(matrix_size), R12(matrix_size), R32(matrix_size), W1(matrix_size), &
                P(matrix_size, matrix_size)
        complex(knd), intent(in) :: Delta(matrix_size, matrix_size, scatterer%number_of_layers), &
                Q01(matrix_size, matrix_size, scatterer%number_of_layers), &
                Q11(matrix_size, matrix_size, scatterer%number_of_layers), &
                Q01Q11(matrix_size, matrix_size, scatterer%number_of_layers), &
                Kappa(matrix_size, matrix_size, scatterer%number_of_layers), &
                Gamma11(matrix_size, matrix_size, scatterer%number_of_layers), &
                Epsilon(matrix_size, matrix_size, scatterer%number_of_layers)
        complex(knd) :: initial_corrector(2 * matrix_size), solution_corrector(2 * matrix_size)
        complex(knd) :: mu(0:scatterer%number_of_layers), eps(0:scatterer%number_of_layers)
        complex(knd) :: c1
        real(knd) :: k1, ddd
        complex(knd) :: big_matr(4 * matrix_size, 2 * matrix_size), tmp(4 * matrix_size, 4 * matrix_size)
        complex(knd) :: big_matr2(4 * matrix_size, 2 * matrix_size)
        integer :: nol, j

!        if (.not. outside_layer%functions_calculated) then
!            call log_message('Outside layer relating functions arent calculated!', ERROR, FILE_DESCRIPTOR(ERROR))
!            return
!        end if
!
!        if (.not. inside_layer%functions_calculated) then
!            call log_message('Inside layer relating functions arent calculated!', ERROR, FILE_DESCRIPTOR(ERROR))
!            return
!        end if

        if (scatterer%number_of_layers + 1 < size(calculation_point%mu)) then
            call log_message('Not enough mu data for layers', ERROR, FILE_DESCRIPTOR(ERROR))
            return
        end if

        allocate(A11(2 * matrix_size, 2 * matrix_size), A31(2 * matrix_size, 2 * matrix_size), &
                A31inv(2 * matrix_size, 2 * matrix_size))

        mu = calculation_point%mu
        eps = calculation_point%eps
        if (.not. is_te_mode) then
            mu = calculation_point%eps
            eps = calculation_point%mu
        endif

        k1 = calculation_point%k
        c1 = scatterer%c0(1)

        nol = scatterer%number_of_layers
        ! setting up core
        R11 = layers(0, nol)%r1d(1:matrix_size) / layers(0, nol)%r1(1:matrix_size)
        R31 = layers(0, nol)%r3d(1:matrix_size) / layers(0, nol)%r3(1:matrix_size)
        R12 = layers(1, nol)%r1d(1:matrix_size) / layers(1, nol)%r1(1:matrix_size)
        W1 = -1.0_knd / (R31 - R11)

        big_matr = 0
!        P = 1
        call get_identity_matrix(P, matrix_size)
!        write(*,*) 'nol = ', nol
        ! A31
        call set_full_matrix(scatterer%spheroidal_type, scatterer%ksi(nol), eps(nol-1: nol), mu(nol-1: nol), &
                W1, R31, R12, P, &
                Delta(:,:,nol), Q01(:,:,nol), Q11(:,:,nol), Q01Q11(:,:,nol), &
                Kappa(:,:,nol), Gamma11(:,:,nol), Epsilon(:,:,nol), &
                matrix_size, big_matr(1 : 2*matrix_size,:))
        ! A11
        call set_full_matrix(scatterer%spheroidal_type, scatterer%ksi(nol), eps(nol-1: nol), mu(nol-1: nol), &
                W1, R11, R12, P, &
                Delta(:,:,nol), Q01(:,:,nol), Q11(:,:,nol), Q01Q11(:,:,nol), &
                Kappa(:,:,nol), Gamma11(:,:,nol), Epsilon(:,:,nol), &
                matrix_size, big_matr(2*matrix_size + 1 : 4*matrix_size, :))
        big_matr(1 : 2*matrix_size,:) = -big_matr(1 : 2*matrix_size,:)

        do j = nol - 1, 1, -1
            R11 = layers(0, j)%r1d(1:matrix_size) / layers(0, j)%r1(1:matrix_size)
            R31 = layers(0, j)%r3d(1:matrix_size) / layers(0, j)%r3(1:matrix_size)
            R12 = layers(1, j)%r1d(1:matrix_size) / layers(1, j)%r1(1:matrix_size)
            R32 = layers(1, j)%r3d(1:matrix_size) / layers(1, j)%r3(1:matrix_size)
            W1 = -1.0_knd / (R31 - R11)
            ddd = scatterer%d(j) / scatterer%d(j + 1)
            call set_pi_matrix(layers(1,j), layers(0,j + 1), matrix_size, P)
            do n = 1, matrix_size
                do l = 1, matrix_size
                    P(n,l) = P(n,l) * layers(1,j)%r1(n) / layers(0,j + 1)%r1(l)
                end do
            end do

            tmp = 0

            ! A31
            call set_full_matrix(scatterer%spheroidal_type, scatterer%ksi(j), eps(j-1:j), mu(j-1:j), &
                    W1, R31, R12, P, &
                    Delta(:,:,j), Q01(:,:,j), Q11(:,:,j), Q01Q11(:,:,j), Kappa(:,:,j), Gamma11(:,:,j), Epsilon(:,:,j), &
                    matrix_size, tmp(1:2*matrix_size, 1:2*matrix_size), ddd)
            ! A11
            call set_full_matrix(scatterer%spheroidal_type, scatterer%ksi(j), eps(j-1:j), mu(j-1:j), &
                    W1, R11, R12, P, &
                    Delta(:,:,j), Q01(:,:,j), Q11(:,:,j), Q01Q11(:,:,j), Kappa(:,:,j), Gamma11(:,:,j), Epsilon(:,:,j), &
                    matrix_size, tmp(2*matrix_size + 1: 4*matrix_size, 1:2*matrix_size), ddd)

            call set_pi_matrix(layers(1,j), layers(0,j + 1), matrix_size, P)
            do n = 1, matrix_size
                do l = 1, matrix_size
                    P(n,l) = P(n,l) * layers(1,j)%r3(n) / layers(0,j + 1)%r3(l)
                end do
            end do

            ! A33
            call set_full_matrix(scatterer%spheroidal_type, scatterer%ksi(j), eps(j-1:j), mu(j-1:j), &
                    W1, R31, R32, P, &
                    Delta(:,:,j), Q01(:,:,j), Q11(:,:,j), Q01Q11(:,:,j), Kappa(:,:,j), Gamma11(:,:,j), Epsilon(:,:,j), &
                    matrix_size, tmp(1:2*matrix_size, 2*matrix_size + 1: 4*matrix_size), ddd)
            ! A13
            call set_full_matrix(scatterer%spheroidal_type, scatterer%ksi(j), eps(j-1:j), mu(j-1:j), &
                    W1, R11, R32, P, &
                    Delta(:,:,j), Q01(:,:,j), Q11(:,:,j), Q01Q11(:,:,j), Kappa(:,:,j), Gamma11(:,:,j), Epsilon(:,:,j), &
                    matrix_size, tmp(2*matrix_size + 1:4*matrix_size, 2*matrix_size + 1: 4*matrix_size), ddd)
            tmp(1 : 2*matrix_size,:) = -tmp(1 : 2*matrix_size,:)

            big_matr = matmul(tmp, big_matr)

        enddo

        A11 = 0
        A31 = 0
        A31inv = 0
        A31 = big_matr(1:2*matrix_size, :)
        A11 = big_matr(2*matrix_size + 1:4*matrix_size, :)

        initial_corrector(1:matrix_size) = 1q0 / (layers(0, 1)%r1(1:matrix_size) * k1)
        initial_corrector((matrix_size + 1):(2 * matrix_size)) = 1q0 / (layers(0, 1)%r1(1:matrix_size) * c1)
        solution_corrector(1:matrix_size) = 1q0 / (layers(0, 1)%r3(1:matrix_size) * k1)
        solution_corrector((matrix_size + 1):(2 * matrix_size)) = 1q0 / (layers(0, 1)%r3(1:matrix_size) * c1)

        call multiply_by_diag_left(A31, 2 * matrix_size, initial_corrector)
        call multiply_by_diag_left(A11, 2 * matrix_size, solution_corrector)

        call quick_inverse_matrix(A31, 2 * matrix_size, A31inv)
        tmatrix = matmul(A11, A31inv)

!        call log_matrix(FILE_DESCRIPTOR(WARNING), 'nonaxisymm_spheroidal_tmatr', tmatrix, .false., 2 * matrix_size)

        deallocate(A11, A31, A31inv)

    end subroutine calculate_nonaxisymmetric_tmatrix

end module spheroidal_tmatrix