! Created by drakosha on 18.02.2021.

module spherical
    use regime
    use legendre_functions
    use constants
    use spheroidal
    use spherical_integrals
    use spheroidal_scatterer
    use matrix
    use constants
    use wavelength_point
    use logging

    implicit none

contains

    subroutine set_connecting_matrix(layer, nsize, lsize, res)
        type(SpheroidalCalculation), intent(in) :: layer
        integer :: m, lsize, nsize, i, j
        complex(knd), intent(out) :: res(nsize, lsize)

        type(LegendreCalculation) :: legendre

        if (nsize > layer%lnum) then
            write(*,*) 'ERROR: nsize = ', nsize, 'lnum = ', layer%lnum
            nsize = layer%lnum
        end if
        call legendre%set(layer%m, layer%m + layer%lnum - 1, 2.0_knd)
        call legendre%calculate()
        do i = 1, nsize
            do j = 1, lsize
!                res(i, j) = cmplx(0q0, 1q0, knd) ** (j - i) * layer%legendre(j - 1, i)
                res(i, j) = cmplx(0q0, 1q0, knd) ** (j - i) * layer%legendre(j - 1, i) * sqrt(1q0 / legendre%coef(j))
            end do
        end do

    end subroutine set_connecting_matrix

    subroutine set_pi_matrix(first, second, matrix_size, P)
        type(SpheroidalCalculation), intent(in) :: first, second
        type(LegendreCalculation) :: leg
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: P(matrix_size, matrix_size)

        complex(knd) :: tmp1(matrix_size, matrix_size), tmp2(matrix_size, matrix_size), value
        integer :: n, l, s, maxd, lnum
        lnum = min(first%maxd, second%maxd)
        call leg%set(first%m, first%m + lnum + 1, 1q0)
        call leg%calculate()
        P = 0
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(n,2) /= mod(l,2)) then
                    cycle
                end if
                do s = mod(n + 1, 2), lnum, 2
                    value = cmplx(0q0, 1q0) ** (l-n) * first%legendre(s, n) * second%legendre(s, l) * &
                            1q0 / leg%coef(s+1)
                    P(n,l) = P(n,l) + value
                    if (abs(value) < DEFAULT_PAIR_INTEGRAL_ACCURACY * abs(P(n,l))) then
                        exit
                    end if
                end do
            end do
        end do
    end subroutine set_pi_matrix

    subroutine set_spherical_tmatrix(spheroidal_size, spheroidal_tmatrix, spherical_size, spherical_tmatrix, &
            connecting_matrix, legendre)
        integer, intent(in) :: spheroidal_size, spherical_size
        complex(knd), intent(in) :: spheroidal_tmatrix(spheroidal_size, spheroidal_size), &
                connecting_matrix(spheroidal_size, spherical_size)
        complex(knd), intent(out) :: spherical_tmatrix(spherical_size, spherical_size)
        type(LegendreCalculation), intent(in) :: legendre

        complex(knd) :: h(spherical_size)
        integer :: n

        spherical_tmatrix = matmul(transpose(connecting_matrix), matmul(spheroidal_tmatrix, connecting_matrix))
        call multiply_by_diag_right(spherical_tmatrix, spherical_size, cmplx(1q0, 0q0) * sqrt(legendre%coef))
        call multiply_by_diag_left(spherical_tmatrix, spherical_size, cmplx(1q0, 0q0) / sqrt(legendre%coef))

    end subroutine set_spherical_tmatrix

    subroutine set_nonaxisymmetric_spherical_tmatrix(n, spheroidal_tmatrix, m, spherical_tmatrix, connecting_matrix, legendre)
        integer, intent(in) :: n, m
        complex(knd), intent(in) :: spheroidal_tmatrix(2 * n, 2 * n), connecting_matrix(n, m)
        complex(knd), intent(out) :: spherical_tmatrix(2 * m, 2 * m)
        type(LegendreCalculation), intent(in) :: legendre

        complex(knd) :: f(0:m, 0:m), orth(m, m), revf(m, m)

        integer :: i, j
        do i = 0, 1
            do j = 0, 1
                call set_spherical_tmatrix(n, spheroidal_tmatrix((i * n + 1):((i + 1) * n), (j * n + 1):((j + 1) * n)), &
                        m, spherical_tmatrix((i * m + 1):((i + 1) * m), (j * m + 1):((j + 1) * m)), connecting_matrix, &
                legendre)
            end do
        end do

    end subroutine set_nonaxisymmetric_spherical_tmatrix

    subroutine calculate_nonaxisymmetric_spherical_initial_te(legendre, calculation_point, scatterer, matrix_size, initial)
        type(LegendreCalculation), intent(in) :: legendre
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalScatterer), intent(in) :: scatterer
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(2 * matrix_size)

        integer :: i

        initial = 0

        do i = 1, matrix_size
            initial(i) = IDEG(mod(i + legendre%m + 2, 4)) * legendre%pr(1, i) / sqrt(legendre%coef(i))
        enddo

        initial = initial * (4q0 / (calculation_point%k * scatterer%alpha%angle_sin))

    end subroutine calculate_nonaxisymmetric_spherical_initial_te

    subroutine calculate_nonaxisymmetric_spherical_initial_tm(legendre, calculation_point, scatterer, matrix_size, initial)
        type(LegendreCalculation), intent(in) :: legendre
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalScatterer), intent(in) :: scatterer
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(2 * matrix_size)

        call calculate_nonaxisymmetric_spherical_initial_te(legendre, calculation_point, scatterer, matrix_size, initial)
        initial = -initial

    end subroutine calculate_nonaxisymmetric_spherical_initial_tm

    subroutine calculate_axisymmetric_spherical_initial_te(legendre, matrix_size, initial)
        type(LegendreCalculation), intent(in) :: legendre
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(matrix_size)
        integer :: i

        initial = 0

        do i = 1, matrix_size
        initial(i) = -2q0 * IDEG(mod(i + legendre%m + 3, 4)) * legendre%pr(1, i) / sqrt(legendre%coef(i))

        end do
    end subroutine calculate_axisymmetric_spherical_initial_te

    subroutine calculate_axisymmetric_spherical_orthogonal_initial_te(legendre, calculation_point, scatterer, matrix_size, initial)
        type(LegendreCalculation), intent(in) :: legendre
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalScatterer), intent(in) :: scatterer
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(matrix_size)
        integer :: i

        initial = 0

        do i = 1, matrix_size
            initial(i) = -IDEG(mod(i + legendre%m + 3, 4)) * legendre%pdr(1, i) &
            * (2q0 * i + 1q0) / (calculation_point%k * scatterer%alpha%angle_sin)
        end do
    end subroutine calculate_axisymmetric_spherical_orthogonal_initial_te

    subroutine calculate_axisymmetric_spherical_initial_tm(legendre, matrix_size, initial)
        type(LegendreCalculation), intent(in) :: legendre
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(matrix_size)

        call calculate_axisymmetric_spherical_initial_te(legendre, matrix_size, initial)
        initial = -initial

    end subroutine calculate_axisymmetric_spherical_initial_tm

    real(knd) function get_axisymmetric_spherical_extinction(legendre, calculation_point, &
            scatterer, matrix_size, solution) result(ext)
        type(LegendreCalculation), intent(in) :: legendre
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalScatterer), intent(in) :: scatterer
        integer, intent(in) :: matrix_size
        complex(knd), intent(in) :: solution(matrix_size)
        integer :: i
        ext = 0q0

        do i = 1, matrix_size
            ext = ext + real(NEGIDEG(mod(i + legendre%m + 3, 4)) * solution(i) * legendre%pr(1, i) * sqrt(legendre%coef(i)), knd)
        enddo

        ext = ext * 4q0 * PI / calculation_point%k ** 2
    end function get_axisymmetric_spherical_extinction

    real(knd) function get_nonaxisymmetric_spherical_extinction(legendre, calculation_point, scatterer, &
            matrix_size, solution) result(ext)
        type(LegendreCalculation), intent(in) :: legendre
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalScatterer), intent(in) :: scatterer
        integer, intent(in) :: matrix_size
        complex(knd), intent(in) :: solution(2 * matrix_size)

        integer :: i
        ext = 0q0

        do i = 1, matrix_size
            ext = ext + real(NEGIDEG(mod(i + legendre%m + 2, 4)) * (calculation_point%k * solution(i) * legendre%pr(1, i) + &
                    solution(i + matrix_size) * IDEG(1) * legendre%pdr(1, i)), knd) * sqrt(legendre%coef(i))! / (this%legendre%coef(i))
        enddo

        ext = -ext * 4q0 * PI / calculation_point%k ** 2 * scatterer%alpha%angle_sin
        if (legendre%m == 0) then
            ext = ext / 2
        end if
    end function get_nonaxisymmetric_spherical_extinction

    real(knd) function get_axisymmetric_spherical_scattering(legendre, calculation_point, matrix_size, solution) result(sca)
        type(LegendreCalculation), intent(in) :: legendre
        type(WavelengthPoint), intent(in) :: calculation_point
        integer, intent(in) :: matrix_size
        complex(knd), intent(in) :: solution(matrix_size)

        integer :: i

        sca = 0

        do i = 1, matrix_size
            sca = sca + abs(solution(i)) ** 2! / legendre%coef(i)
        end do
        sca = sca * 2q0 * PI / calculation_point%k**2

    end function get_axisymmetric_spherical_scattering

    real(knd) function get_axisymmetric_spherical_orthogonal_scattering(calculation_point, matrix_size, solution) result(sca)
        type(WavelengthPoint), intent(in) :: calculation_point
        integer, intent(in) :: matrix_size
        complex(knd), intent(in) :: solution(matrix_size)

        integer :: i

        sca = 0

        do i = 1, matrix_size
            sca = sca + abs(solution(i)) ** 2
        end do
        sca = sca * 2q0 * PI / calculation_point%k**2

    end function get_axisymmetric_spherical_orthogonal_scattering

    real(knd) function get_nonaxisymmetric_spherical_scattering(legendre, calculation_point, &
            scatterer, matrix_size, solution) result(sca)
        type(LegendreCalculation), intent(in) :: legendre
        type(WavelengthPoint), intent(in) :: calculation_point
        type(SpheroidalScatterer), intent(in) :: scatterer
        integer, intent(in) :: matrix_size
        complex(knd), intent(in) :: solution(2 * matrix_size)

        integer :: i, j, m
        complex(knd) :: ideg1, Omega(matrix_size, matrix_size), &
                Kappa(matrix_size, matrix_size), Tau(matrix_size, matrix_size)
        real(knd) :: k1, coef

        k1 = calculation_point%k

        sca = 0
        ideg1 = cmplx(0q0, 1q0, knd)

        call calculate_spherical_omega(legendre, Omega, matrix_size)
        call calculate_spherical_kappa(legendre, Kappa, matrix_size)
        call calculate_spherical_tau(legendre, Tau, matrix_size)
        do i = 1, matrix_size
            do j = 1, matrix_size
                !coef = qsqrt(legendre%coef(i) * legendre%coef(j))
                sca = sca + real(ideg1**(j - i) * (&
                        k1**2 * solution(i) * conjg(solution(j)) * Omega(i, j) + &
                                IDEG(1) * k1 * (&
                                        solution(i + matrix_size) * conjg(solution(j)) * &
                                                Kappa(j, i) - &
                                                conjg(solution(j + matrix_size)) * solution(i) *&
                                                        Kappa(i, j)) + &
                                solution(i + matrix_size) * conjg(solution(j + matrix_size)) * &
                                        Tau(i, j)), knd)! / coef
            end do
        end do
        sca = sca * PI / k1 ** 2
        if (legendre%m == 0) then
            sca = sca / 2q0
        end if

    end function get_nonaxisymmetric_spherical_scattering

    subroutine set_spherical_orth_tmatrix(nonorth_tmatrix, matrix_size, orth_tmatrix)
        integer :: matrix_size, n, l
        complex(knd) :: nonorth_tmatrix(2 * matrix_size, 2 * matrix_size), f(matrix_size, matrix_size), &
                orth_tmatrix(2 * matrix_size, 2 * matrix_size), rev(matrix_size, matrix_size)
        do n = 1, matrix_size
            do l = 1, matrix_size
                f(n, l) = (kroneker(n+1, l) + kroneker(n-1, l)) / sqrt((2 * n + 1q0)*(2q0*l+1q0))
            end do
        end do

        call inverse_matrix(f, matrix_size, rev)

        orth_tmatrix=0
        orth_tmatrix = matmul(matmul(f, nonorth_tmatrix(1:matrix_size, 1:matrix_size)) + &
                nonorth_tmatrix(matrix_size + 1:matrix_size, 1:matrix_size), rev)

    end subroutine set_spherical_orth_tmatrix

end module spherical