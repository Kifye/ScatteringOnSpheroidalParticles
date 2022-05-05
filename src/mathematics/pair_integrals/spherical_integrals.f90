! Created by drakosha on 18.02.2021.

module spherical_integrals
    use regime
    use elemfunc
    use legendre_functions
    implicit none
contains

    real(knd) function kroneker(n, m)

        integer :: n, m

        if (n == m) then
            kroneker = 1q0
        else
            kroneker = 0q0
        end if
    end function kroneker

    subroutine calculate_spherical_kappa(legendre, Kappa, matrix_size)
        type(LegendreCalculation) :: legendre
        complex(knd), intent(out) :: Kappa(matrix_size, matrix_size)
        integer m, n, l, matrix_size

        m = legendre%m
        Kappa = 0
        !write(*,*) 'spherical kappa'
        do n = m, m + matrix_size - 1
            do l = m, m + matrix_size - 1
                Kappa(l - m + 1, n - m + 1) = -(kroneker(n + 1, l) * (l - 1) * (l - m) / (2 * l - 1) - &
                        kroneker(n - 1, l) * (l + 2) * (l + m + 1) / (2 * l + 3)) / &
                        sqrt(legendre%coef(l - m + 1) / legendre%coef(n - m + 1)) ! qsqrt((2q0 * l + 1) / (2q0 * n + 1))
                !write(*,*) 'n = ', l - m + 1, 'l = ', n - m + 1, 'val = ', Kappa(l - m + 1, n - m + 1)
            enddo
        enddo
        !write(*,*)

    end subroutine calculate_spherical_kappa

    subroutine calculate_spherical_omega(legendre, Omega, matrix_size)
        type(LegendreCalculation) :: legendre
        complex(knd), intent(out) :: Omega(matrix_size, matrix_size)
        integer m, n, l, matrix_size
        m = legendre%m
        Omega = 0

        do n = m, m + matrix_size - 1
            do l = m, m + matrix_size - 1
                Omega(l - m + 1, n - m + 1) = (kroneker(n, l) * 2q0 * (l * l + l + m * m - 1) / (2 * l - 1) / (2 * l + 3) - &
                        kroneker(n + 2, l) * (l - m) * (l - m - 1) / (2 * l - 1) / (2 * l - 3) - &
                        kroneker(n - 2, l) * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 5)) / &
                        sqrt(legendre%coef(l - m + 1) / legendre%coef(n - m + 1))!qsqrt((2q0 * l + 1) / (2q0 * n + 1))
            enddo
        end do
    end subroutine calculate_spherical_omega

    subroutine calculate_spherical_tau(legendre, Tau, matrix_size)
        type(LegendreCalculation) :: legendre
        complex(knd), intent(out) :: Tau(matrix_size, matrix_size)
        integer m, n, l, matrix_size

        m = legendre%m
        Tau = 0
        do n = m, m + matrix_size - 1
            do l = m, m + matrix_size - 1
                !write(*,*) 'n = ', n, 'l = ', l, 'val = ', kroneker(n, l) * l * (l + 1)
                Tau(l - m + 1, n - m + 1) = kroneker(n, l) * l * (l + 1)
            enddo
        end do
    end subroutine calculate_spherical_tau


end module spherical_integrals