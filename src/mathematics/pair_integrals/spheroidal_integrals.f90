!  contains functions for calculation of pair integrals of spheroidal functions
!  given by 2 variables of the type SpheroidalCalculation from module
!  spheroidal as matrices of a requested size
module spheroidal_integrals
    use regime
    use spheroidal
    use constants
    use logging
    implicit none

    real(knd), private, parameter :: accuracy = DEFAULT_PAIR_INTEGRAL_ACCURACY
    integer, private, parameter :: bucket = 20
contains
    logical function update_with_accuracy(value, delta, accuracy)
        complex(knd), intent(inout) :: value
        complex(knd), intent(in) :: delta
        real(knd), intent(in) :: accuracy

        value = value + delta
        update_with_accuracy = (abs(delta / value) < accuracy)
    end function update_with_accuracy

    !  fills the array common_multiplier with real numbers
    !  (r+2m)/r!*2/(2r+2m+1) for r = 0..accuracy
    !  common_multiplier is allocated if necessary
    subroutine fill_common_multiplier(m, accuracy, common_multiplier)
        integer :: m, accuracy, r
        real(knd), allocatable, dimension(:) :: common_multiplier
        real(knd) :: factorial

        if (allocated(common_multiplier) .and. size(common_multiplier) < accuracy + 1) then
            deallocate(common_multiplier)
        end if
        if (.not. allocated(common_multiplier)) then
            allocate(common_multiplier(0:accuracy))
        end if

        factorial = 1q0
        do r = 2, 2 * m
            factorial = factorial * r
        end do
        do r = 0, accuracy
            common_multiplier(r) = factorial * 2q0 / (2 * r + 2 * m + 1)
            factorial = factorial * (r + 2 * m + 1) / (r + 1)
        enddo
    end subroutine fill_common_multiplier

    !  Delta_{n,l}^{(m)}(c_i, c_j) = \int_{-1}^{1} S_{m,n}(c_i, \eta) S_{m,n}(c_j, \eta) d\eta =
    !  \sum'_{r=0}^{accuracy} d_r^{m,n}(c_i) d_r^{m,l}(c_j) common_multiplier(m, r)
    !  in sum' only r = n-l (mod 2) are used,
    !  angular spheroidal functions and legendre coefficients are normalized
    !  and taken from SpheroidalCalculation variables
    !  common_multiplier is calculated in subroutine fill_common_multiplier
    subroutine calculate_delta(first, second, Delta, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: Delta(matrix_size, matrix_size)

        real(knd), allocatable, dimension(:) :: common_multiplier
        integer n, l, r, k
        real :: start, finish
        complex(knd) :: value

        Delta = 0q0
        if (first%m /= second%m) then
            write(*, *) 'different m in delta'
            return
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'delta too few elements in first', first%lnum
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if
        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)

        !write(*,*) 'common mult = ', common_multiplier(1:100)
        call cpu_time(start)
        do l = 1, matrix_size
            do n = 2 - mod(l, 2), matrix_size, 2
                !  sum from mod(n + 1, 2) because here n = n0 - m + 1 =>
                !  mod(n0-m, 2) = mod(n0 - m + 1 + 1, 2) = mod(n + 1, 2)
                do r = mod(n + 1, 2), min(first%maxd, second%maxd), 2 * bucket
                    !if (n == 1 .and. l == 1) then
                    !    write(*,*) 'add to delta', first%legendre(r, n) * second%legendre(r, l) * &
                    !            common_multiplier(r)
                    !end if
                    do k = r, r + 2 * (bucket - 1), 2
                        value = first%legendre(k, n) * second%legendre(k, l) * &
                                common_multiplier(k)
                        Delta(n, l) = Delta(n, l) + value
                    end do
                    if (abs(value) < accuracy * abs(Delta(n, l))) then
!                         write(*,*) 'n = ', n, 'l = ', l, 'required r = ', r
                         exit
                    end if
                enddo
            enddo
        enddo
        call cpu_time(finish)
!        write(*,*) 'tim = ', finish - start
!        call log_matrix(FILE_DESCRIPTOR(WARNING), 'Delta', Delta, .false., matrix_size)
        !write(*,*)
!        write(*,*) 'Delta60 = ', Delta(:,60)
        deallocate(common_multiplier)

    end subroutine calculate_delta

    subroutine calculate_pi(first, second, Pi, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: Pi(matrix_size, matrix_size)

        real(knd), allocatable, dimension(:) :: common_multiplier
        integer n, l, r, k
        real :: start, finish
        complex(knd) :: value

        Pi = 0q0
        if (first%m /= second%m) then
            write(*, *) 'different m in delta'
            return
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'delta too few elements in first', first%lnum
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if
        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)

        !write(*,*) 'common mult = ', common_multiplier(1:100)
        call cpu_time(start)
        do l = 1, matrix_size
            do n = 2 - mod(l, 2), matrix_size, 2
                !  sum from mod(n + 1, 2) because here n = n0 - m + 1 =>
                !  mod(n0-m, 2) = mod(n0 - m + 1 + 1, 2) = mod(n + 1, 2)
                do r = mod(n + 1, 2), min(first%maxd, second%maxd), 2 * bucket
                    !if (n == 1 .and. l == 1) then
                    !    write(*,*) 'add to delta', first%legendre(r, n) * second%legendre(r, l) * &
                    !            common_multiplier(r)
                    !end if
                    do k = r, r + 2 * (bucket - 1), 2
                        value = first%legendre(k, n) * second%legendre(k, l) * &
                                common_multiplier(k)
                        Pi(n, l) = Pi(n, l) + value
                    end do
                    if (abs(value) < accuracy * abs(Pi(n, l))) then
                        !                         write(*,*) 'n = ', n, 'l = ', l, 'required r = ', r
                        exit
                    end if
                enddo
            enddo
        enddo
        call cpu_time(finish)
        write(*,*) 'tim = ', finish - start
        !call log_matrix('Delta', matrix_size, matrix_size, Delta, .false., matrix_size, matrix_size, FD_INFO)
        !write(*,*)
        deallocate(common_multiplier)

    end subroutine calculate_pi

    real(knd) function gamma_c_lower(m, r)
        integer m, r
        gamma_c_lower = real(r, knd) / (2 * r + 2 * m - 1)
    end function gamma_c_lower
    real(knd) function gamma_c_upper(m, r)
        integer m, r
        gamma_c_upper = real(2 * m + 1 + r, knd) / (2 * r + 2 * m + 3)
    end function gamma_c_upper

    subroutine calculate_gamma(first, second, Gamma, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Gamma(matrix_size, matrix_size)
        integer n, l, r, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in gamma'
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'too few elements in first'
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if

        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)
        Gamma = 0q0

        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 0) then
                    cycle
                endif
                ix = mod(n + 1, 2)
                if (ix == 0) then
                    Gamma(n, l) = first%legendre(ix, n) * second%legendre(ix + 1, l) * gamma_c_upper(first%m, ix) * &
                            common_multiplier(ix)
                    ix = 2
                endif

                !write(*,*) first%maxd, second%maxd
                do r = ix, min(first%maxd, second%maxd) - 1, 2
                    if (update_with_accuracy(Gamma(n, l),first%legendre(r, n) * common_multiplier(r) * &
                            (second%legendre(r - 1, l) * gamma_c_lower(first%m, r) + second%legendre(r + 1, l) * &
                                    gamma_c_upper(first%m, r)),accuracy)) then
                        exit
                    end if
                enddo
            enddo
        enddo

        deallocate(common_multiplier)
    end subroutine calculate_gamma

    real(knd) function kappa_c_lower(m, r)
        integer m, r
        kappa_c_lower = -real(r * (r + m - 1), knd) / (2 * r + 2 * m - 1)
    end function kappa_c_lower
    real(knd) function kappa_c_upper(m, r)
        integer m, r
        kappa_c_upper = real((r + m + 2) * (2 * m + 1 + r), knd) / (2 * r + 2 * m + 3)
    end function kappa_c_upper

    subroutine calculate_kappa(first, second, Kappa, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Kappa(matrix_size, matrix_size)
        integer n, l, r, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in gamma'
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'too few elements in first'
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if
        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)
        Kappa = 0q0

        !write(*,*) 'spheroidal kappa'
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 0) then
                    cycle
                endif
                ix = mod(n + 1, 2)
                if (ix == 0) then
                    Kappa(n, l) = first%legendre(ix, n) * second%legendre(ix + 1, l) * kappa_c_upper(first%m, ix) * &
                            common_multiplier(ix)
                    ix = 2
                endif

                !write(*,*) first%maxd, second%maxd
                do r = ix, min(first%maxd, second%maxd) - 1, 2
!                    write(*,*) 'current_update: first_leg = ', first%legendre(r, n), 'cm = ', common_multiplier(r), &
!                            'second leg = ', second%legendre(r - 1, l), second%legendre(r + 1, l), &
!                            'kappa f = ', kappa_c_lower(first%m, r), kappa_c_upper(first%m, r)
                    if (update_with_accuracy(Kappa(n, l), first%legendre(r, n) * common_multiplier(r) * &
                            (second%legendre(r - 1, l) * kappa_c_lower(first%m, r) + second%legendre(r + 1, l) * &
                                    kappa_c_upper(first%m, r)), accuracy)) then
!                        write(*,*) 'n = ', n, 'l = ', l, 'required r = ', r
                        exit
                    end if
                enddo
                !write(*,*) 'n = ', n, 'l = ', l, 'val = ', Kappa(n,l)
            enddo

        enddo

        deallocate(common_multiplier)
    end subroutine calculate_kappa

    ! Sigma
    real(knd) function sigma_c_lower(m, r)
        integer m, r
        sigma_c_lower = -real((r - 1) * r * (r + m - 1), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function sigma_c_lower
    real(knd) function sigma_c_middle(m, r)
        integer m, r
        sigma_c_middle = real(3 * (r + m) * (r + m + 1) - m * m - 2, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function sigma_c_middle
    real(knd) function sigma_c_upper(m, r)
        integer m, r
        sigma_c_upper = real((r + m + 2) * (r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function sigma_c_upper

    subroutine calculate_sigma(first, second, Sigma, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Sigma(matrix_size, matrix_size)
        integer n, l, r, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in epsilon'
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'too few elements in first'
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if
        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)
        Sigma = 0q0

        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                ix = mod(n + 1, 2)

                Sigma(n, l) = first%legendre(ix, n) * common_multiplier(ix) * &
                        (second%legendre(ix, l) * sigma_c_middle(first%m, ix) + &
                                second%legendre(ix + 2, l) * sigma_c_upper(first%m, ix))

                do r = ix + 2, min(first%maxd, second%maxd) - 2, 2
                    if (update_with_accuracy(Sigma(n, l), first%legendre(r, n) * common_multiplier(r) * &
                            (second%legendre(r - 2, l) * sigma_c_lower(first%m, r) + &
                                    second%legendre(r, l) * sigma_c_middle(first%m, r) + &
                                    second%legendre(r + 2, l) * sigma_c_upper(first%m, r)), accuracy)) then
                        exit
                    end if
                enddo
            enddo
        enddo

        deallocate(common_multiplier)
    end subroutine calculate_sigma

    real(knd) function epsilon_c_lower(m, r)
        integer m, r
        epsilon_c_lower = -real((r - 1) * r * (r + m - 2), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function epsilon_c_lower
    real(knd) function epsilon_c_middle(m, r)
        integer m, r
        epsilon_c_middle = real((r + m) * (r + m + 1) - 3 * m * m, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function epsilon_c_middle
    real(knd) function epsilon_c_upper(m, r)
        integer m, r
        epsilon_c_upper = real((r + m + 3) * (r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function epsilon_c_upper

    subroutine calculate_epsilon(first, second, Epsilon, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Epsilon(matrix_size, matrix_size)
        integer n, l, r, matrix_size, ix

        if (first%m /= second%m) then
            write(*, *) 'different m in epsilon'
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'too few elements in first'
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if
        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)
        Epsilon = 0q0

        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                ix = mod(n + 1, 2)

                Epsilon(n, l) = first%legendre(ix, n) * common_multiplier(ix) * &
                        (second%legendre(ix, l) * epsilon_c_middle(first%m, ix) + &
                                second%legendre(ix + 2, l) * epsilon_c_upper(first%m, ix))

                do r = ix + 2, min(first%maxd, second%maxd) - 2, 2
                    if (update_with_accuracy(Epsilon(n, l), first%legendre(r, n) * common_multiplier(r) * &
                            (second%legendre(r - 2, l) * epsilon_c_lower(first%m, r) + &
                                    second%legendre(r, l) * epsilon_c_middle(first%m, r) + &
                                    second%legendre(r + 2, l) * epsilon_c_upper(first%m, r)), accuracy)) then
                        exit
                    end if
                enddo
            enddo
        enddo

        deallocate(common_multiplier)
    end subroutine calculate_epsilon


    ! Omega
    real(knd) function omega_c_lower(m, r)
        integer m, r
        omega_c_lower = -real((r - 1) * r, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function omega_c_lower
    real(knd) function omega_c_middle(m, r)
        integer m, r
        omega_c_middle = real(2 * ((r + m) * (r + m + 1) + m * m - 1), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function omega_c_middle
    real(knd) function omega_c_upper(m, r)
        integer m, r
        omega_c_upper = -real((r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function omega_c_upper

    subroutine calculate_omega(first, second, Omega, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: Omega(matrix_size, matrix_size)

        integer :: ix, n, l, r
        if (first%m /= second%m) then
            write(*, *) 'different m in epsilon'
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'too few elements in first m lnum = ', first%m, first%lnum
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if
        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)
        Omega = 0q0

        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                ix = mod(n + 1, 2)

                Omega(n, l) = first%legendre(ix, n) * common_multiplier(ix) * &
                        (second%legendre(ix, l) * omega_c_middle(first%m, ix) + &
                                second%legendre(ix + 2, l) * omega_c_upper(first%m, ix))

                do r = ix + 2, min(first%maxd, second%maxd) - 2, 2
                    if (update_with_accuracy(Omega(n, l), first%legendre(r, n) * common_multiplier(r) * &
                            (second%legendre(r - 2, l) * omega_c_lower(first%m, r) + &
                                    second%legendre(r, l) * omega_c_middle(first%m, r) + &
                                    second%legendre(r + 2, l) * omega_c_upper(first%m, r)), accuracy)) then
                        exit
                    end if
                enddo
            enddo
        enddo
        deallocate(common_multiplier)
    end subroutine calculate_omega

    real(knd) function tau_c_middle(m, r)
        integer m, r
        tau_c_middle = real((r + m) * (r + m + 1), knd)
    end function tau_c_middle

    subroutine calculate_tau(first, second, Tau, matrix_size)
        class(SpheroidalCalculation), intent(in) :: first, second
        real(knd), allocatable, dimension(:) :: common_multiplier
        complex(knd), intent(out) :: Tau(matrix_size, matrix_size)
        integer n, l, r, matrix_size

        if (first%m /= second%m) then
            write(*, *) 'different m in delta'
        end if
        if (first%lnum < matrix_size) then
            write(*, *) 'too few elements in first'
            return
        end if
        if (second%lnum < matrix_size) then
            write(*, *) 'too few elements in second'
            return
        end if
        !write(*,*) 'in tau m = ', first%m, 'lnum = ', first%lnum, second%lnum, 'ms = ', matrix_size
        !write(*,*) 'size = ', size(Tau), size(Tau(1,:)), size(Tau(:, 1))
        call fill_common_multiplier(first%m, MAXIMUM_D_COEFF_NUMBER, common_multiplier)
        Tau = 0q0

        !write(*,*) 'spheroidal', first%c, second%c
        do n = 1, matrix_size
            do l = 1, matrix_size
                if (mod(abs(n - l), 2) == 1) then
                    cycle
                endif
                do r = mod(n + 1, 2), min(first%maxd, second%maxd), 2
                    if (update_with_accuracy(Tau(n, l), first%legendre(r, n) * second%legendre(r, l) * &
                            common_multiplier(r) * tau_c_middle(first%m, r), accuracy)) then
                        exit
                    end if
                enddo
                !write(*,*) 'n = ', n, 'l = ', l, 'val = ', Tau(n, l)
            enddo
        enddo
        !write(*,*)

        deallocate(common_multiplier)
    end subroutine calculate_tau

end module spheroidal_integrals