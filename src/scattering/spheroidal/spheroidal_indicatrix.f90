module spheroidal_indicatrix
    use constants
    use spheroidal
    use regime
    use spheroidal_scatterer
    use wavelength_point
    use angle
    implicit none

contains
    !  prints the scattering indicatrix to a given file descriptor
    !  iteration over
    !    - buckets of theta with size BUCKET_SIZE from the module constants
    !      (the last bucket can be smaller depending on if BUCKET_SIZE does not
    !      divide ntheta
    ! inside get_scattering_matrix_for_bucket
    !    - m from minm to maxm using solution_te and solution_tm  obtained when
    !      the qfactors were calculated
    !    - theta inside the current bucket
    !    - phi
    subroutine print_scattering_indicatrix(scatterer, calculation_point, minm, maxm, &
            matrix_size, solution_te, solution_tm, &
            ntheta, thetas, nphi, phis, fd)
        type(SpheroidalScatterer), intent(in) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        integer, intent(in) :: minm, maxm, matrix_size, ntheta, nphi, fd
        complex(knd), intent(in) :: solution_te(2 * matrix_size, minm:maxm), solution_tm(2 * matrix_size, minm:maxm)

        type(SpheroidalCalculation) :: direction_calculation
        real(knd) :: dtheta, dphi
        real(knd) :: arg(BUCKET_SIZE)
        integer :: bucket_start, bucket_end, bucket_length, i, m, j
        type(AngleType) :: thetas(ntheta + 1), phis(nphi + 1, minm:maxm)
        real(knd) :: bucket(4,4,nphi + 1,BUCKET_SIZE)

        do bucket_start = 1, ntheta + 1, BUCKET_SIZE
            bucket_end = min(ntheta + 1, bucket_start + BUCKET_SIZE - 1)
            bucket_length = bucket_end - bucket_start + 1
            bucket(:,:,:,1:bucket_length) = get_scatterring_matrix_for_bucket(scatterer, calculation_point, &
                    minm, maxm, matrix_size, &
            solution_te, solution_tm, bucket_length, thetas(bucket_start:bucket_end), nphi + 1, phis)
            !write(*,*) 'bucket = ', bucket(:,:,:,1:bucket_length)
            do i = 1, bucket_length
                do j = 1, nphi + 1
                    write(fd, *) thetas(bucket_start + i - 1)%value * 180q0 / PI, &
                            phis(j, min(1, maxm))%value * 180q0 / PI, bucket(:,:,j,i)
                end do
            end do
        end do

    end subroutine print_scattering_indicatrix

    function get_scatterring_matrix_for_bucket(scatterer, calculation_point, minm, maxm, &
            matrix_size, solution_te, solution_tm, &
            ntheta, thetas, nphi, phis) result(scat)
        type(SpheroidalScatterer) :: scatterer
        type(WavelengthPoint), intent(in) :: calculation_point
        integer, intent(in) :: minm, maxm, matrix_size, ntheta, nphi
        complex(knd), intent(in) :: solution_te(2 * matrix_size, minm:maxm), solution_tm(2 * matrix_size, minm:maxm)
        type(AngleType), intent(in) :: thetas(ntheta), phis(nphi, minm:maxm)

        type(SpheroidalCalculation) :: direction_calculation
        real(knd) :: arg(ntheta)

        real(knd) :: scat(4,4,nphi,ntheta)
        complex(knd) :: ampl(2,2,nphi,ntheta), part(2,2)
        integer :: m, i, j

        scat = 0
        ampl = 0
        arg = (/(thetas(i)%angle_cos, i = 1, ntheta)/)
        do m = minm, maxm
            call direction_calculation%calculate(m, matrix_size, scatterer%c0(1), scatterer%ksi(1), &
                    ntheta, arg, scatterer%spheroidal_type, .false., .true., .false.)
            do i = 1, ntheta
                do j = 1, nphi
                    part = 0
                    part = get_nonaxisymmetric_amplitude_matrix(&
                            matrix_size, solution_te(:,m), solution_tm(:,m), calculation_point, &
                            thetas(i), phis(j, m), direction_calculation, i)
                    ampl(:,:,j,i) = ampl(:,:,j,i) + part
                    write(*,*) 'm = ', m, ' part = ', part, ' dpart = ', part / ampl(:,:,j,i)
                end do
            end do
        end do
        do i = 1, ntheta
            do j = 1, nphi
                scat(:,:,j,i) = get_scattering_matrix(ampl(:,:,j,i))
            end do
        end do
        !write(*,*) 'scat = ', scat
    end function get_scatterring_matrix_for_bucket

    complex(knd) function get_axisymmetric_amplitude_te_1_1(matrix_size, solution, direction_calculation, argnum) result(T11)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(matrix_size)

        integer :: i
        !write(*, *) 'get_axisymmetric_amplitude_1_1', ' m = ', this%m
        T11 = 0

        do i = 1, matrix_size
            T11 = T11 + NEGIDEG(mod(i, 4)) * solution(i) * direction_calculation%s1(i, argnum)
        end do

    end function get_axisymmetric_amplitude_te_1_1

    complex(knd) function get_axisymmetric_amplitude_tm_2_2(matrix_size, solution, direction_calculation, argnum) result(T22)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(matrix_size)

        !write(*, *) 'get_axisymmetric_amplitude_2_2', ' m = ', this%m
        T22 = -get_axisymmetric_amplitude_te_1_1(matrix_size, solution, direction_calculation, argnum)

    end function get_axisymmetric_amplitude_tm_2_2

    complex(knd) function get_nonaxisymmetric_amplitude_1_1(matrix_size, solution, calculation_point, &
            theta, mphi, direction_calculation, argnum) result(T11)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)

        integer :: i

        !write(*, *) 'get_nonaxisymmetric_amplitude_1_1', ' m = ', this%m

        T11 = 0

        do i = 1, matrix_size
            T11 = T11 + NEGIDEG(mod(i + direction_calculation%m + 2, 4)) * &
                    (calculation_point%k * solution(i) * direction_calculation%s1(i, argnum) + &
                            IDEG(1) * solution(i + matrix_size) * direction_calculation%s1d(i, argnum))
        end do

        T11 = T11 * theta%angle_sin * mphi%angle_cos

    end function get_nonaxisymmetric_amplitude_1_1

    complex(knd) function get_nonaxisymmetric_amplitude_1_2(matrix_size, solution, &
            theta, mphi, direction_calculation, argnum) result(T12)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)

        integer :: i
        !write(*,*) 'get_nonaxisymmetric_amplitude_1_2', ' m = ', this%m

        T12 = 0
        do i = 1, matrix_size
            T12 = T12 - NEGIDEG(mod(i + direction_calculation%m + 3, 4)) * solution(i + matrix_size) * &
                    direction_calculation%m * direction_calculation%s1(i, argnum)
        end do

        T12 = T12 / theta%angle_sin * mphi%angle_sin

    end function get_nonaxisymmetric_amplitude_1_2

    complex(knd) function get_nonaxisymmetric_amplitude_2_1(matrix_size, solution, &
            theta, mphi, direction_calculation, argnum) result(T21)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)

        !write(*, *) 'get_nonaxisymmetric_amplitude_2_1', ' m = ', this%m
        T21 = get_nonaxisymmetric_amplitude_1_2(matrix_size, solution, &
                theta, mphi, direction_calculation, argnum)

    end function get_nonaxisymmetric_amplitude_2_1

    complex(knd) function get_nonaxisymmetric_amplitude_2_2(matrix_size, solution, calculation_point, &
            theta, mphi, direction_calculation, argnum) result(T22)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)
        !write(*, *) 'get_nonaxisymmetric_amplitude_2_2', ' m = ', this%m

        T22 = -get_nonaxisymmetric_amplitude_1_1(matrix_size, solution, calculation_point, &
                theta, mphi, direction_calculation, argnum)

    end function get_nonaxisymmetric_amplitude_2_2

    function get_nonaxisymmetric_amplitude_matrix(matrix_size, solution_te, solution_tm, calculation_point, &
            theta, mphi, direction_calculation, argnum) result(ampl)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution_te(2 * matrix_size), solution_tm(2 * matrix_size)

        complex(knd) :: ampl(2,2)

        ampl(1,1) = get_nonaxisymmetric_amplitude_1_1(matrix_size, solution_tm, calculation_point, &
                theta, mphi, direction_calculation, argnum)
        ampl(1,2)= get_nonaxisymmetric_amplitude_1_2(matrix_size, solution_te, &
                theta, mphi, direction_calculation, argnum)
        ampl(2,1) = get_nonaxisymmetric_amplitude_2_1(matrix_size, solution_tm, &
                theta, mphi, direction_calculation, argnum)
        ampl(2,2) = get_nonaxisymmetric_amplitude_2_2(matrix_size, solution_te, calculation_point, &
                theta, mphi, direction_calculation, argnum)
        if (direction_calculation%m == 0) then
            ampl = ampl * 0.5_knd
        end if
        !write(*,*) 'theta = ', theta%value, 'phi = ', mphi%value, 'ampl = ', ampl

    end function get_nonaxisymmetric_amplitude_matrix

    function get_scattering_matrix(ampl) result(scat)
        complex(knd), intent(in) :: ampl(2,2)

        real(knd) :: scat(4,4)

        write(*,*) 'ampl0 = ', ampl
        scat(1,1) = (abs(ampl(1,1))**2 + abs(ampl(1,2))**2 + abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(1,2) = (-abs(ampl(1,1))**2 - abs(ampl(1,2))**2 + abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(1,3) = real(ampl(2,2) * conjg(ampl(1,2)) + ampl(1,1) * conjg(ampl(2,1)),knd)
        scat(1,4) = imag(ampl(2,2) * conjg(ampl(1,2)) - ampl(1,1) * conjg(ampl(2,1)))
        scat(2,1) = (-abs(ampl(1,1))**2 + abs(ampl(1,2))**2 - abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(2,2) = (abs(ampl(1,1))**2 - abs(ampl(1,2))**2 - abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(2,3) = real(ampl(2,2) * conjg(ampl(1,2)) - ampl(1,1) * conjg(ampl(2,1)),knd)
        scat(2,4) = imag(ampl(2,2) * conjg(ampl(1,2)) + ampl(1,1) * conjg(ampl(2,1)))
        scat(3,1) = real(ampl(2,2) * conjg(ampl(2,1)) + ampl(1,1) * conjg(ampl(1,2)),knd)
        scat(3,2) = real(ampl(2,2) * conjg(ampl(2,1)) - ampl(1,1) * conjg(ampl(1,2)),knd)
        scat(3,3) = real(ampl(1,1) * conjg(ampl(2,2)) + ampl(1,2) * conjg(ampl(2,1)),knd)
        scat(3,4) = imag(ampl(2,2) * conjg(ampl(1,1)) + ampl(2,1) * conjg(ampl(1,2)))
        scat(4,1) = imag(ampl(2,1) * conjg(ampl(2,2)) + ampl(1,1) * conjg(ampl(1,2)))
        scat(4,2) = imag(ampl(2,1) * conjg(ampl(2,2)) - ampl(1,1) * conjg(ampl(1,2)))
        scat(4,3) = imag(ampl(1,1) * conjg(ampl(2,2)) - ampl(1,2) * conjg(ampl(2,1)))
        scat(4,4) = real(ampl(1,1) * conjg(ampl(2,2)) - ampl(1,2) * conjg(ampl(2,1)),knd)
        write(*,*) 'ampl = ', ampl
        write(*,*) 'scat = ', scat

    end function get_scattering_matrix

end module spheroidal_indicatrix