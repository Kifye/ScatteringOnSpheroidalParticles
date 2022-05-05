! Created by drakosha on 15.07.2021.

module spheroidal_initial
    use spheroidal
    use wavelength_point
    use regime
    use angle
    use constants
    use logging
    implicit none
contains
    ! Initializers
    subroutine set_axisymmetric_initial_te(outside_layer, matrix_size, initial, argnum)
        type(SpheroidalCalculation), intent(in) :: outside_layer
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(matrix_size)
        integer, intent(in) :: argnum

        integer :: i

        do i = 1, matrix_size
            initial(i) = -2q0 * IDEG(mod(i, 4)) * outside_layer%s1(i, argnum)
        enddo

        !call log_array('outside_layer%s1', matrix_size, outside_layer%s1(:, argnum), 10, FD_INFO)
        !write(*,*) 'lnum = ', outside_layer%lnum
        !call log_array('symmetric_initial_te', matrix_size, initial, 10, FD_INFO)
    end subroutine set_axisymmetric_initial_te

    subroutine set_nonaxisymmetric_initial_te(outside_layer, calculation_point, alpha, matrix_size, initial, argnum)
        type(SpheroidalCalculation), intent(in) :: outside_layer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: alpha
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(2 * matrix_size)
        integer, intent(in) :: argnum

        integer :: i, j
        real(knd) :: coefficient

        initial = 0

        if (abs(alpha%value) > 1q-6) then
            coefficient = outside_layer%spheroidal_type * 4.0q0 / (calculation_point%k * alpha%angle_sin)
            do i = 1, matrix_size
                initial(i) = coefficient * IDEG(mod(i + outside_layer%m + 2, 4)) * outside_layer%s1(i, argnum)
            enddo
        elseif (outside_layer%m == 1) then
            coefficient = outside_layer%spheroidal_type * 4.0q0 / calculation_point%k
            do i = 1, matrix_size
                do j = 1 - mod(i, 2), outside_layer%maxd, 2
                    initial(i) = initial(i) + (j + 1) * (j + 2) / 2q0 * outside_layer%legendre(j, i)
                end do
                initial(i) = coefficient * IDEG(mod(i + outside_layer%m + 2, 4)) * initial(i)
            end do
        endif

    end subroutine set_nonaxisymmetric_initial_te

    subroutine set_axisymmetric_initial_tm(outside_layer, calculation_point, matrix_size, initial, argnum)
        type(SpheroidalCalculation), intent(in) :: outside_layer
        type(WavelengthPoint), intent(in) :: calculation_point
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(matrix_size)
        integer, intent(in) :: argnum

        integer :: i
        real(knd) :: eps_to_mu

        eps_to_mu = sqrt(calculation_point%eps(0) / calculation_point%mu(0))

        call set_axisymmetric_initial_te(outside_layer, matrix_size, initial, argnum)
        initial = -initial * eps_to_mu

    end subroutine set_axisymmetric_initial_tm

    subroutine set_nonaxisymmetric_initial_tm(outside_layer, calculation_point, alpha, matrix_size, initial, argnum)
        type(SpheroidalCalculation), intent(in) :: outside_layer
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: alpha
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: initial(2 * matrix_size)
        integer, intent(in) :: argnum

        integer :: i
        real(knd) :: eps_to_mu

        eps_to_mu = sqrt(calculation_point%eps(0) / calculation_point%mu(0))

        call set_nonaxisymmetric_initial_te(outside_layer, calculation_point, alpha, matrix_size, initial, argnum)
        initial = -outside_layer%spheroidal_type * initial * eps_to_mu

    end subroutine set_nonaxisymmetric_initial_tm

end module spheroidal_initial