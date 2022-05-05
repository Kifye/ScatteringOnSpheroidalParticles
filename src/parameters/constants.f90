!  contains general constants and functions for obtaining parameter values
!  like size of matrices and number of legendre coefficients
module constants
    use regime

    implicit none

    real(knd), parameter :: PI = 4q0 * atan(1q0)
    real(knd), parameter :: BASE_ACCURACY = 1q-16
    real(knd), parameter :: MIN_M_RATIO = 1q-16
    real(knd), parameter :: DEFAULT_PAIR_INTEGRAL_ACCURACY = 1q-24;
    real(4), parameter :: REVERSE_MATRIX_SIZE_RATIO = 2.0 !116.0 / 102.0
    integer, parameter :: MAXIMUM_MATRIX_SIZE = 400
    integer, parameter :: MAXIMUM_D_COEFF_NUMBER = 100000
    integer, parameter :: BUCKET_SIZE = 1000
    complex(knd), parameter :: IDEG(0:3) = (/ &
            cmplx(1q0, 0q0, knd), cmplx(0q0, 1q0, knd), &
            cmplx(-1q0, 0q0, knd), cmplx(0q0, -1q0, knd)/)
    complex(knd), parameter :: NEGIDEG(0:3) = (/ &
            cmplx(1q0, 0q0, knd), cmplx(0q0, -1q0, knd), &
                    cmplx(-1q0, 0q0, knd), cmplx(0q0, 1q0, knd)/)

    enum, bind(c)
        ! Spheroidal type
        enumerator :: OBLATE = -1, PROLATE = 1
        !  levels of logs
        !  DO NOT CHANGE THE ORDER! PARAMETER ARRAYS FOR LOGS ARE DEFINED AS
        !  ARRAY(ERROR:DETAIL)
        enumerator :: ERROR, WARNING, DEBUG, INFO, DETAIL
        ! Scattering modes
        enumerator TE, TM
        ! Types of scattering factors
        enumerator :: qfactors, cfactors, normalized_cfactors
        enumerator :: no_absorbtion, consequential
    end enum
contains

    integer function get_full_matrix_size(matrix_size)
        integer, intent(in) :: matrix_size
        get_full_matrix_size = aint(REVERSE_MATRIX_SIZE_RATIO * matrix_size)
    end function get_full_matrix_size

    integer function get_full_function_size(matrix_size)
        integer, intent(in) :: matrix_size
        get_full_function_size = aint(REVERSE_MATRIX_SIZE_RATIO * matrix_size)
    end function get_full_function_size

end module constants