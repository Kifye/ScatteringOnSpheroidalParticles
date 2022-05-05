!  Wrapper for subroutines for spheroidal function calculation by Anri van Buren
!  from modules vb_prolate and vb_oblate
!  the kind parameter for complex function values is set in module regime as knd
module spheroidal

    use regime
    use constants
    use logging

    use complex_prolate_swf
    use complex_oblate_swf
    !use vb_prolate
    !use vb_oblate
    use param

    implicit none
    private
    public :: get_ksi, get_c, theoretical_wronskian_radial

    !  Represents the spheroidal function
    !  contains the tabulated radial functions of 1st, 2nd, 3rd kind and their derivatives, R^(i)_{m, n}(c, \ksi)
    !  angular functions of the 1st kind and it's derivative, S^(1)_{m, n}(c, \eta)
    !  and legendre coefficients for the angular function
    !  normalized with Mexxie-Shafke narmalization scheme to 1
    !  its calculate method calls subroutines calculating spheroidal functions by Anri van Buren:
    !  cprofcn (module vb_prolate) and coblfcn (module vb_oblate)
    type, public :: SpheroidalCalculation
        !  spheroidal_type = 1 for prolate s.f. and -1 for oblate s.f.
        integer :: spheroidal_type
        ! spheroidal function parameter m
        integer :: m
        !  lnum = n - m + 1 >= 0
        integer :: lnum
        ! spheroidal function parameter c
        complex(knd) :: c

        !  ksi = \ksi - 1 for prolate functions and
        !  ksi = \ksi for oblate functions
        !  so ksi \in [0..)
        real(knd) :: ksi
        ! radial functions and their derivatives for n - m + 1 = 1,.., lnum, arrays [1:lnum]
        complex(knd), allocatable, dimension(:) :: r1, r1d, r2, r2d, r3, r3d

        !  number of values of \eta for angular function
        integer :: narg
        !  array[1:narg] with parameters \eta of angular function
        real(knd), allocatable, dimension(:) :: arg
        !  angular function and its derivative for \eta \in arg, n - m = 0,.., lnum - 1
        !  arrays[1:lnum][1:narg]
        complex(knd), allocatable, dimension(:, :) :: s1, s1d

        !  number of legendre coefficients calculated (including zeros)
        integer :: maxd
        !  legendre coefficients for angular functions, arrays [0:maxd][1:lnum]
        complex(knd), allocatable, dimension(:, :) :: legendre

        logical :: functions_calculated
        logical :: legendre_coefficients_calculated
    contains
        !  Preparation
        !  (re)allocates arrays of function values to required sizes if necessary
        procedure, private :: check_and_reallocate_arrays

        !  Calculation
        !  calculates spheroidal functions with giver parameters
        procedure, public :: calculate

        !  get array of wronskian values
        procedure :: get_wronskian

        !  Free memory if it was allocated (noexception)
        procedure :: delete_radial_functions, delete_angular_functions, delete_legendre_coefficients, &
                delete_calculation_full

        !  destructor, calls delete_calculation_full
        final :: delete_calculation
    end type SpheroidalCalculation

    character(len=MESSAGE_LENGTH), parameter :: spheroidal_function_log_format = &
            '(a,a,a,a,i4,a,i4,a,2F15.7,a,f15.7,a,i4,a,10000f15.7)'
    integer, parameter :: fd_spheroidal_calculation = FILE_DESCRIPTOR(INFO)
    integer, parameter :: fd_spheroidal_calculation_detail = FILE_DESCRIPTOR(DETAIL)

contains
    subroutine set(this, m, n, c, ksi, narg, arg, f)
        class(SpheroidalCalculation) :: this
        integer, intent(in) :: m, n, narg
        real(knd), intent(in) :: ksi, arg(narg)
        complex(knd), intent(in) :: c
        integer, optional, intent(in) :: f

        if (present(f)) then
            this%spheroidal_type = f
        else
            this%spheroidal_type = 1
        endif

        this%m = m
        this%c = c
        if (this%spheroidal_type == 1) then
            this%ksi = ksi - 1q0
        else
            this%ksi = ksi
        endif
        this%narg = narg
        this%lnum = n - m + 1

        if (.not.(allocated(this%r1) .and. (n - m + 1 == this%lnum))) then
            !write(*,*) 'allocating spheroidal calculation'
            if (allocated(this%r1)) then
                call delete_calculation(this)
            endif

            this%lnum = n - m + 1

            allocate(this%r1(this%lnum), this%r1d(this%lnum), &
                    this%r2(this%lnum), this%r2d(this%lnum), this%r3(this%lnum), this%r3d(this%lnum))
            allocate(this%s1(this%lnum, this%narg), this%s1d(this%lnum, this%narg))
            allocate(this%arg(narg))
        endif
        !allocate(this%arg(this%narg))
        this%arg = arg
        this%functions_calculated = .false.

    end subroutine set


    ! checks if the arrays for function values and angular function argument values need to be reallocated
    ! assignes and reallocates them only if necessary
    subroutine check_and_reallocate_arrays(this, lnum, narg, allocate_radial, allocate_angular)
        class(SpheroidalCalculation), intent(inout) :: this
        integer, intent(in) :: lnum, narg
        logical, intent(in) :: allocate_radial, allocate_angular

        if (allocate_radial .and. allocated(this%r1) .and. lnum /= this%lnum) then
            deallocate(this%r1, this%r1d, this%r2, this%r2d, this%r3, this%r3d, this%s1, this%s1d)
        end if

        if (allocated(this%arg) .and. narg /= this%narg) then
            deallocate(this%arg)
        end if

        if (allocate_angular .and. allocated(this%s1) .and. (narg /= this%narg .or. lnum /= this%lnum)) then
            deallocate(this%s1, this%s1d)
        end if

        this%lnum = lnum
        this%narg = narg
        if (allocate_radial .and. .not. allocated(this%r1)) then
            allocate(this%r1(lnum), this%r1d(lnum), this%r2(lnum), this%r2d(lnum), this%r3(lnum), this%r3d(lnum))
        end if

        if (.not. allocated(this%arg)) then
            allocate(this%arg(narg))
        end if
        if (allocate_angular .and. .not. allocated(this%s1)) then
            allocate(this%s1(lnum, narg), this%s1d(lnum, narg))
        end if
    end subroutine check_and_reallocate_arrays

    subroutine check_expected_accuracy(lnum, narg, naccr, naccs, naccds)
        integer, intent(in) :: lnum, narg, naccr(lnum), naccs(lnum, narg), naccds(lnum, narg)
        character(len = MESSAGE_LENGTH) :: log
        integer :: i, j

        do i = 1, lnum
            if (naccr(i) < EXPECTED_FUNCTION_ACCURACY) then
                write(log, *) 'Radial spheroidal function value for n - m + 1 = ', i, ' only has accuracy ', naccr(i)
                call log_message(log, WARNING, FILE_DESCRIPTOR(WARNING))
            end if
        end do
        do j = 1, narg
            do i = 1, lnum
                if (naccs(i, j) < EXPECTED_FUNCTION_ACCURACY) then
                    write(log, *) 'Angular spheroidal function value for argnum ', j, ', n - m + 1 = ', i, &
                            ' only has accuracy ', naccs(i, j)
                    call log_message(log, WARNING, FILE_DESCRIPTOR(WARNING))
                end if
                if (naccds(i, j) < EXPECTED_FUNCTION_ACCURACY) then
                    write(log, *) 'Angular spheroidal function derivative value for argnum ', j, ', n - m + 1 = ', i, &
                            ' only has accuracy ', naccds(i, j)
                    call log_message(log, WARNING, FILE_DESCRIPTOR(WARNING))
                end if
            end do
        end do
    end subroutine check_expected_accuracy

    subroutine calculate(this, m, lnum, c, ksi, narg, arg, f, &
            need_radial_functions, need_angular_functions, need_legendre_coefficients)
        class(SpheroidalCalculation), intent(inout) :: this
        integer, intent(in) :: m, lnum, narg, f
        real(knd), intent(in) :: ksi, arg(narg)
        complex(knd), intent(in) :: c
        logical, optional, intent(in) :: need_radial_functions, need_angular_functions, need_legendre_coefficients
        integer :: legendre_exp(0:MAXIMUM_D_COEFF_NUMBER)
        complex(knd) :: value
        real(knd) :: rad

        !  cprofcn and coblfcn calculate exponent and mantissa in different arrays
        !  these arrays are used to store them, their types have kind = vb_kind
        !  the arrays from SpheroidalCalculation contain the full value
        complex(vb_knd), allocatable, dimension(:) :: r1, r1d, r2, r2d
        integer, allocatable, dimension(:) :: r1_exp, r1d_exp, r2_exp, r2d_exp
        complex(vb_knd), allocatable, dimension(:, :) :: s1, s1d
        integer, allocatable, dimension(:, :) :: s1_exp, s1d_exp
        !  legendre coeffictient ratios,
        !  the array is allocated inside cprofcn, coblfcn and deallocated right after
        !  this%legendre is initialized from it
        integer :: maxd, mmaxval
        complex(vb_knd), allocatable, dimension(:, :) :: enr
        !  expected accuracies for function values
        integer, allocatable, dimension(:) :: naccr
        integer, allocatable, dimension(:, :) :: naccs, naccds
        real(vb_knd) :: vb_arg(narg)

        integer :: i, j, mode_radial, mode_angular
        logical :: need_radial, need_angular, need_legendre

        !  process optional parameters. used for indicatrix calculation
        !  where only angular function values are needed
        need_radial = .true.
        need_angular = .true.
        need_legendre = .true.
        if (present(need_radial_functions)) then
            need_radial = need_radial_functions
        end if
        if (present(need_angular_functions)) then
            need_angular = need_angular_functions
        end if
        if (present(need_legendre_coefficients)) then
            need_legendre = need_legendre_coefficients
        end if
        if (need_radial) then
            mode_radial = 2
        else
            mode_radial = 0
        end if
        if (need_angular) then
            mode_angular = 2
        else
            mode_angular = 0
        end if

        allocate(r1(lnum), r1_exp(lnum), r1d(lnum), r1d_exp(lnum), &
                r2(lnum), r2_exp(lnum), r2d(lnum), r2d_exp(lnum))
        r1 = 0; r1d = 0; r2 = 0; r2d = 0; r1_exp = 0; r1d_exp = 0; r2_exp = 0; r2d_exp = 0

        allocate(s1(lnum, narg), s1_exp(lnum, narg), s1d(lnum, narg), s1d_exp(lnum, narg))
        s1 = 0; s1d = 0; s1_exp = 0; s1d_exp = 0

        allocate(naccr(lnum), naccs(lnum, narg), naccds(lnum, narg))

        call log_message('Start spheroidal calculation', INFO, FILE_DESCRIPTOR(INFO))
        call log_spheroidal_calculation(f, m, lnum, c, ksi, narg, arg)

        !  call functions by van Buren
        !  use different arrays to support different kind parameters for van Buren's
        !  procedures and our procedures
        if (f == 1) then
            call cprofcn_new(cmplx(c%re, c%im, vb_knd), m, lnum, mode_radial, real(ksi, vb_knd) - 1, &
                    mode_angular, 1, narg, real(arg, vb_knd), &
                    r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp, naccr, &
                    s1, s1_exp, s1d, s1d_exp, naccs, naccds, need_legendre, maxd, enr)
!            call cprofcn(cmplx(c%re, c%im, vb_knd), m, lnum, 2, real(ksi, vb_knd) - 1, &
!                    r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp, &
!                    2, 1, narg, real(arg, vb_knd), s1, s1_exp, s1d, s1d_exp, enr, maxd)
        else
            call coblfcn_new(cmplx(c%re, c%im, vb_knd), m, lnum, mode_radial, real(ksi, vb_knd), &
                    mode_angular, 1, narg, real(arg, vb_knd), &
                    r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp, naccr, &
                    s1, s1_exp, s1d, s1d_exp, naccs, naccds, need_legendre, maxd, enr)
!            call coblfcn(c, m, lnum, 2, ksi, &
!                    r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp, &
!                    2, 1, narg, arg, s1, s1_exp, s1d, s1d_exp, enr, maxd)
        end if

        !write(*,*) 'enr1 = ', enr(1:1000, 1)
        !write(*,*) 's1_base = ', s1(1:10, 1)

        !  set arguments
        this%spheroidal_type = f
        this%m = m
        this%c = c
        this%ksi = ksi
        call this%check_and_reallocate_arrays(lnum, narg, need_radial, need_angular)
        this%arg = arg

        !write(*,*) 'r1 = ', r1(1)
        !  set function values
        if (need_radial) then
            this%r1 = 0
            this%r1d = 0
            this%r2 = 0
            this%r2d = 0
            this%r3 = 0
            this%r3d = 0

            this%r1 = r1 * (10q0**r1_exp)
            !write(*,*) 'this%r1 = ', this%r1
            this%r1d = r1d * (10q0**r1d_exp)
            this%r2 = r2 * (10q0**r2_exp)
            this%r2d = r2d * (10q0**r2d_exp)
            this%r3 = this%r1 + cmplx(0q0, 1q0, knd) * this%r2
            this%r3d = this%r1d + cmplx(0q0, 1q0, knd) * this%r2d
        end if
        if (need_angular) then
            this%s1 = 0
            this%s1d = 0
            this%s1 = s1 * (10q0**s1_exp)
            this%s1d = s1d * (10q0**s1d_exp)
        end if
        this%functions_calculated = .true.

        !call check_expected_accuracy(lnum, narg, naccr, naccs, naccds)
        deallocate(r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp)
        deallocate(s1, s1_exp, s1d, s1d_exp)
        deallocate(naccr, naccs, naccds)
        call log_message('Calculated spheroidal functions', INFO, FILE_DESCRIPTOR(INFO))

        !  set legendre coefficients

        if (need_legendre) then
            !  the enr array obtained from cprofcn and coblfcn has the ratios of
            !  legendre coefficients enr(i) = legendre(2 * i + 2 + ix) / legendre(2 * i + ix)
            !  ix = (n - m) mod 2
            if (allocated(this%legendre)) then
                deallocate(this%legendre)
            endif
            !  set legendre
            this%maxd = min(2 * maxd + 1, MAXIMUM_D_COEFF_NUMBER)
            allocate(this%legendre(0:this%maxd, this%lnum))
            this%legendre = 0
            rad = radix(real(enr(1,1)))
!            write(*,*) 'rad = ', rad
            mmaxval = 0
            do i = 1, this%lnum
                legendre_exp = 0
                this%legendre(mod(i + 1, 2), i) = 1.0q0
                !write(*,*) 'i = ', i, 'enr = ', enr(:100,i)
                do j = mod(abs(i + 1), 2) + 2, this%maxd, 2
                    value = this%legendre(j - 2, i) * enr(j / 2, i)
!                    write(*,*) 'real_exp = ', exponent(real(value, knd)), 'imag_exp = ', exponent(imag(value))
                    legendre_exp(j) = exponent(real(value, knd)) + legendre_exp(j - 2)
                    this%legendre(j, i) = value / rad ** exponent(real(value, knd))
                    mmaxval = max(mmaxval, legendre_exp(j))
                    !                this%legendre(j, i) = this%legendre(j - 2, i) * enr(i, j / 2) ! for vb_...
                enddo
!                write(*,*) 'maxval = ', mmaxval
!                write(*,*) 'legexp = ', legendre_exp(:10)
!                legendre_exp(:this%maxd) = legendre_exp(:this%maxd) - mmaxval
!                write(*,*) 'legexp = ', legendre_exp(:10)
                do j = mod(abs(i + 1), 2), this%maxd, 2
                    legendre_exp(j) = legendre_exp(j) - mmaxval
                this%legendre(j,i) = this%legendre(j,i) * rad ** legendre_exp(j)

                    end do
            enddo

            deallocate(enr)

            !  normalize legendre
            do i = 1, this%lnum
                call normalize(this%legendre(:, i), this%maxd, this%m)
            enddo
!            write(*,*) 'leg = ', this%legendre(:10,1)
            this%legendre_coefficients_calculated = .true.
            call log_message('Calculated legendre coefficients', INFO, FILE_DESCRIPTOR(INFO))
            call log_message('Finished spheroidal calculation', INFO, FILE_DESCRIPTOR(INFO))
            !call log_array('s1', this%lnum, this%s1(:,1), this%lnum, FILE_DESCRIPTOR(INFO))
            !call log_array('s1d', this%lnum, this%s1d(:,1), this%lnum, FILE_DESCRIPTOR(INFO))
        endif
        !call log_spheroidal_calculation_detailed(this)
!        write(*,*) 'end d120 = ', this%legendre(:200,120)
    end subroutine calculate

    function get_wronskian(this) result(WR)
        class(SpheroidalCalculation) :: this
        complex(knd) :: WR(this%lnum)

        !write(*, *) this%r1 * this%r2d
        !write(*, *) this%r2 * this%r1d
        !write(*,*) this%r1 * this%r2d -
        WR = this%r1 * this%r2d - this%r2 * this%r1d
    end function get_wronskian

    complex(knd) function theoretical_wronskian_radial(c, ksi, f) result(wr)
        complex(knd) :: c
        real(knd) :: ksi
        integer :: f
        if (f == 1) then
            wr = 1.0q0 / c / (ksi * ksi - 1.0q0)
        else
            wr = 1.0q0 / c / (ksi * ksi + 1.0q0)
        end if
    end function theoretical_wronskian_radial


    complex(knd) function get_c(xv, ab, f, ri)
        real(knd) :: xv, ab
        integer :: f
        complex(knd), optional :: ri

        if (.not. present(ri)) then
            ri = 1q0
        end if

        if (f == 1) then
            get_c = (1q0 / ab)**(1q0 / 3q0)
        else
            get_c = (1q0 / ab)**(2q0 / 3q0)
        endif
        get_c = xv * sqrt(ab**2q0 - 1q0) * get_c * ri

    end function get_c

    real(knd) function get_ksi(ab, f)
        real(knd) :: ab
        integer :: f

        if (f == 1) then
            get_ksi = ab / sqrt(ab * ab - 1q0)
        else
            get_ksi = 1q0 / sqrt(ab * ab - 1q0)
        endif

    end function get_ksi


    subroutine delete_calculation(this)
        type(SpheroidalCalculation), intent(inout) :: this

        call this%delete_calculation_full()

    end subroutine delete_calculation

    subroutine delete_radial_functions(this)
        class(SpheroidalCalculation), intent(inout) :: this

        this%functions_calculated = .false.
        if (allocated(this%r1)) then
            deallocate(this%r1, this%r1d, this%r2, this%r2d, this%r3, this%r3d)
        endif

    end subroutine delete_radial_functions

    subroutine delete_angular_functions(this)
        class(SpheroidalCalculation), intent(inout) :: this

        this%functions_calculated = .false.
        if (allocated(this%arg)) then
            deallocate(this%arg, this%s1, this%s1d)
        endif

    end subroutine delete_angular_functions

    subroutine delete_legendre_coefficients(this)
        class(SpheroidalCalculation), intent(inout) :: this

        this%legendre_coefficients_calculated = .false.
        if (allocated(this%legendre)) then
            deallocate(this%legendre)
        endif

    end subroutine delete_legendre_coefficients

    subroutine delete_calculation_full(this)
        class(SpheroidalCalculation), intent(inout) :: this

        call this%delete_radial_functions()
        call this%delete_angular_functions()
        call this%delete_legendre_coefficients()

        call log_message('Deleted full spheroidal calculation', INFO, fd_spheroidal_calculation)
    end subroutine delete_calculation_full

    !  normalizes the array of legendre coefficients d of the size maxd so that
    !  \sum_{r=0}^\infty d[r]^2 * 2/(2r+2m+1) * (r+2m)!/r! = 1
    subroutine normalize(d, maxd, m)
        integer :: n, r, maxd, m
        complex(knd) d(0:maxd), norm
        real(knd) :: factorial

        factorial = 1q0
        do r = 2, 2 * m
            factorial = factorial * r
        end do

        norm = 0q0
        do r = 0, maxd
            norm = norm + d(r) ** 2 * factorial * 2q0 / (2 * r + 2 * m + 1)
            factorial = factorial * (r + 2 * m + 1) / (r + 1)
        enddo

        !write(*,*) 'norm = ', norm
        norm = sqrt(norm)
        d = d / norm
    end subroutine normalize

    ! LOGGING
    !  logs only parameters for which spheroidal calculation is starting
    subroutine log_spheroidal_calculation(f, m, lnum, c, ksi, narg, arg)
        integer, intent(in) :: m, lnum, narg, f
        real(knd), intent(in) :: ksi, arg(narg)
        complex(knd), intent(in) :: c

        character(len=MESSAGE_LENGTH) :: log
        character(len=10) :: spheroidal_type

        if (f == 1) then
            spheroidal_type = "prolate"
        else
            spheroidal_type = "oblate"
        end if

        log = ''
        write(log, spheroidal_function_log_format) "Calculate ", spheroidal_type, &
                " spheroidal functions with parameters: ", "m = ", m, ", lnum = ", lnum, &
                ", c = ", c, ", ksi = ", ksi, &
                ", with ", narg, " values of eta: ", arg

        call log_message(log, INFO, fd_spheroidal_calculation)
    end subroutine log_spheroidal_calculation

    !  logs a full calculation result
    subroutine log_spheroidal_calculation_detailed(calculation)
        type(SpheroidalCalculation) :: calculation
        character(len=MESSAGE_LENGTH) :: log
        character(len=10) :: spheroidal_type
        integer :: i, length

        length = min(calculation%lnum, 10)

        if (calculation%spheroidal_type == 1) then
            spheroidal_type = "prolate"
        else
            spheroidal_type = "oblate"
        end if

        call log_message("Spheroidal calculation completed for parameters", &
                DETAIL, fd_spheroidal_calculation_detail)

        write(log, *) 'Spheroidal type', spheroidal_type
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'm = ', calculation%m
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'lnum = ', calculation%lnum
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)

        write(log, *) 'c = ', calculation%c
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'ksi = ', calculation%ksi
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'narg = ', calculation%narg
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'arg = ', calculation%arg
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)

        call log_message("Radial functions:", DETAIL, fd_spheroidal_calculation_detail)

        write(log, *) 'r1 = ', calculation%r1(1:length)
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'r1d = ', calculation%r1d(1:length)
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'r2 = ', calculation%r2(1:length)
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'r2d = ', calculation%r2d(1:length)
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'r3 = ', calculation%r3(1:length)
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'r3d = ', calculation%r3d(1:length)
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)

        call log_message("Angular functions:", DETAIL, fd_spheroidal_calculation_detail)

        do i = 1, calculation%narg
            write(log, *) 'argnum = ', i, 'arg = ', calculation%arg(i)
            call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
            write(log, *) 's1 = ', calculation%s1(1:length, i)
            call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
            write(log, *) 's1d = ', calculation%s1d(1:length, i)
            call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        end do

        call log_message("Legendre coefficients:", DETAIL, fd_spheroidal_calculation_detail)
        write(log, *) 'maxd = ', calculation%maxd
        call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        do i = 1, calculation%lnum
            write(log, *) 'n = ', i + calculation%m - 1, 'legendre = ', calculation%legendre(0:length, i)
            call log_message(log, DETAIL, fd_spheroidal_calculation_detail)
        end do

        call log_new_line(DETAIL, fd_spheroidal_calculation_detail)
    end subroutine log_spheroidal_calculation_detailed
end module spheroidal