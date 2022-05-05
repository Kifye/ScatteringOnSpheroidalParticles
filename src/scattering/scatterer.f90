! Created by drakosha on 15.07.2021.

module spheroidal_scatterer
    use angle
    use regime
    use constants
    implicit none
    private

    ! represents geometric properties of the scatterer
    type, public :: SpheroidalScatterer

        integer :: spheroidal_type

        !  parameters used during calculation obtained from shape

        integer :: number_of_layers
        !  radial coordinates of the scatterer's layers
        real(knd), allocatable, dimension(:) :: ksi   !  1..number_of_layers
        real(knd), allocatable, dimension(:) :: d   !  1..number_of_layers
        !  angle between the spheroidal symmetry axis and the direction of incident light propagation
        type(AngleType) :: alpha
        !  precalculated factor from extinction and scattering coefficient
        complex(knd), allocatable, dimension(:) :: c0   !  1..number_of_layers
        !  argument of the spheroidal function corresponding to the external field (layer 0)
        real(knd), allocatable, dimension(:) :: common_factor   !  1..number_of_layers

    contains

        procedure, public :: set

        final :: delete_scatterer
    end type SpheroidalScatterer

    type, public :: SpheroidalShape
        ! shape of the scatterrer

        integer :: spheroidal_type
        !  radius of the equivalent sphere
        real(knd) :: rv
        !  major semiaxis
        real(knd) :: a
        !  minor semiaxes
        real(knd) :: b
        !  a / b
        real(knd) :: ab

        type(AngleType) :: alpha
    contains
        procedure, public :: set => set_spheroidal_shape
    end type SpheroidalShape

contains
    !  for the given rv and ab returns a - the major semiaxis
    real(knd) function get_a(f, rv, ab) result(a)
        integer :: f
        real(knd) :: rv, ab

        if (f == 1) then
            a = rv * (ab ** (2q0 / 3q0))
        else
            a = rv * (ab ** (1q0 / 3q0))
        endif
    end function get_a

    subroutine set(this, f, xv, ab, alpha, number_of_layers)

        class(SpheroidalScatterer) :: this
        real(knd) :: xv(number_of_layers), ab(number_of_layers), alpha, a(number_of_layers), b(number_of_layers)
        integer :: f, number_of_layers

        this%spheroidal_type = f

        if (allocated(this%ksi) .and. number_of_layers /= this%number_of_layers) then
            deallocate(this%ksi, this%c0, this%common_factor, this%d)
        end if
        if (.not. allocated(this%ksi)) then
            allocate(this%ksi(1:number_of_layers), this%c0(1:number_of_layers), this%common_factor(1:number_of_layers),&
            this%d(1:number_of_layers))
        end if
        this%number_of_layers = number_of_layers

        call this%alpha%set(alpha)


        if (f == 1) then
            this%ksi = ab / sqrt(ab * ab - 1q0)
            this%c0 = (1q0 / ab)**(1q0 / 3q0)
        else
            this%ksi = 1q0 / sqrt(ab * ab - 1q0)
            this%c0 = (1q0 / ab)**(2q0 / 3q0)
        endif
        this%c0 = xv * sqrt(ab**2q0 - 1q0) * this%c0
        b = xv / ab**(1q0/3q0)
        a = ab * b
        this%d = (a * a - b * b) ** 0.5q0


        this%common_factor = 1q0 / (abs(this%c0)**2q0 * &
                sqrt((this%ksi**2 - this%spheroidal_type) * (this%ksi**2 - this%spheroidal_type * this%alpha%angle_cos**2)))

    end subroutine set

    subroutine set_spheroidal_shape(this, f, rv, ab, alpha)

        class(SpheroidalShape) :: this
        real(knd) :: rv, ab, alpha
        integer :: f

        this%spheroidal_type = f
        this%rv = rv
        this%a = get_a(f, rv, ab)
        this%b = this%a / ab
        this%ab = ab
        call this%alpha%set(alpha)

    end subroutine set_spheroidal_shape

    subroutine delete_scatterer(this)

        type(SpheroidalScatterer), intent(inout) :: this

        if (allocated(this%ksi)) then
            deallocate(this%ksi, this%c0, this%common_factor, this%d)
        end if

    end subroutine delete_scatterer
end module spheroidal_scatterer