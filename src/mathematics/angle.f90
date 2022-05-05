module angle
    use regime
    use constants
    implicit none
    private
    !  structure intended to optimize cos and sin calculation
    type, public :: AngleType
        real(knd) :: value
        real(knd) :: angle_cos
        real(knd) :: angle_sin
    contains
        procedure, public :: set => set_angle
        procedure, public :: in_degrees
    end type AngleType

    !  Represents spherical coordinates in the system where:
    !  1. OZ axis = symmtery axis of the spheroidal scatterer
    !  2. incident light is in XOZ plane
    !  theta is the angle between OZ and a given direction
    !  phi is the angle between XOZ plane and a given direction
    !  positive towards the OY axis
    type, public :: SphericalDirection
        type(AngleType) :: theta
        type(AngleType) :: phi
    contains
        procedure, public :: set => set_spherical_direction
    end type SphericalDirection

contains

    subroutine set_angle(this, value)
        class(AngleType), intent(out) :: this
        real(knd), intent(in) :: value

        this%value = value
        this%angle_cos = cos(value)
        this%angle_sin = sin(value)
    end subroutine set_angle

    real(knd) function in_degrees(this)
        class(AngleType), intent(in) :: this

        in_degrees = this%value * 180.0_knd / PI
    end function in_degrees

    subroutine set_spherical_direction(direction, theta, phi)
        class(SphericalDirection), intent(out) :: direction
        real(knd), intent(in) :: theta, phi

        call direction%theta%set(theta)
        call direction%phi%set(phi)
    end subroutine set_spherical_direction
end module angle