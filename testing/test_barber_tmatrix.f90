program test_barber_tmatrix
    use spheroidal
    use wavelength_point
    use logging
    use communication
    use spheroidal_scatterer
    use spheroidal_scattering
    use scattering_functions
    use tuning
    use constants
    implicit none

    type(SpheroidalCalculation) :: calculation
    type(WavelengthPoint) :: point
    type(ScatteringQuery), target :: query
    type(ScatteringResult) :: response
    type(SpheroidalScatterer) :: scatterer
    type(ScatteringFactors) :: result, &
            spheroidal_te_uv, spheroidal_te_pq, spheroidal_tm_uv, spheroidal_tm_pq, &
            spherical_te_uv, spherical_te_pq, spherical_tm_uv, spherical_tm_pq
    type(SpheroidalShape) :: shape
    integer :: no_abs_ms, cons_ms, f, i, j, k, matrix_size, nol, minm, maxm, m
    real(knd) :: a, b, lambda, alpha, xa
    real(knd), allocatable :: xv(:), rv(:), ab(:)
    complex(knd), allocatable :: ri(:)
    character(1024) :: filename
    real(knd), allocatable, dimension(:) :: ab_s, rv_s, xv_s
    complex(knd), allocatable, dimension(:) :: ri_s
    character(256) :: name

    call log_init()
    filename = '../input.txt'
    if (command_argument_count() > 0) then
        call get_command_argument(1, filename)
    end if
    call read_input('../input.txt', f, nol, rv, xv, ab, alpha, lambda, ri, matrix_size, minm, maxm)

    response = GetSphericalTmatrix(f, nol, xv, ab, alpha, &
            lambda, &
            ri, matrix_size, minm, maxm)
    call log_message('TE mode', WARNING, FILE_DESCRIPTOR(WARNING))
    do m = minm, response%actual_maxm
        write(name, *) 'TUVBAR m = ', m
        call log_matrix(FILE_DESCRIPTOR(WARNING), name, &
                response%result_spherical_te%uv_tmatrix(:2*matrix_size,:2*matrix_size,m), .false., 2*matrix_size)
    end do

    call log_message('TM mode', WARNING, FILE_DESCRIPTOR(WARNING))
    do m = minm, response%actual_maxm
        write(name, *) 'TUVBAR m = ', m
        call log_matrix(FILE_DESCRIPTOR(WARNING), name, &
                response%result_spherical_tm%uv_tmatrix(:2*matrix_size,:2*matrix_size,m), .false., 2*matrix_size)
    end do

    ! nonconfocal
    !    response = GetSpheroidalQFactors(1, 2, (/0.727415757314481q0, 0.5099122256638764q0/), (/ 2.0q0, 3.0q0 /), 0.523598q0, &
    !            2.0q0*PI, &
    !            (/cmplx(1q0, 0q0, knd), cmplx(1.2q0, 0.0q0, knd), cmplx(1.4q0, -0.0q0, knd) /), matrix_size, 0)

    spheroidal_te_uv = response%result_te%uv_factors
    spheroidal_te_pq = response%result_te%pq_factors
    spheroidal_tm_uv = response%result_tm%uv_factors
    spheroidal_tm_pq = response%result_tm%pq_factors

    write(*,*) 'Q_norm_spheroidal_te_pq'
    call spheroidal_te_pq%log(WARNING, FILE_DESCRIPTOR(WARNING))
    write(*,*) 'Q_norm_spheroidal_te_uv'
    call spheroidal_te_uv%log(WARNING, FILE_DESCRIPTOR(WARNING) )
    write(*,*) 'Q_norm_spheroidal_tm_pq'
    call spheroidal_tm_pq%log(WARNING, FILE_DESCRIPTOR(WARNING))
    write(*,*) 'Q_norm_spheroidal_tm_uv'
    call spheroidal_tm_uv%log(WARNING, FILE_DESCRIPTOR(WARNING))

    deallocate(xv, rv, ab, ri)
    call log_close()

end program test_barber_tmatrix