! Created by drakosha on 07.10.2021.

program test_c_of_n
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
            spherical_te_uv, spherical_te_pq, spherical_tm_uv, spherical_tm_pq, &
            spherical_orth_te_uv, spherical_orth_te_pq, spherical_orth_tm_uv, spherical_orth_tm_pq
    type(SpheroidalShape) :: shape
    integer :: no_abs_ms, cons_ms, f, i, j, k, matrix_size, nol, matrix_size0, minm, maxm
    real(knd), allocatable :: ab(:), xv(:), rv(:)
    complex(knd), allocatable :: ri(:)
    real(knd) :: a(2), b(2), lambda, alpha, xa, d(2), dv
    complex(knd) :: ca, cb, cd, cx, cab, ct, cq, alp, beta
    character(1024) :: filename
    real(knd), allocatable, dimension(:) :: ab_s, rv_s, xv_s
    complex(knd), allocatable, dimension(:) :: ri_s

    call log_init()
    filename = '../input.txt'
    if (command_argument_count() > 0) then
        call get_command_argument(1, filename)
    end if
    call read_input('../input.txt', f, nol, rv, xv, ab, alpha, lambda, ri, matrix_size0, minm, maxm)
    98 format('#',1A5,' ',20A25)
    if (f == 1) then
        b = (rv**3 / ab) ** (1q0 / 3q0)
    else
        b = (rv**3 / ab**2) ** (1q0 / 3q0)
    end if
    a = b * ab
    write(*,*) 'V_c / V = ', a(2)*b(2)**2 / (a(1)*b(1)**2)
    99 format('#',1A5,' ',8A25)
    write(*,99) 'N', &
            'Q_spheroidal_te_uv_ext', 'Q_spheroidal_te_pq_ext', &
            'Q_spheroidal_te_uv_sca', 'Q_spheroidal_te_pq_sca',&
            'Q_spheroidal_tm_uv_ext', 'Q_spheroidal_tm_pq_ext', &
            'Q_spheroidal_tm_uv_sca', 'Q_spheroidal_tm_pq_sca'
    do matrix_size = matrix_size0, 60, 8
        response = GetSpheroidalQFactors(f, nol, xv, ab, alpha, &
                lambda, &
                ri, matrix_size, minm,maxm)

        spheroidal_te_uv = response%result_te%uv_factors
        spheroidal_te_pq = response%result_te%pq_factors
        spheroidal_tm_uv = response%result_tm%uv_factors
        spheroidal_tm_pq = response%result_tm%pq_factors

        100 format(' ',1I5,' ',8E25.15)
        write(*,100) matrix_size, &
                spheroidal_te_uv%ext, spheroidal_te_pq%ext, &
                spheroidal_te_uv%sca, spheroidal_te_pq%sca,&
                spheroidal_tm_uv%ext, spheroidal_tm_pq%ext, &
                spheroidal_tm_uv%sca, spheroidal_tm_pq%sca
    end do
    deallocate(ab, xv, rv, ri)
    call log_close()

end program test_c_of_n