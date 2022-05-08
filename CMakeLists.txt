cmake_minimum_required(VERSION 3.16)
project(ScatteringOnSpheroids Fortran)


enable_language(Fortran)

SET(SOURCES_PARAMETERS
        src/parameters/constants.f90
        src/parameters/regime.f90)
SET(SOURCES_UTILS
        src/utils/logging.f90
        src/utils/communication.f90
        src/utils/geometry.f90
        src/utils/tuning.f90
        src/utils/mfunctions.f90
        )
SET(SOURCES_MATH
        src/mathematics/angle.f90
        src/mathematics/matrix.f90)
SET(SOURCES_SPHEROIDAL_FUNCTIONS
        src/mathematics/functions/param.f90
        src/mathematics/functions/complex_prolate_swf.f90
        src/mathematics/functions/complex_oblate_swf.f90
        src/mathematics/functions/spheroidal.f90
        src/mathematics/pair_integrals/spheroidal_integrals.f90)
SET(SOURCES_SPHERICAL
        src/mathematics/functions/legendre_functions.f90
        src/mathematics/pair_integrals/spherical_integrals.f90
        src/scattering/spherical/spherical_scattering.f90)
SET(SOURCES_SCATTERING
        src/scattering/wavelength_point.f90
        src/scattering/scatterer.f90
        src/scattering/spheroidal/spheroidal_initial.f90
        src/scattering/spheroidal/spheroidal_tmatrix.f90
        src/scattering/spheroidal/spheroidal_scattering.f90
        src/scattering/spheroidal/spheroidal_indicatrix.f90)
SET(SOURCES_ALL
        ${SOURCES_PARAMETERS} ${SOURCES_UTILS} ${SOURCES_MATH} ${SOURCES_SPHEROIDAL_FUNCTIONS}
        ${SOURCES_SPHERICAL} ${SOURCES_SCATTERING} src/scattering_functions.f90)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")
add_executable(ScatteringOnSpheroids main.f90 ${SOURCES_ALL})
add_executable(test_c_of_n testing/test_c_of_n.f90 ${SOURCES_ALL})
add_executable(test_barber_tmatrix testing/test_barber_tmatrix.f90 ${SOURCES_ALL})
