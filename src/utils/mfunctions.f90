! Created by drakosha on 20.10.2021.

module mfunctions
contains
    function mz_p(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), mz_p

        mz_p = cos(m * phi) * (cos(theta) * sin(theta) * pdr(n - m) * bes(n) / ro - sin(theta) * pr(n - m) * besd(n))

    end function mz_p

    function mz_t(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), mz_t

        mz_t = -m * sin(m * phi) * cos(theta) / sin(theta) * pr(n - m) * bes(n) / ro

    end function mz_t

    function mz_r(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), mz_r

        mz_r = -m * sin(m * phi) * pr(n - m) * bes(n) / ro

    end function mz_r

    function mr_p(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), mr_p

        mr_p = cos(m * phi) * pdr(n - m) * bes(n) * sin(theta)

    end function mr_p

    function mr_t(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), mr_t

        mr_t = -m * sin(m * phi) * pr(n - m) / sin(theta) * bes(n)

    end function mr_t

    function mr_r(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), mr_r

        mr_r = 0

    end function mr_r

    function nr_r(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), nr_r

        nr_r = sin(m * phi) * n * (n + 1) * pr(n - m) * bes(n) / ro

    end function nr_r

    function nr_t(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), nr_t

        nr_t = -sin(m * phi) * sin(theta) * pdr(n - m) * (besd(n) + bes(n) / ro)

    end function nr_t

    function nr_p(m, n, l, theta, phi, ro, size, pr, pdr, bes, besd)
        use regime
        implicit none
        integer :: m, n, l, size
        real(knd) :: theta, phi
        complex(knd) :: ro, pr(0:size), pdr(0:size), bes(0:size), besd(0:size), nr_p

        nr_p = m * cos(m * phi) * pr(n - m) / sin(theta) * (besd(n) + bes(n) / ro)

    end function nr_p

end module mfunctions