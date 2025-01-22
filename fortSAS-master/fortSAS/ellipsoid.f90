module ellipsoid
use math
implicit none
contains
  double precision function p_ellipsoid(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: L, R, alpha, SLDellipsoid, SLDmatrix
    double precision :: M
    L = p(1)
    R = p(2)
    SLDellipsoid = p(3)
    SLDmatrix = p(4)
    alpha = p(5)

    M = (L**2*(cos(alpha*deg))**2 + R**2*(sin(alpha*deg))**2)**0.5

    p_ellipsoid = 4d0*pi*R**2*L*(SLDellipsoid - SLDmatrix)*(sin(q*M)-q*M*cos(q*M))/(q*M)**3
  end function p_ellipsoid

  double precision function ff_ellipsoid(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    ff_ellipsoid = abs(p_ellipsoid(q, p, Np))**2
  end function ff_ellipsoid
  
  double precision function ff_ellipsoid_orientation_average(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        
        call integrate_size_distribution(q, p, Np, &
                            5, 0d0, 90d0, 1d0, &
                            ff_ellipsoid, sinx, &
                            ff_ellipsoid_orientation_average)

  end function ff_ellipsoid_orientation_average

  subroutine formfactor(q, L, R, SLDellipsoid, SLDmatrix, &
                    sigL, sigR,&
                    Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: L, R, SLDellipsoid, SLDmatrix
    double precision, intent(in) :: sigL, sigR
    integer, intent(in) :: Nq

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    double precision, dimension(Nq) :: ff_amplitude
    integer :: iq
    double precision :: qval
    double precision :: Lmin, Lmax, Rmin, Rmax
    double precision :: ff_norm, ff_normR
    !double precision :: alpha_probability
    p = (/L, R, SLDellipsoid, SLDmatrix, 0d0/)
    call get_cutoff_lognormal(L, sigL, Lmin, Lmax)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    !call get_cutoff_gaussian(alpha, sigalpha, alphamin, alphamax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
        call integrate_two_size_distributions(q(iq), p, Np, &
                        1, Lmin, Lmax, sigL, &
                        2, Rmin, Rmax, sigR, &
                        ff_ellipsoid_orientation_average, lognormal, lognormal, &
                        ff_intensity(iq))
    end do
    !$omp end do
    !$omp end parallel

    !call integrate_size_distribution(1d0, p, Np, &
    !                     5, 0d0, 90d0, 1d0, &
    !                     func_one, sinx,&
    !                     ff_norm)
    !call integrate_two_size_distributions(1d0, p, Np, &
    !                     1, Lmin, Lmax, sigL, &
    !                     2, Rmin, Rmax, sigR, &
    !                     func_one, lognormal,lognormal,&
    !                     ff_normR)

    !call integrate_prob_distribution(alpha, sigalpha, alphamin, alphamax,&
    !                sindeg_gaussian, alpha_probability)
    ff_intensity = ff_intensity!/ff_norm/ff_normR
  end subroutine formfactor

  subroutine magnetic_formfactor(q, L, R, alpha, SLDellipsoid, SLDmatrix, &
                    sigL, sigR, sigalpha, mag_SLDellipsoid, mag_SLDmatrix,&
                    xi, sin2alpha, plus_or_minus, Nq, ff_out)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: L, R, alpha, SLDellipsoid, SLDmatrix
    double precision, intent(in) :: sigL, sigR, sigalpha
    double precision, intent(in) :: mag_SLDellipsoid, mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq
    double precision, dimension(Nq), intent(out) :: ff_out

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p

    integer :: iq
    double precision :: qval
    double precision :: Lmin, Lmax
    double precision :: Rmin, Rmax
    double precision :: alphamin, alphamax
    double precision :: alpha_probability

    call get_cutoff_lognormal(L, sigL, Lmin, Lmax)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    call get_cutoff_gaussian(alpha, sigalpha, alphamin, alphamax)

    p = (/L, R, alpha, SLDellipsoid, SLDmatrix/)

    !$omp parallel
    !$omp do
    do iq=1, Nq
        call integrate_three_size_distributions(q(iq), p, Np, &
                        1, Lmin, Lmax, sigL, &
                        2, Rmin, Rmax, sigR, &
                        3, alphamin, alphamax, sigalpha, &
                        ff_ellipsoid, lognormal, lognormal, sindeg_gaussian,&
                        ff_out(iq))
    end do
    !$omp end do
    !$omp end parallel

    call integrate_prob_distribution(alpha, sigalpha, alphamin, alphamax,&
                    sindeg_gaussian, alpha_probability)
    ff_out = ff_out / alpha_probability

    contains
        double precision function call_magnetic(q, p, Np)
            double precision, intent(in) :: q
            double precision, dimension(Np), intent(in) :: p
            integer, intent(in) :: Np

            double precision, dimension(Np) :: p_mag
            p_mag = (/p(1), p(2), p(3), mag_SLDellipsoid, mag_SLDmatrix/)
            call magnetic_scattering(q, xi, sin2alpha, plus_or_minus, p,&
                p_mag, p_ellipsoid, p_ellipsoid,&
                Nq, Np, Np, call_magnetic)
        end function call_magnetic
    end subroutine magnetic_formfactor

  subroutine sld(L, SLDellipsoid, SLDmatrix, x, sld_array)
    double precision, intent(in) :: L, SLDellipsoid, SLDmatrix

    integer, parameter :: Nx=4
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

      x(1) = 0
      x(2) = L
      x(3) = L
      x(4) = L + 100

      sld_array(1) = SLDellipsoid
      sld_array(2) = SLDellipsoid
      sld_array(3) = SLDmatrix
      sld_array(4) = SLDmatrix
  end subroutine sld
end module ellipsoid
