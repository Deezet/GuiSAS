module sphere
use math
implicit none
contains
  double precision function p_sphere(q, p, Np)
    ! Sphere Formfacotr Amplitude
    ! Parameters: R, SLDsphere, SLDmatrix
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: qR, R, SLDsphere, SLDmatrix
    R = p(1)
    SLDsphere = p(2)
    SLDmatrix = p(3)

    qR = q*R
    if (qR /= 0) then
      p_sphere = 4d0 * pi * R**3 * (SLDsphere - SLDmatrix) * &
             (sin(qR) - qR*cos(qR))/qR**3
    else
      p_sphere = 4d0/3d0 * pi * R**3 * (SLDsphere - SLDmatrix)
    end if
  end function p_sphere

  double precision function ff_sphere(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
    ff_sphere = abs(p_sphere(q, p, Np))**2
  end function ff_sphere

  subroutine formfactor(&
    q, R, SLDsphere, SLDmatrix, sigR, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, SLDsphere, SLDmatrix
    double precision, intent(in) :: sigR
    integer, intent(in) :: Nq

    integer, parameter :: Np=3
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq
    double precision :: Rmin, Rmax

    p = (/R, SLDsphere, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)

    do iq=1, Nq
      call integrate_size_distribution(q(iq), p, Np, &
              1, Rmin, Rmax, sigR, &
              ff_sphere, lognormal, ff_intensity(iq))
    end do
  end subroutine formfactor

  subroutine magnetic_formfactor(&
    q, R, SLDsphere, SLDmatrix, sigR, mag_SLDsphere, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus, Nq, ff_out)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, SLDsphere, SLDmatrix
    double precision, intent(in) :: sigR
    double precision, intent(in) :: mag_SLDsphere, mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq

    integer, parameter :: Np=3
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_out

    integer :: iq
    double precision :: Rmin, Rmax

    p = (/R, SLDsphere, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)

    do iq=1, Nq
      call integrate_size_distribution(q(iq), p, Np, &
            1, Rmin, Rmax, sigR, &
            call_magnetic, lognormal,&
            ff_out(iq))
    end do
    contains
      double precision function call_magnetic(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np

        double precision, dimension(Np) :: p_mag
        p_mag = (/p(1), mag_SLDsphere, mag_SLDmatrix/)
        call magnetic_scattering(&
          q, xi, sin2alpha, plus_or_minus, p, p_mag,&
          p_sphere, p_sphere, Nq, Np, Np, call_magnetic)
      end function call_magnetic
  end subroutine magnetic_formfactor

  subroutine sld(R, SLDsphere, SLDmatrix, x, sld_array)
    double precision, intent(in) :: R, SLDsphere, SLDmatrix

    integer, parameter :: Nx=4
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = R + 100

    sld_array(1) = SLDsphere
    sld_array(2) = SLDsphere
    sld_array(3) = SLDmatrix
    sld_array(4) = SLDmatrix
  end subroutine sld
end module sphere