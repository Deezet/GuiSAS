module sphere_cs
use math
implicit none
contains
  double precision function p_sphere(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: qR, R, SLDsphere, SLDmatrix
    R = p(1)
    SLDsphere = p(2)
    SLDmatrix = p(3)

    qR = q*R
    if (qR /= 0) then
      p_sphere = 4d0 * pi * R**3 * (SLDsphere - SLDmatrix) *&
             (sin(qR) - qR*cos(qR))/qR**3
    else
      p_sphere = 4d0/3d0 * pi * R**3 * (SLDsphere - SLDmatrix)
    end if
  end function p_sphere

  double precision function p_sphere_cs(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: R, d, SLDcore, SLDshell, SLDmatrix
    double precision :: ff_amp_core, ff_amp_shell
    R = p(1)
    d = p(2)
    SLDcore = p(3)
    SLDshell = p(4)
    SLDmatrix = p(5)
    ff_amp_shell = p_sphere(q, (/R+d, SLDshell, SLDmatrix/), 3)
    ff_amp_core = p_sphere(q, (/R, SLDcore, SLDshell/), 3)
    p_sphere_cs = ff_amp_shell + ff_amp_core
  end function p_sphere_cs

  double precision function ff_sphere_cs(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    ff_sphere_cs = abs(p_sphere_cs(q, p, Np))**2
  end function ff_sphere_cs

  subroutine formfactor(&
    q, R, D, SLDcore, SLDshell, SLDmatrix, sigR, sigD, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, D, SLDcore, SLDshell, SLDmatrix
    double precision, intent(in) :: sigR, sigD
    integer, intent(in) :: Nq

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq
    double precision :: Rmin, Rmax
    double precision :: Dmin, Dmax

    p = (/R, d, SLDcore, SLDshell, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    call get_cutoff_lognormal(D, sigD, Dmin, Dmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(&
        q(iq), p, Np, &
        1, Rmin, Rmax, sigR, &
        2, Dmin, Dmax, sigD, &
        ff_sphere_cs, lognormal, lognormal,&
        ff_intensity(iq)&
      )
    end do
    !$omp end do
    !$omp end parallel
  end subroutine formfactor

  subroutine magnetic_formfactor(&
    q, R, D, SLDcore, SLDshell, SLDmatrix, sigR, sigD,&
    dmag_deadlayer, mag_SLDcore, mag_SLDshell, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, D, SLDcore, SLDshell, SLDmatrix
    double precision, intent(in) :: sigR, sigD, dmag_deadlayer
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell
    double precision, intent(in) :: mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq
    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p

    integer :: iq
    double precision :: Rmin, Rmax
    double precision :: Dmin, Dmax

    p = (/R, d, SLDcore, SLDshell, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    call get_cutoff_lognormal(D, sigD, Dmin, Dmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(q(iq), p, Np, &
            1, Rmin, Rmax, sigR, &
            2, Dmin, Dmax, sigD, &
            call_magnetic, lognormal, lognormal,&
            ff_intensity(iq))
    end do
    !$omp end do
    !$omp end parallel

    contains
      double precision function call_magnetic(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np

        double precision, dimension(Np) :: p_mag
        double precision :: reducedRadius
        reducedRadius = p(1)-dmag_deadlayer
        if (reducedRadius < 0) then
          reducedRadius = 0
        end if

        p_mag = (/reducedRadius, p(2), mag_SLDcore, mag_SLDshell, mag_SLDmatrix/)
        call magnetic_scattering(&
          q, xi, sin2alpha, plus_or_minus, p, p_mag,&
          p_sphere_cs, p_sphere_cs, Nq, Np, Np, call_magnetic)
      end function call_magnetic
  end subroutine magnetic_formfactor

  subroutine sld(R, d, SLDcore, SLDshell, SLDmatrix, x, sld_array)
    double precision, intent(in) :: R, d, SLDcore, SLDshell, SLDmatrix

    integer, parameter :: Nx=6
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = R + d
    x(5) = R + d
    x(6) = R + d + 100

    sld_array(1) = SLDcore
    sld_array(2) = SLDcore
    sld_array(3) = SLDshell
    sld_array(4) = SLDshell
    sld_array(5) = SLDmatrix
    sld_array(6) = SLDmatrix
  end subroutine sld
end module sphere_cs
