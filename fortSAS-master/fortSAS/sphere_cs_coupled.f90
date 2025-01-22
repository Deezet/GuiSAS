module sphere_cs_coupled
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

  double precision function p_sphere_cs_coupled(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: RpD, d, SLDcore, SLDshell, SLDmatrix
    double precision :: ff_amp_core, ff_amp_shell
    double precision :: R
    RpD = p(1)
    d = p(2)
    SLDcore = p(3)
    SLDshell = p(4)
    SLDmatrix = p(5)
    R = RpD - d
    if (R < 0) then
      R = 0
    end if
    ff_amp_core = p_sphere(&
      q, (/R, SLDcore, SLDshell/), 3 &
    )

    ff_amp_shell = p_sphere(&
      q, (/RpD, SLDshell, SLDmatrix/) , 3 &
    )
    p_sphere_cs_coupled = ff_amp_shell + ff_amp_core
  end function p_sphere_cs_coupled

  double precision function ff_sphere_cs_coupled(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    ff_sphere_cs_coupled = abs(p_sphere_cs_coupled(q, p, Np))**2
  end function ff_sphere_cs_coupled

  subroutine formfactor(&
    q, RpD, D, SLDcore, SLDshell, SLDmatrix, sigRpD, sigD,&
    Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: RpD, D, SLDcore, SLDshell, SLDmatrix
    double precision, intent(in) :: sigRpD, sigD
    integer, intent(in) :: Nq

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq
    double precision :: RpDmin, RpDmax
    double precision :: Dmin, Dmax

    p = (/RpD, d, SLDcore, SLDshell, SLDmatrix/)
    call get_cutoff_lognormal(RpD, sigRpD, RpDmin, RpDmax)
    call get_cutoff_lognormal(D, sigD, Dmin, Dmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(&
        q(iq), p, Np, &
        1, RpDmin, RpDmax, sigRpD, &
        2, Dmin, Dmax, sigD, &
        ff_sphere_cs_coupled, lognormal, lognormal,&
        ff_intensity(iq) &
      )
    end do
    !$omp end do
    !$omp end parallel
  end subroutine formfactor

  subroutine magnetic_formfactor(&
    q, RpD, D, SLDcore, SLDshell, SLDmatrix, sigRpD, sigD,&
    dmag_deadlayer, mag_SLDcore, mag_SLDshell, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus, Nq, ff_out)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: RpD, D, SLDcore, SLDshell, SLDmatrix
    double precision, intent(in) :: sigRpD, sigD, dmag_deadlayer
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell, mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq
    double precision, dimension(Nq), intent(out) :: ff_out

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p

    integer :: iq
    double precision :: RpDmin, RpDmax
    double precision :: Dmin, Dmax

    p = (/RpD, d, SLDcore, SLDshell, SLDmatrix/)
    call get_cutoff_lognormal(RpD, sigRpD, RpDmin, RpDmax)
    call get_cutoff_lognormal(D, sigD, Dmin, Dmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(&
        q(iq), p, Np, &
        1, RpDmin, RpDmax, sigRpD, &
        2, Dmin, Dmax, sigD, &
        call_magnetic, lognormal, lognormal,&
        ff_out(iq) &
      )
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

        p_mag = (/&
        reducedRadius,&
          p(2),&
          mag_SLDcore,&
          mag_SLDshell,&
          mag_SLDmatrix&
        /)
        call magnetic_scattering(&
          q, xi, sin2alpha, plus_or_minus, p,&
          p_mag, p_sphere_cs_coupled, p_sphere_cs_coupled,&
          Nq, Np, Np, call_magnetic &
        )
      end function call_magnetic
  end subroutine magnetic_formfactor

  subroutine p_crossterm(&
    q, particle_size, d_shell,&
    SLDcore, SLDshell, SLDmatrix, &
    sig_size, dmag_deadlayer, mag_SLDcore, mag_SLDshell, &
    xi, sin2alpha, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: particle_size, d_shell
    double precision, intent(in) :: SLDcore, SLDshell, SLDmatrix
    double precision, intent(in) :: sig_size, dmag_deadlayer
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell
    double precision, intent(in) :: xi, sin2alpha
    integer, intent(in) :: Nq

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: size_min, size_max

    p = (/particle_size, d_shell, SLDcore, SLDshell, SLDmatrix/)
    call get_cutoff_lognormal(particle_size, sig_size, size_min, size_max)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_size_distribution(q(iq), p, Np, &
          1, size_min, size_max, sig_size, &
          call_magnetic, lognormal, ff_intensity(iq))
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

        p_mag = (/&
          reducedRadius, p(2), mag_SLDcore, mag_SLDshell, 0d0&
        /)

        call crossterm(&
          q, xi, sin2alpha, p,&
          p_mag, p_sphere_cs_coupled, p_sphere_cs_coupled,&
          Np, Np, call_magnetic &
        )
      end function call_magnetic
  end subroutine p_crossterm

  subroutine sld(RpD, d, SLDcore, SLDshell, SLDmatrix, x, sld_array)
    double precision, intent(in) :: RpD, d, SLDcore, SLDshell, SLDmatrix

    integer, parameter :: Nx=6
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    double precision :: R
    R = RpD - d
    if (R < 0) then
      R = 0
    end if
    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = RpD
    x(5) = RpD
    x(6) = RpD + 100

    sld_array(1) = SLDcore
    sld_array(2) = SLDcore
    sld_array(3) = SLDshell
    sld_array(4) = SLDshell
    sld_array(5) = SLDmatrix
    sld_array(6) = SLDmatrix
  end subroutine sld
end module sphere_cs_coupled
