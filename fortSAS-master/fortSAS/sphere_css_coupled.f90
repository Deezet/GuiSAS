module sphere_css_coupled
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
      p_sphere = 4d0 * pi * R**3 * (SLDsphere - SLDmatrix) * &
             (sin(qR) - qR*cos(qR))/qR**3
    else
      p_sphere = 4d0/3d0 * pi * R**3 * (SLDsphere - SLDmatrix)
    end if
  end function p_sphere

  double precision function p_sphere_css_coupled(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: particle_size, dshell, R, dsurfactant
    double precision :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
    double precision :: p_amp_core, p_amp_shell, p_amp_surfactant
    particle_size = p(1)
    dshell = p(2)
    dsurfactant = p(3)
    SLDcore = p(4)
    SLDshell = p(5)
    SLDsurfactant = p(6)
    SLDmatrix = p(7)

    R = particle_size - dshell
    if (R < 0) then
      R = 0
    end if
    p_amp_core = p_sphere( &
      q, (/R, SLDcore, SLDshell/), 3 &
    )

    p_amp_shell = p_sphere( &
      q, (/particle_size, SLDshell, SLDsurfactant/), 3 &
    )

    p_amp_surfactant = p_sphere( &
      q, (/particle_size+dsurfactant, SLDsurfactant, SLDmatrix/), 3 &
    )

    p_sphere_css_coupled = p_amp_surfactant + p_amp_shell + p_amp_core
  end function p_sphere_css_coupled

  double precision function ff_sphere_css_coupled(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
    ff_sphere_css_coupled = &
            abs(p_sphere_css_coupled(q, p, Np))**2
  end function ff_sphere_css_coupled

  subroutine formfactor(&
    q, particle_size, dshell, dsurfactant, &
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sigRpD, sigDshell, Nq, ff_intensity)

    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: particle_size, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant
    double precision, intent(in) :: SLDmatrix
    double precision, intent(in) :: sigRpD, sigDshell
    integer, intent(in) :: Nq

    integer, parameter :: Np=7
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_intensity


    integer :: iq
    double precision :: RpDmin, RpDmax
    double precision :: Dshell_min, Dshell_max

    p = (/particle_size, dshell, dsurfactant,&
        SLDcore, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(particle_size, sigRpD, RpDmin, RpDmax)
    call get_cutoff_lognormal(dshell, sigDshell,  Dshell_min, Dshell_max)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(&
        q(iq), p, Np, &
        1, RpDmin, RpDmax, sigRpD, &
        2, dshell_min, dshell_max, sigDshell, &
        ff_sphere_css_coupled, lognormal, lognormal,&
        ff_intensity(iq) &
      )
    end do
    !$omp end do
    !$omp end parallel
  end subroutine formfactor

  subroutine magnetic_formfactor( &
    q, particle_size, dshell, dsurfactant, &
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sigRpD, sigDshell, &
    mag_SLDcore, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus, Nq, ff_intensity)

    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: particle_size, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
    double precision, intent(in) :: sigRpD, sigDshell
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell, mag_SLDsurfactant
    double precision, intent(in) :: mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq

    integer, parameter :: Np=7
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: RpDmin, RpDmax
    double precision :: Dshell_min, Dshell_max

    p = (/particle_size, dshell, dsurfactant,&
        SLDcore, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(particle_size, sigRpD, RpDmin, RpDmax)
    call get_cutoff_lognormal(dshell, sigDshell, Dshell_min, Dshell_max)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(&
        q(iq), p, Np, &
        1, RpDmin, RpDmax, sigRpD, &
        2, Dshell_min, Dshell_max, sigDshell, &
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

        p_mag = (/p(1), p(2), p(3),&
              mag_SLDcore, mag_SLDshell, mag_SLDsurfactant,&
              mag_SLDmatrix/)
        call magnetic_scattering(q, xi, sin2alpha, plus_or_minus, p,&
          p_mag, p_sphere_css_coupled, p_sphere_css_coupled,&
          Nq, Np, Np, call_magnetic)
      end function call_magnetic
  end subroutine magnetic_formfactor

  subroutine p_crossterm(&
    q, particle_size, d_shell, dsurfactant, &
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sig_size, &
    mag_SLDcore, mag_SLDshell, mag_SLDsurfactant,&
    xi, sin2alpha, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: particle_size, d_shell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
    double precision, intent(in) :: sig_size
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell, mag_SLDsurfactant
    double precision, intent(in) :: xi, sin2alpha
    integer, intent(in) :: Nq

    integer, parameter :: Np=7
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: size_min, size_max

    p = (/particle_size, d_shell, dsurfactant,&
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix/)
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

        p_mag = (/&
        p(1), p(2), p(3), mag_SLDcore, mag_SLDshell, mag_SLDsurfactant, 0d0&
        /)

        call crossterm(&
          q, xi, sin2alpha, p,&
          p_mag, p_sphere_css_coupled, p_sphere_css_coupled,&
          Np, Np, call_magnetic &
        )
      end function call_magnetic
  end subroutine p_crossterm

  subroutine sld(&
    particle_size, dshell, dsurfactant, &
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, x, sld_array)

    double precision, intent(in) :: particle_size, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix

    integer, parameter :: Nx=8
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    double precision :: R

    R = particle_size - dshell
    if (R < 0) then
      R = 0
    end if
    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = particle_size
    x(5) = particle_size
    x(6) = particle_size + dsurfactant
    x(7) = particle_size + dsurfactant
    x(8) = particle_size + dsurfactant + 100

    sld_array(1) = SLDcore
    sld_array(2) = SLDcore
    sld_array(3) = SLDshell
    sld_array(4) = SLDshell
    sld_array(5) = SLDsurfactant
    sld_array(6) = SLDsurfactant
    sld_array(7) = SLDmatrix
    sld_array(8) = SLDmatrix
  end subroutine sld
end module sphere_css_coupled
