module sphere_css_dead
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

  double precision function p_multi(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: R, d_shell, d_dead, d_surfactant
    double precision :: SLDcore, SLDshell, SLDdead, SLDsurfactant, SLDmatrix
    double precision :: p_amp_core, p_amp_shell, p_amp_dead, p_amp_surfactant

    R = p(1)
    d_shell = p(2)
    d_dead = p(3)
    d_surfactant = p(4)
    SLDcore = p(5)
    SLDshell = p(6)
    SLDdead = p(7)
    SLDsurfactant = p(8)
    SLDmatrix = p(9)

    if(d_dead > d_shell) then
      d_dead = d_shell
    end if

    p_amp_core = p_sphere(&
      q, (/R, SLDcore, SLDshell/), 3)
    p_amp_shell = p_sphere(&
      q, (/R + d_shell - d_dead, SLDshell, SLDdead/), 3)
    p_amp_dead = p_sphere(&
      q, (/R + d_shell, SLDdead, SLDsurfactant/), 3)
    p_amp_surfactant = p_sphere(&
      q, (/R + d_shell + d_surfactant, SLDsurfactant, SLDmatrix/), 3)
    p_multi = p_amp_core + p_amp_shell + p_amp_dead + p_amp_surfactant
  end function p_multi

  double precision function ff_multi_coupled(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
    ff_multi_coupled = abs(p_multi(q, p, Np))**2
  end function ff_multi_coupled

  subroutine formfactor(&
    q, R, d_shell, d_dead, d_surfactant, &
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sigR, Nq, ff_intensity)

    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, d_shell, d_dead, d_surfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
    double precision, intent(in) :: sigR
    integer, intent(in) :: Nq

    integer, parameter :: Np=9
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: Rmin, Rmax

    p = (/R, d_shell, d_dead, d_surfactant,&
        SLDcore, SLDshell, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_size_distribution(&
        q(iq), p, Np, &
        1, Rmin, Rmax, sigR, &
        ff_multi_coupled, lognormal,&
        ff_intensity(iq))
    end do
    !$omp end do
    !$omp end parallel
  end subroutine formfactor

  subroutine magnetic_formfactor(&
    q, R, d_shell, d_dead, d_surfactant,&
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sigR, mag_SLDcore, mag_SLDshell, &
    xi, sin2alpha, plus_or_minus, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, d_shell, d_dead, d_surfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant
    double precision, intent(in) :: SLDmatrix, sigR
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq

    integer, parameter :: Np=9
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: Rmin, Rmax

    p = (/R, d_shell, d_dead, d_surfactant,&
        SLDcore, SLDshell, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_size_distribution(&
        q(iq), p, Np, &
        1, Rmin, Rmax, sigR, &
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

        p_mag = (/&
          p(1), p(2), p(3), p(4),&
          mag_SLDcore, mag_SLDshell, 0d0, 0d0, 0d0&
        /)

        call magnetic_scattering(&
          q, xi, sin2alpha, plus_or_minus, p,&
          p_mag, p_multi, p_multi,&
          Nq, Np, Np, call_magnetic &
        )
      end function call_magnetic
  end subroutine magnetic_formfactor

  subroutine p_crossterm(&
    q, R, d_shell, d_dead, d_surfactant,&
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sigR, mag_SLDcore, mag_SLDshell, &
    xi, sin2alpha, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, d_shell, d_dead, d_surfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant
    double precision, intent(in) :: SLDmatrix, sigR
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell
    double precision, intent(in) :: xi, sin2alpha
    integer, intent(in) :: Nq

    integer, parameter :: Np=9
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: Rmin, Rmax

    p = (/R, d_shell, d_dead, d_surfactant,&
      SLDcore, SLDshell, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_size_distribution(q(iq), p, Np, &
          1, Rmin, Rmax, sigR, &
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

        p_mag = (/&
          p(1), p(2), p(3), p(4),&
          mag_SLDcore, mag_SLDshell, 0d0, 0d0, 0d0&
        /)

        call crossterm(&
          q, xi, sin2alpha, p,&
          p_mag, p_multi, p_multi,&
          Np, Np, call_magnetic &
        )
      end function call_magnetic
  end subroutine p_crossterm

  subroutine sld(&
    R, d_shell, d_surfactant, &
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, x, sld_array)
    double precision, intent(in) :: R, d_shell, d_surfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix

    integer, parameter :: Nx=8
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = R + d_shell
    x(5) = R + d_shell
    x(6) = R + d_shell + d_surfactant
    x(7) = R + d_shell + d_surfactant
    x(8) = R + d_shell + d_surfactant + 100


    sld_array(1) = SLDcore
    sld_array(2) = SLDcore
    sld_array(3) = SLDshell
    sld_array(4) = SLDshell
    sld_array(5) = SLDsurfactant
    sld_array(6) = SLDsurfactant
    sld_array(7) = SLDmatrix
    sld_array(8) = SLDmatrix
  end subroutine sld

  subroutine magnetic_sld(&
    R, d_shell, d_dead, d_surfactant, &
    mag_SLDcore, mag_SLDshell, x, magsld_array)
    double precision, intent(in) :: R, d_shell, d_dead, d_surfactant
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell

    integer, parameter :: Nx=8
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: magsld_array

    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = R + d_shell - d_dead
    x(5) = R + d_shell - d_dead
    x(6) = R + d_shell - d_dead + d_surfactant
    x(7) = R + d_shell - d_dead + d_surfactant
    x(8) = R + d_shell - d_dead + d_surfactant + 100


    magsld_array(1) = mag_SLDcore
    magsld_array(2) = mag_SLDcore
    magsld_array(3) = mag_SLDshell
    magsld_array(4) = mag_SLDshell
    magsld_array(5) = 0
    magsld_array(6) = 0
    magsld_array(7) = 0
    magsld_array(8) = 0
  end subroutine magnetic_sld
end module sphere_css_dead
