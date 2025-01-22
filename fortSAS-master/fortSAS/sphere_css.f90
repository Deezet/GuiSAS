module sphere_css
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

  double precision function p_sphere_css(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: R, dshell, dsurfactant
    double precision :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
    double precision :: p_amp_core, p_amp_shell, p_amp_surfactant
    R = p(1)
    dshell = p(2)
    dsurfactant = p(3)
    SLDcore = p(4)
    SLDshell = p(5)
    SLDsurfactant = p(6)
    SLDmatrix = p(7)
    p_amp_surfactant = p_sphere(q, &
            (/R+dshell+dsurfactant, SLDsurfactant, SLDmatrix/), 3)
    p_amp_shell = p_sphere(q, (/R+dshell, SLDshell, SLDsurfactant/), 3)
    p_amp_core = p_sphere(q, (/R, SLDcore, SLDshell/), 3)
    p_sphere_css = p_amp_surfactant + p_amp_shell + p_amp_core
  end function p_sphere_css

  double precision function ff_sphere_css(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
    ff_sphere_css = &
            abs(p_sphere_css(q, p, Np))**2
  end function ff_sphere_css

  subroutine formfactor(&
    q, R, dshell, dsurfactant, SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sigR, sigDshell, sigDsurfactant, Nq, ff_intensity)

    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
    double precision, intent(in) :: sigR, sigDshell, sigDsurfactant
    integer, intent(in) :: Nq

    integer, parameter :: Np=7
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq
    double precision :: Rmin, Rmax
    double precision :: Dshell_min, Dshell_max
    double precision :: Dsurfactant_min, Dsurfactant_max

    p = (/R, dshell, dsurfactant, SLDcore, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    call get_cutoff_lognormal(dshell, sigDshell, Dshell_min, Dshell_max)
    call get_cutoff_lognormal(dsurfactant, sigDsurfactant, &
                              Dsurfactant_min, Dsurfactant_max)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_three_size_distributions(&
        q(iq), p, Np, &
        1, Rmin, Rmax, sigR, &
        2, Dshell_min, Dshell_max, sigDshell, &
        3, Dsurfactant_min, Dsurfactant_max, sigDsurfactant, &
        ff_sphere_css, lognormal, lognormal, lognormal,&
        ff_intensity(iq) &
      )
    end do
    !$omp end do
    !$omp end parallel
  end subroutine formfactor

  subroutine magnetic_formfactor(&
    q, R, dshell, dsurfactant, SLDcore, SLDshell, SLDsurfactant, SLDmatrix, &
    sigR, sigDshell, sigDsurfactant, &
    mag_SLDcore, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
    double precision, intent(in) :: sigR, sigDshell, sigDsurfactant
    double precision, intent(in) :: mag_SLDcore, mag_SLDshell, mag_SLDsurfactant
    double precision, intent(in) :: mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq

    integer, parameter :: Np=7
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: Rmin, Rmax
    double precision :: Dshell_min, Dshell_max
    double precision :: Dsurfactant_min, Dsurfactant_max

    p = (/R, dshell, dsurfactant, SLDcore, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    call get_cutoff_lognormal(dshell, sigDshell, Dshell_min, Dshell_max)
    call get_cutoff_lognormal(dsurfactant, sigDsurfactant, &
                              Dsurfactant_min, Dsurfactant_max)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_three_size_distributions(&
        q(iq), p, Np, &
        1, Rmin, Rmax, sigR, &
        2, Dshell_min, Dshell_max, sigDshell, &
        3, Dsurfactant_min, Dsurfactant_max, sigDsurfactant, &
        call_magnetic, lognormal, lognormal, lognormal,&
        ff_intensity(iq) &
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

        p_mag = (/&
          p(1), p(2), p(3),&
          mag_SLDcore, mag_SLDshell, mag_SLDsurfactant,&
          mag_SLDmatrix &
        /)
        call magnetic_scattering(&
          q, xi, sin2alpha, plus_or_minus, p,&
          p_mag, p_sphere_css, p_sphere_css,&
          Nq, Np, Np, call_magnetic &
        )
      end function call_magnetic
  end subroutine magnetic_formfactor

  subroutine sld(&
    R, dshell, dsurfactant, &
    SLDcore, SLDshell, SLDsurfactant, SLDmatrix, x, sld_array)
    double precision, intent(in) :: R, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDshell, SLDsurfactant
    double precision, intent(in) :: SLDmatrix

    integer, parameter :: Nx=8
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = R + dshell
    x(5) = R + dshell
    x(6) = R + dshell + dsurfactant
    x(7) = R + dshell + dsurfactant
    x(8) = R + dshell + dsurfactant + 100

    sld_array(1) = SLDcore
    sld_array(2) = SLDcore
    sld_array(3) = SLDshell
    sld_array(4) = SLDshell
    sld_array(5) = SLDsurfactant
    sld_array(6) = SLDsurfactant
    sld_array(7) = SLDmatrix
    sld_array(8) = SLDmatrix
  end subroutine sld
end module sphere_css
