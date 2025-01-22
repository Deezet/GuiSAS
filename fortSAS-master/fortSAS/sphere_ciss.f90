module sphere_ciss
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

    double precision :: R, dint, dshell, dsurfactant
    double precision :: SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix
    double precision :: p_amp_core, p_amp_int, p_amp_shell, p_amp_surfactant
    R = p(1)
    dint = p(2)
    dshell = p(3)
    dsurfactant = p(4)
    SLDcore = p(5)
    SLDint = p(6)
    SLDshell = p(7)
    SLDsurfactant = p(8)
    SLDmatrix = p(9)
    p_amp_surfactant = p_sphere(q, &
            (/R+dint+dshell+dsurfactant, SLDsurfactant, SLDmatrix/), 3)
    p_amp_shell = p_sphere(q, (/R+dint+dshell, SLDshell, SLDsurfactant/), 3)
    p_amp_int = p_sphere(q, (/R+dint, SLDint, SLDshell/), 3)
    p_amp_core = p_sphere(q, (/R, SLDcore, SLDint/), 3)
    p_sphere_css = p_amp_surfactant + p_amp_shell + p_amp_int + p_amp_core
  end function p_sphere_css

  double precision function ff_sphere_css(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
    ff_sphere_css = &
            abs(p_sphere_css(q, p, Np))**2
  end function ff_sphere_css

  subroutine formfactor(&
    q, R, dint, dshell, dsurfactant, SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix, &
    sigR, sigDint, sigDshell, Nq, ff_intensity)

    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, dint, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix
    double precision, intent(in) :: sigR, sigDint, sigDshell
    integer, intent(in) :: Nq

    integer, parameter :: Np=9
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq
    double precision :: Rmin, Rmax
    double precision :: Dint_min, Dint_max
    double precision :: Dshell_min, Dshell_max

    p = (/R, dint, dshell, dsurfactant, SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    call get_cutoff_lognormal(dint, sigDint, Dint_min, Dint_max)
    call get_cutoff_lognormal(dshell, sigDshell, Dshell_min, Dshell_max)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_three_size_distributions(&
        q(iq), p, Np, &
        1, Rmin, Rmax, sigR, &
        2, Dint_min, Dint_max, sigDint, &
        3, Dshell_min, Dshell_max, sigDshell, &
        ff_sphere_css, lognormal, lognormal, lognormal,&
        ff_intensity(iq) &
      )
    end do
    !$omp end do
    !$omp end parallel
  end subroutine formfactor

  subroutine magnetic_formfactor(&
    q, R, dint, dshell, dsurfactant, SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix, &
    sigR, sigDint, sigDshell, &
    mag_SLDcore, mag_SLDint, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R, dint, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix
    double precision, intent(in) :: sigR, sigDint, sigDshell
    double precision, intent(in) :: mag_SLDcore, mag_SLDint, mag_SLDshell, mag_SLDsurfactant
    double precision, intent(in) :: mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq

    integer, parameter :: Np=9
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_intensity

    integer :: iq
    double precision :: Rmin, Rmax
    double precision :: Dint_min, Dint_max
    double precision :: Dshell_min, Dshell_max

    p = (/R, dint, dshell, dsurfactant, SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix/)
    call get_cutoff_lognormal(R, sigR, Rmin, Rmax)
    call get_cutoff_lognormal(dint, sigDint, Dint_min, Dint_max )
    call get_cutoff_lognormal(dshell, sigDshell, &
                              Dshell_min, Dshell_max)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_three_size_distributions(&
        q(iq), p, Np, &
        1, Rmin, Rmax, sigR, &
        2, Dint_min, Dint_max, sigDshell, &
        3, Dshell_min, Dshell_max, sigDshell, &
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
          p(1), p(2), p(3), &
          mag_SLDcore, mag_SLDint, mag_SLDshell, mag_SLDsurfactant,&
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
    R, dint, dshell, dsurfactant, &
    SLDcore, SLDint, SLDshell, SLDsurfactant, SLDmatrix, x, sld_array)
    double precision, intent(in) :: R, dint, dshell, dsurfactant
    double precision, intent(in) :: SLDcore, SLDint, SLDshell, SLDsurfactant
    double precision, intent(in) :: SLDmatrix

    integer, parameter :: Nx=10
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    x(1) = 0
    x(2) = R
    x(3) = R
    x(4) = R + dint
    x(5) = R + dint
    x(6) = R + dint + dshell
    x(7) = R + dint + dshell
    x(8) = R + dint + dshell + dsurfactant
    x(9) = R + dint + dshell + dsurfactant
    x(10) = R + dint + dshell + dsurfactant + 100

    sld_array(1) = SLDcore
    sld_array(2) = SLDcore
    sld_array(3) = SLDint
    sld_array(4) = SLDint
    sld_array(5) = SLDshell
    sld_array(6) = SLDshell
    sld_array(7) = SLDsurfactant
    sld_array(8) = SLDsurfactant
    sld_array(9) = SLDmatrix
    sld_array(10) = SLDmatrix
  end subroutine sld
end module sphere_ciss