module sphere_linhulls
use math
implicit none
contains
  double precision function p_sphere_lin_hull(q, p, Np)
    ! Sphere Formfacotr Amplitude
    ! Parameters: R, SLDsphere, SLDmatrix
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: R1, R2, R3
    double precision :: dR, dsurf, SLDcore, SLDhull, SLDsurf, SLDmatrix
    double precision :: m, qR1, qR2, qR3
    R1 = p(1)
    dR = p(2)
    dsurf = p(3)
    SLDcore = p(4)
    SLDhull = p(5)
    SLDsurf = p(6)
    SLDmatrix = p(7)

    R2 = R1 + dR
    R3 = R2 + dsurf
    if (dR /= 0) then
      m = (SLDhull - SLDcore) / dR
    else
      m = 0d0
    end if
    qR1 = q*R1
    qR2 = q*R2
    qR3 = q*R3
    if (q /= 0) then
      p_sphere_lin_hull = 4d0*pi/q**3*(SLDcore - SLDsurf)*&
            (sin(qR2) - qR2*cos(qR2)) +&
            4d0*pi/q**3*(SLDsurf - SLDmatrix)*&
            (sin(qR3) - qR3*cos(qR3)) +&
            4d0*pi*m/q**4*&
            ((2*qR2 - qR1)*sin(qR2)+&
             (qR2*(qR1-qR2) + 2d0)*cos(qR2) - qR1*sin(qR1) -&
             2*cos(qR1))
    else
      p_sphere_lin_hull = 4d0/3d0 *pi*(SLDcore - SLDsurf)*R2**3/3 +&
                4d0/3d0 *pi*(SLDsurf - SLDmatrix)*R3**3/3 +&
                4d0/3d0 *pi*m*(R2**4/4d0 + R1**4/12d0 + R1*R2**3/6d0)
    end if
  end function p_sphere_lin_hull

  double precision function ff_sphere_lin_hull(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
    ff_sphere_lin_hull = abs(p_sphere_lin_hull(q, p, Np))**2
  end function ff_sphere_lin_hull

  subroutine formfactor(q, R1, dR, dsurf,&
            SLDcore, SLDhull, SLDsurf, SLDmatrix,&
            sigR1, sigdR,&
            Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R1, dR, dsurf
    double precision, intent(in) :: SLDcore, SLDhull, SLDsurf, SLDmatrix
    double precision, intent(in) :: sigR1, sigdR
    integer, intent(in) :: Nq

    integer, parameter :: Np=7
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq
    double precision :: R1min, R1max, dRmin, dRmax
    p = (/R1, dR, dsurf, SLDcore, SLDhull, SLDsurf, SLDmatrix/)
    call get_cutoff_lognormal(R1, sigR1, R1min, R1max)
    call get_cutoff_lognormal(dR, sigdR, dRmin, dRmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(q(iq), p, Np, &
              1, R1min, R1max, sigR1, &
              2, dRmin, dRmax, sigdR, &
              ff_sphere_lin_hull, lognormal, lognormal, ff_intensity(iq))
    end do
    !$omp end do
    !$omp end parallel
  end subroutine formfactor

  subroutine magnetic_formfactor(&
    q, R1, dR, dsurf, SLDcore, SLDhull, SLDsurf, SLDmatrix,&
    sigR1, sigdR, mag_SLDcore, mag_SLDhull, mag_SLDsurf, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus, Nq, ff_out)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R1, dR, dsurf
    double precision, intent(in) :: SLDcore, SLDhull, SLDsurf, SLDmatrix
    double precision, intent(in) :: sigR1, sigdR
    double precision, intent(in) :: mag_SLDcore, mag_SLDhull, mag_SLDsurf, mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq

    integer, parameter :: Np=7
    double precision, dimension(Np) :: p

    double precision, dimension(Nq), intent(out) :: ff_out

    integer :: iq
    double precision :: R1min, R1max, dRmin, dRmax
    p = (/R1, dR, dsurf, SLDcore, SLDhull, SLDsurf, SLDmatrix/)
    call get_cutoff_lognormal(R1, sigR1, R1min, R1max)
    call get_cutoff_lognormal(dR, sigdR, dRmin, dRmax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
      call integrate_two_size_distributions(q(iq), p, Np, &
            1, R1min, R1max, sigR1, &
            2, dRmin, dRmax, sigdR, &
            call_magnetic, lognormal, lognormal,&
            ff_out(iq))
    end do
    !$omp end do
    !$omp end parallel
    contains
      double precision function call_magnetic(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np

        double precision, dimension(Np) :: p_mag
        p_mag = (/p(1), p(2), dsurf,&
             mag_SLDcore, mag_SLDhull, mag_SLDsurf, mag_SLDmatrix/)
        call magnetic_scattering(q, xi, sin2alpha, plus_or_minus, p,&
          p_mag, p_sphere_lin_hull, p_sphere_lin_hull,&
          Nq, Np, Np, call_magnetic)
      end function call_magnetic
  end subroutine magnetic_formfactor

  subroutine sld(&
    R1, dR, dsurf, &
    SLDcore, SLDhull, SLDsurf, SLDmatrix, x, sld_array)
    double precision, intent(in) :: R1, dR, dsurf
    double precision, intent(in) :: SLDcore, SLDhull, SLDsurf
    double precision, intent(in) :: SLDmatrix

    integer, parameter :: Nx=8
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    x(1) = 0
    x(2) = R1
    x(3) = R1 + dR
    x(4) = R1 + dR
    x(5) = R1 + dR
    x(6) = R1 + dR + dsurf
    x(7) = R1 + dR + dsurf
    x(8) = R1 + dR + dsurf + 100

    sld_array(1) = SLDcore
    sld_array(2) = SLDcore
    sld_array(3) = SLDhull
    sld_array(4) = SLDhull
    sld_array(5) = SLDsurf
    sld_array(6) = SLDsurf
    sld_array(7) = SLDmatrix
    sld_array(8) = SLDmatrix
  end subroutine sld

  subroutine magnetic_sld(&
    R1, dR, dsurf, &
    mag_SLDcore, mag_SLDhull, x, magsld_array)
    double precision, intent(in) :: R1, dR, dsurf
    double precision, intent(in) :: mag_SLDcore, mag_SLDhull

    integer, parameter :: Nx=8
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: magsld_array

    x(1) = 0
    x(2) = R1
    x(3) = R1 + dR
    x(4) = R1 + dR
    x(5) = R1 + dR
    x(6) = R1 + dR + dsurf
    x(7) = R1 + dR + dsurf
    x(8) = R1 + dR + dsurf + 100

    magsld_array(1) = mag_SLDcore
    magsld_array(2) = mag_SLDcore
    magsld_array(3) = mag_SLDhull
    magsld_array(4) = mag_SLDhull
    magsld_array(5) = 0
    magsld_array(6) = 0
    magsld_array(7) = 0
    magsld_array(8) = 0
  end subroutine magnetic_sld
end module sphere_linhulls
