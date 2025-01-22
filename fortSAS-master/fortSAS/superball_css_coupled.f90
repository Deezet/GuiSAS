module superball_css_coupled
  use math
  implicit none
  contains
    double precision function amplitude_superball_oriented(q, p, Np)
      double precision, intent(in) :: q
      double precision, dimension(Np), intent(in) :: p
      integer, intent(in) :: Np

      double precision, dimension(Np) :: help_p
      double precision :: R, p_val, alpha, beta
      double precision :: SLDsuperball, SLDmatrix

      double precision :: qx, qy, qz
      double precision :: real_ff_superball

      R = p(1)
      p_val = p(2)
      SLDsuperball = p(3)
      SLDmatrix = p(4)
      alpha = p(5)
      beta = p(6)
      qx = q * cos(alpha)*sin(beta)
      qy = q * sin(alpha)*sin(beta)
      qz = q * cos(beta)

      ! following parameters are being passed on:
      help_p = (/qx, qy, qz, R, 2d0*p_val/)
      ! Calculate real part of formfactor amplitude using dqag
      call twodim_integral_variable_bounds(-1d0, 1d0, &
                          ymin_func, ymax_func, &
                          help_p, 5, &
                          ff_real_part, &
                          real_ff_superball)
      ! Return absolute value squared times contrast
      amplitude_superball_oriented =&
        (SLDsuperball - SLDmatrix) * 2d0 * R**2 / qz * real_ff_superball
      contains
        double precision function ymin_func(x, p, Np)
          double precision, intent(in) :: x
          double precision, dimension(Np), intent(in) :: p
          integer, intent(in) :: Np
          ymin_func = -(1d0-abs(x)**p(5))**(1d0/p(5))
        end function

        double precision function ymax_func(x, p, Np)
          double precision, intent(in) :: x
          double precision, dimension(Np), intent(in) :: p
          integer, intent(in) :: Np
          ymax_func = (1d0-abs(x)**p(5))**(1d0/p(5))
        end function

        double precision function ff_real_part(x,y, p, Np)
          double precision, intent(in) :: x, y
          double precision, dimension(Np), intent(in) :: p
          integer, intent(in) :: Np

          double precision :: qx, qy, qz, R, n_exp, zval
          qx = p(1)
          qy = p(2)
          qz = p(3)
          R = p(4)
          n_exp = p(5)
          zval = abs(1d0 - abs(x)**n_exp - abs(y)**n_exp)**(1d0/n_exp)
          ff_real_part = sin(qz*R*zval)*cos(qy*R*y + qx*R*x)
        end function
    end function amplitude_superball_oriented

    double precision function p_superball_css(q, p, Np)
      double precision, intent(in) :: q
      double precision, dimension(Np), intent(in) :: p
      integer, intent(in) :: Np

      double precision :: particleSize, R, dShell, dSurfactant, p_val
      double precision :: SLDcore, SLDshell, SLDsurfactant, SLDmatrix
      double precision :: alpha, beta
      double precision :: ff_amp_core, ff_amp_shell, ff_amp_surfactant
      particleSize = p(1)
      dShell = p(2)
      dSurfactant = p(3)
      p_val = p(4)
      SLDcore = p(5)
      SLDshell = p(6)
      SLDsurfactant = p(7)
      SLDmatrix = p(8)
      alpha = p(9)
      beta = p(10)
      R = particleSize - dShell
      if (R < 0) then
        R = 0
      end if

      ff_amp_surfactant = amplitude_superball_oriented(q,&
        (/R+dShell+dSurfactant, p_val, SLDsurfactant, SLDmatrix, alpha, beta/),&
        6)
      ff_amp_shell = amplitude_superball_oriented(q, (/R+dShell, p_val, SLDshell, SLDsurfactant, alpha, beta/), 6)
      ff_amp_core = amplitude_superball_oriented(q, (/R, p_val, SLDcore, SLDshell, alpha, beta/), 6)
      p_superball_css = ff_amp_surfactant + ff_amp_shell + ff_amp_core
    end function p_superball_css

    subroutine oriented_formfactor(&
      q, particleSize, dShell, dSurfactant, p_val,&
      SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, alpha, beta, Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: particleSize, dShell, dSurfactant, p_val
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix
      double precision, intent(in) :: alpha, beta
      integer, intent(in) :: Nq

      integer, parameter :: Np=10
      double precision, dimension(Np) :: p
      double precision, dimension(Nq), intent(out) :: ff_intensity
      integer :: iq

      p = (/particleSize, dShell, dSurfactant, p_val, SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, alpha, beta/)

      !$omp parallel
      !$omp do
      do iq=1, Nq
        ff_intensity(iq) = abs(p_superball_css(q(iq), p, Np))**2
      end do
      !$omp end do
      !$omp end parallel
    end subroutine oriented_formfactor

    subroutine formfactor(&
      q, particleSize, dShell, dSurfactant, p_val, SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, sigR, &
      x_herm, w_herm, x_leg, w_leg, N_gh_order, N_gl_order, &
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: particleSize, dShell, dSurfactant, p_val
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix
      double precision, intent(in) :: sigR
      double precision, dimension(N_gh_order), intent(in) :: x_herm, w_herm
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
      integer, intent(in) :: Nq, N_gh_order, N_gl_order
      double precision, dimension(Nq), intent(out) :: ff_intensity

      integer, parameter :: Np=8

      double precision :: current_particleSize, current_phi, current_th, current_sinth
      double precision :: w_R_i, w_phi_i, w_theta_i
      double precision :: intsum_norm
      double precision, dimension(Nq) :: hff, intsum_ff

      integer :: iR, ith, iphi
      double precision :: integrand_prefacs

      intsum_ff = 0d0
      intsum_norm = 0d0
      n_integration_cuts = 1
      if (sigR > tolerance) then
        do iR=1, N_gh_order
          ! gaussian hermite quadrature
          write(*,'(1a1,i0,A2,$)') char(13), iR * 100/N_gh_order, " %"

          current_particleSize = particleSize*dexp(sq2 * x_herm(iR) * sigR)
          w_R_i = w_herm(iR)
          do ith=1, N_gl_order
            current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
            current_sinth = sin(current_th)
            w_theta_i = pi/4d0 * w_leg(ith)

            do iphi=1, N_gl_order
              current_phi = pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
              w_phi_i =  pi/4d0 * w_leg(iphi)
              call oriented_formfactor(q, current_particleSize, dShell, dSurfactant, p_val,&
                SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, &
                current_phi, current_th, Nq, hff)

              integrand_prefacs = w_R_i * w_phi_i * w_theta_i * current_sinth

              intsum_ff = intsum_ff + integrand_prefacs * hff
              intsum_norm = intsum_norm + integrand_prefacs
            end do
          end do
        end do
      else
        ! no size distribution
        do ith=1, N_gl_order
          current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
          current_sinth = sin(current_th * deg)
          w_theta_i = pi/4d0 * w_leg(ith)

          do iphi=1, N_gl_order
            current_phi =  pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
            w_phi_i = pi/4d0 * w_leg(iphi)
            call oriented_formfactor(q, particleSize, dShell, dSurfactant, p_val,&
              SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, &
              current_phi, current_th, Nq, hff)

            integrand_prefacs = w_phi_i * w_theta_i * current_sinth

            intsum_ff = intsum_ff + integrand_prefacs * hff
            intsum_norm = intsum_norm + integrand_prefacs
          end do
        end do
      end if
      ff_intensity = intsum_ff/intsum_norm
      print *, " Done."
    end subroutine formfactor

    subroutine magnetic_oriented_formfactor(&
      q, particleSize, dShell, dSurfactant, p_val, SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, alpha, beta, &
      mag_SLDsuperball, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix, xi, sin2alpha, plus_or_minus,&
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: particleSize, dShell, dSurfactant, p_val
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: mag_SLDsuperball, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix
      double precision, intent(in) :: xi, sin2alpha, plus_or_minus
      integer, intent(in) :: Nq

      integer, parameter :: Np=10
      double precision, dimension(Np) :: p
      double precision, dimension(Nq), intent(out) :: ff_intensity
      integer :: iq

      p = (/particleSize, dShell, dSurfactant, p_val, SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, alpha, beta/)

      !$omp parallel
      !$omp do
      do iq=1, Nq
        ff_intensity(iq) = call_magnetic(q(iq), p, Np)
      end do
      !$omp end do
      !$omp end parallel

      contains
      double precision function call_magnetic(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np

        double precision, dimension(Np) :: p_mag
        p_mag = (/p(1), p(2), p(3), p(4), mag_SLDsuperball, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix, alpha, beta/)
        call magnetic_scattering(&
          q, xi, sin2alpha, plus_or_minus, p, p_mag,&
          p_superball_css,&
          p_superball_css,&
          Nq, Np, Np, call_magnetic)
      end function call_magnetic
    end subroutine magnetic_oriented_formfactor

    subroutine magnetic_formfactor(&
      q, particleSize, dShell, dSurfactant, p_val, SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, sigR,&
      mag_SLDsuperball, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix, xi, sin2alpha, plus_or_minus, &
      x_herm, w_herm, x_leg, w_leg, N_gh_order, N_gl_order, &
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: particleSize, dShell, dSurfactant, p_val
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix
      double precision, intent(in) :: sigR
      double precision, intent(in) :: mag_SLDsuperball, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix
      double precision, intent(in) :: xi, sin2alpha, plus_or_minus
      double precision, dimension(N_gh_order), intent(in) :: x_herm, w_herm
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
      integer, intent(in) :: Nq, N_gh_order, N_gl_order
      double precision, dimension(Nq), intent(out) :: ff_intensity

      integer, parameter :: Np=6

      double precision :: current_particleSize, current_phi, current_th, current_sinth
      double precision :: w_R_i, w_phi_i, w_theta_i

      double precision :: intsum_norm
      double precision, dimension(Nq) :: hff, intsum_ff

      integer :: iR, ith, iphi
      double precision :: integrand_prefacs

      intsum_ff = 0d0
      intsum_norm = 0d0
      n_integration_cuts = 1
      if (sigR > tolerance) then
        do iR=1, N_gh_order
          ! gaussian hermite quadrature
          write(*,'(1a1,i0,A2,$)') char(13), iR * 100/N_gh_order, " %"

          current_particleSize = particleSize*dexp(sq2 * x_herm(iR) * sigR)
          w_R_i = w_herm(iR)

          do ith=1, N_gl_order
            current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
            current_sinth = sin(current_th)
            w_theta_i = pi/4d0 * w_leg(ith)

            do iphi=1, N_gl_order
              current_phi =  pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
              w_phi_i = pi/4d0 * w_leg(iphi)

              call magnetic_oriented_formfactor(&
                q, current_particleSize, dShell, dSurfactant, p_val,&
                SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, current_phi, current_th, &
                mag_SLDsuperball, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix, xi, sin2alpha, plus_or_minus, &
                Nq, hff &
              )

              integrand_prefacs = w_R_i * w_phi_i * w_theta_i * current_sinth

              intsum_ff = intsum_ff + integrand_prefacs * hff
              intsum_norm = intsum_norm + integrand_prefacs
            end do
          end do
        end do
      else
        ! no size distribution
        do ith=1, N_gl_order
          current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
          current_sinth = sin(current_th)
          w_theta_i = pi/4d0 * w_leg(ith)

          do iphi=1, N_gl_order
            current_phi =  pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
            w_phi_i = pi/4d0 * w_leg(iphi)

            call magnetic_oriented_formfactor(&
              q, particleSize, dShell, dSurfactant, p_val,&
              SLDsuperball, SLDshell, SLDsurfactant, SLDmatrix, current_phi, current_th, &
              mag_SLDsuperball, mag_SLDshell, mag_SLDsurfactant, mag_SLDmatrix, xi, sin2alpha, plus_or_minus, &
              Nq, hff &
            )

            integrand_prefacs = w_phi_i * w_theta_i * current_sinth

            intsum_ff = intsum_ff + integrand_prefacs * hff
            intsum_norm = intsum_norm + integrand_prefacs
          end do
        end do
      end if
      ff_intensity = intsum_ff/intsum_norm
      print *, " Done."
    end subroutine magnetic_formfactor

    subroutine sld(particleSize, dShell, dSurfactant, SLDcube, SLDshell, SLDsurfactant, SLDmatrix, x, sld_array)
      double precision, intent(in) :: particleSize, dShell, dSurfactant
      double precision, intent(in) :: SLDcube, SLDshell, SLDsurfactant, SLDmatrix
      double precision :: R

      integer, parameter :: Nx=8
      double precision, dimension(Nx), intent(out) :: x
      double precision, dimension(Nx), intent(out) :: sld_array

      R = particleSize - dShell
      if (R < 0) then
        R = 0
      end if

      x(1) = 0
      x(2) = R
      x(3) = R
      x(4) = R + dShell
      x(5) = R + dShell
      x(6) = R + dShell + dSurfactant
      x(7) = R + dShell + dSurfactant
      x(8) = R + dShell + dSurfactant + 100

      sld_array(1) = SLDcube
      sld_array(2) = SLDcube
      sld_array(3) = SLDshell
      sld_array(4) = SLDshell
      sld_array(5) = SLDsurfactant
      sld_array(6) = SLDsurfactant
      sld_array(7) = SLDmatrix
      sld_array(8) = SLDmatrix
    end subroutine sld
  endmodule superball_css_coupled