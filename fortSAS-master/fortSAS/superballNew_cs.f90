module superball_new_cs
  use math
  implicit none
  contains
    double precision function amplitude_superball_oriented&
      (qx, qy, qz, r, pval, N_gl_order, x_leg, w_leg)
      double precision, intent(in) :: qx, qy, qz
      double precision, intent(in) :: r, pval

      integer, intent(in) :: N_gl_order
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg

      integer :: ix, iy
      double precision :: current_x, current_y, outer_integral, inner_integral
      double precision :: w_x, w_y, cos_x, gamma, zeta
      double precision :: twop, inv2p

      twop = 2d0 * pval
      inv2p = 1d0 / twop

      outer_integral = 0d0
      do ix=1, N_gl_order
        current_x = 5d-1 * (x_leg(ix) + 1d0) ! -1 to 1 for 0 to 1
        cos_x = cos(R * qx * current_x)
        gamma = ( 1d0 - abs(current_x)**(twop))**(inv2p)
        w_x = 5d-1 * w_leg(ix)
        inner_integral = 0d0
        do iy=1, N_gl_order
          current_y = 5d-1* gamma * (x_leg(iy) + 1d0) ! -1 to 1 for 0 to gamma
          w_y = 5d-1* gamma * w_leg(iy)
          zeta = ( 1d0 - abs(current_x)**(twop) - abs(current_y)**(twop) )**(inv2p)
          inner_integral = inner_integral + w_y*cos(r*qy*current_y)*sin(r*qz*zeta)
        end do
        outer_integral = outer_integral + w_x * cos_x * inner_integral
      end do
      ! Return absolute value squared times contrast
      amplitude_superball_oriented = 8d0*r**2 / qz * outer_integral
    end function amplitude_superball_oriented

    subroutine oriented_formfactor(&
      q, u_qx, u_qy, u_qz, R, d, p_val, p_shell,&
      SLDsuperball, SLDshell, SLDmatrix,&
      x_leg, w_leg, N_gl_order, &
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: u_qx, u_qy, u_qz
      double precision, intent(in) :: R, d, p_val, p_shell
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDmatrix
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
      integer, intent(in) :: Nq, N_gl_order

      double precision, dimension(Nq), intent(out) :: ff_intensity
      integer :: iq

      double precision :: amp_core, amp_shell

      !$omp parallel private(amp_core, amp_shell)
      !$omp do
      do iq=1, Nq
        amp_core = amplitude_superball_oriented&
          (q(iq)*u_qx, q(iq)*u_qy, q(iq)*u_qz,&
          r, p_val, N_gl_order, x_leg, w_leg)
        amp_shell = amplitude_superball_oriented&
          (q(iq)*u_qx, q(iq)*u_qy, q(iq)*u_qz,&
          r+d, p_shell, N_gl_order, x_leg, w_leg)

        ff_intensity(iq) = abs(&
          (SLDsuperball - SLDshell)*amp_core +&
          (SLDshell - SLDmatrix)*amp_shell)**2
      end do
      !$omp end do
      !$omp end parallel
    end subroutine oriented_formfactor

    subroutine formfactor(&
      q, R, d, p_val, p_shell, SLDsuperball, SLDshell, SLDmatrix, sigR, &
      x_herm, w_herm, x_leg, w_leg, N_gh_order, N_gl_order, &
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: R, d, p_val, p_shell
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDmatrix
      double precision, intent(in) :: sigR
      double precision, dimension(N_gh_order), intent(in) :: x_herm, w_herm
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
      integer, intent(in) :: Nq, N_gh_order, N_gl_order
      double precision, dimension(Nq), intent(out) :: ff_intensity

      double precision :: current_R, current_phi, current_th, sin_theta
      double precision :: u_qx, u_qy, u_qz
      double precision :: w_R_i, w_phi_i, w_theta_i
      double precision :: intsum_norm
      double precision, dimension(Nq) :: hff, intsum_ff

      integer :: iR, ith, iphi
      double precision :: integrand_prefacs

      intsum_ff = 0d0
      intsum_norm = 0d0
      if (sigR > tolerance) then
        do iR=1, N_gh_order
          ! gaussian hermite quadrature
          write(*,'(1a1,i0,A2,$)') char(13), iR * 100/N_gh_order, " %"

          current_R = R*dexp(sq2 * x_herm(iR) * sigR)
          w_R_i = w_herm(iR)
          do ith=1, N_gl_order
            current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
            u_qz = cos(current_th)
            sin_theta = sin(current_th)
            w_theta_i = pi/4d0 * w_leg(ith)

            do iphi=1, N_gl_order
              current_phi = pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
              w_phi_i =  pi/4d0 * w_leg(iphi)
              u_qx = cos(current_phi) * sin_theta
              u_qy = sin(current_phi) * sin_theta
              call oriented_formfactor(q, u_qx, u_qy, u_qz,&
                current_R, d, p_val, p_shell, &
                SLDsuperball, SLDshell, SLDmatrix,&
                x_leg, w_leg, N_gl_order,&
                Nq, hff)

              integrand_prefacs = w_R_i * w_phi_i * w_theta_i * sin_theta

              intsum_ff = intsum_ff + integrand_prefacs * hff
              intsum_norm = intsum_norm + integrand_prefacs
            end do
          end do
        end do
      else
        ! no size distribution
        do ith=1, N_gl_order
          current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
          u_qz = cos(current_th)
          sin_theta = sin(current_th)
          w_theta_i = pi/4d0 * w_leg(ith)

          do iphi=1, N_gl_order
            current_phi =  pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
            w_phi_i = pi/4d0 * w_leg(iphi)
            u_qx = cos(current_phi) * sin_theta
            u_qy = sin(current_phi) * sin_theta
            call oriented_formfactor(q, u_qx, u_qy, u_qz,&
              R, d, p_val, p_shell,&
              SLDsuperball, SLDshell, SLDmatrix,&
              x_leg, w_leg, N_gl_order,&
              Nq, hff)

            integrand_prefacs = w_phi_i * w_theta_i * sin_theta

            intsum_ff = intsum_ff + integrand_prefacs * hff
            intsum_norm = intsum_norm + integrand_prefacs
          end do
        end do
      end if
      ff_intensity = intsum_ff/intsum_norm
      print *, " Done."
    end subroutine formfactor

    subroutine magnetic_oriented_formfactor(&
      q,  u_qx, u_qy, u_qz, R, d, p_val, p_shell, SLDsuperball, SLDshell, SLDmatrix, &
      mag_SLDsuperball, mag_SLDshell, mag_SLDmatrix, xi, sin2alpha, plus_or_minus,&
      x_leg, w_leg, N_gl_order, &
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: u_qx, u_qy, u_qz
      double precision, intent(in) :: R, d, p_val, p_shell
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDmatrix
      double precision, intent(in) :: mag_SLDsuperball, mag_SLDshell, mag_SLDmatrix
      double precision, intent(in) :: xi, sin2alpha, plus_or_minus
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
      integer, intent(in) :: Nq, N_gl_order

      double precision, dimension(Nq), intent(out) :: ff_intensity
      integer :: iq

      double precision :: amp_core, amp_shell, f_nuclear, f_magnetic

      !$omp parallel private(amp_core, amp_shell, f_nuclear, f_magnetic)
      !$omp do
      do iq=1, Nq
        amp_core = amplitude_superball_oriented(&
          q(iq)*u_qx, q(iq)*u_qy, q(iq)*u_qz,&
          r, p_val, N_gl_order, x_leg, w_leg&
        )
        amp_shell = amplitude_superball_oriented(&
          q(iq)*u_qx, q(iq)*u_qy, q(iq)*u_qz,&
          r+d, p_shell, N_gl_order, x_leg, w_leg&
        )

        f_nuclear = &
          (SLDsuperball - SLDshell)*amp_core +&
          (SLDshell - SLDmatrix)*amp_shell
        f_magnetic = &
          (mag_SLDsuperball - mag_SLDshell)*amp_core +&
          (mag_SLDshell - mag_SLDmatrix)*amp_shell

        ff_intensity(iq) = abs(f_nuclear)**2 +&
            (abs(f_magnetic)**2 - plus_or_minus*2d0*xi*f_nuclear*f_magnetic)*sin2alpha
      end do
      !$omp end do
      !$omp end parallel
    end subroutine magnetic_oriented_formfactor

    subroutine magnetic_formfactor(&
      q, R, d, p_val, p_shell, SLDsuperball, SLDshell, SLDmatrix, sigR,&
      mag_SLDsuperball, mag_SLDshell, mag_SLDmatrix, xi, sin2alpha, plus_or_minus, &
      x_herm, w_herm, x_leg, w_leg, N_gh_order, N_gl_order, Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: R, d, p_val, p_shell
      double precision, intent(in) :: SLDsuperball, SLDshell, SLDmatrix
      double precision, intent(in) :: sigR
      double precision, intent(in) :: mag_SLDsuperball, mag_SLDshell, mag_SLDmatrix
      double precision, intent(in) :: xi, sin2alpha, plus_or_minus
      double precision, dimension(N_gh_order), intent(in) :: x_herm, w_herm
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
      integer, intent(in) :: Nq, N_gh_order, N_gl_order
      double precision, dimension(Nq), intent(out) :: ff_intensity

      double precision :: current_R, current_phi, current_th, sin_theta
      double precision :: u_qx, u_qy, u_qz
      double precision :: w_R_i, w_phi_i, w_theta_i
      double precision :: intsum_norm
      double precision, dimension(Nq) :: hff, intsum_ff

      integer :: iR, ith, iphi
      double precision :: integrand_prefacs

      intsum_ff = 0d0
      intsum_norm = 0d0
      if (sigR > tolerance) then
        do iR=1, N_gh_order
          ! gaussian hermite quadrature
          write(*,'(1a1,i0,A2,$)') char(13), iR * 100/N_gh_order, " %"

          current_R = R*dexp(sq2 * x_herm(iR) * sigR)
          w_R_i = w_herm(iR)
          do ith=1, N_gl_order
            current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
            u_qz = cos(current_th)
            sin_theta = sin(current_th)
            w_theta_i = pi/4d0 * w_leg(ith)

            do iphi=1, N_gl_order
              current_phi = pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
              w_phi_i =  pi/4d0 * w_leg(iphi)
              u_qx = cos(current_phi) * sin_theta
              u_qy = sin(current_phi) * sin_theta

              call magnetic_oriented_formfactor(q, u_qx, u_qy, u_qz,&
                current_R, d, p_val, p_shell,&
                SLDsuperball, SLDshell, SLDmatrix,&
                mag_SLDsuperball, mag_SLDshell, mag_SLDmatrix, xi, sin2alpha, plus_or_minus,&
                x_leg, w_leg, N_gl_order,&
                Nq, hff &
              )

              integrand_prefacs = w_R_i * w_phi_i * w_theta_i * sin_theta

              intsum_ff = intsum_ff + integrand_prefacs * hff
              intsum_norm = intsum_norm + integrand_prefacs
            end do
          end do
        end do
      else
        ! no size distribution
        do ith=1, N_gl_order
          current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
          u_qz = cos(current_th)
          sin_theta = sin(current_th)
          w_theta_i = pi/4d0 * w_leg(ith)

          do iphi=1, N_gl_order
            current_phi =  pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
            w_phi_i = pi/4d0 * w_leg(iphi)
            u_qx = cos(current_phi) * sin_theta
            u_qy = sin(current_phi) * sin_theta
            call magnetic_oriented_formfactor(q, u_qx, u_qy, u_qz,&
              R, d, p_val, p_shell,&
              SLDsuperball, SLDshell, SLDmatrix,&
              mag_SLDsuperball, mag_SLDshell, mag_SLDmatrix, xi, sin2alpha, plus_or_minus,&
              x_leg, w_leg, N_gl_order,&
              Nq, hff &
            )

            integrand_prefacs = w_phi_i * w_theta_i * sin_theta

            intsum_ff = intsum_ff + integrand_prefacs * hff
            intsum_norm = intsum_norm + integrand_prefacs
          end do
        end do
      end if
      ff_intensity = intsum_ff/intsum_norm
      print *, " Done."
    end subroutine magnetic_formfactor

    subroutine sld(R, d, SLDcube, SLDshell, SLDmatrix, x, sld_array)
      double precision, intent(in) :: R, d, SLDcube, SLDshell, SLDmatrix

      integer, parameter :: Nx=6
      double precision, dimension(Nx), intent(out) :: x
      double precision, dimension(Nx), intent(out) :: sld_array

      x(1) = 0
      x(2) = R
      x(3) = R
      x(4) = R + d
      x(5) = R + d
      x(6) = R + d + 100

      sld_array(1) = SLDcube
      sld_array(2) = SLDcube
      sld_array(3) = SLDshell
      sld_array(4) = SLDshell
      sld_array(5) = SLDmatrix
      sld_array(6) = SLDmatrix
    end subroutine sld
  endmodule superball_new_cs
