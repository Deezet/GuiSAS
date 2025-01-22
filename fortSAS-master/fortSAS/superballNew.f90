module superball_new
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
      q, u_qx, u_qy, u_qz, R, p_val,&
      x_leg, w_leg, N_gl_order, &
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: u_qx, u_qy, u_qz
      double precision, intent(in) :: R, p_val
      double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
      integer, intent(in) :: Nq, N_gl_order

      double precision, dimension(Nq), intent(out) :: ff_intensity
      integer :: iq

      !$omp parallel
      !$omp do
      do iq=1, Nq
        ff_intensity(iq) = abs(&
        amplitude_superball_oriented&
        (q(iq)*u_qx, q(iq)*u_qy, q(iq)*u_qz,&
        r, p_val, N_gl_order, x_leg, w_leg))**2
      end do
      !$omp end do
      !$omp end parallel
    end subroutine oriented_formfactor

    subroutine formfactor(&
      q, R, p_val, SLDsuperball, SLDmatrix, sigR, &
      x_herm, w_herm, x_leg, w_leg, N_gh_order, N_gl_order, &
      Nq, ff_intensity)
      double precision, dimension(Nq), intent(in) :: q
      double precision, intent(in) :: R, p_val, SLDsuperball, SLDmatrix
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
                current_R, p_val,&
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
              R, p_val,&
              x_leg, w_leg, N_gl_order,&
              Nq, hff)

            integrand_prefacs = w_phi_i * w_theta_i * sin_theta

            intsum_ff = intsum_ff + integrand_prefacs * hff
            intsum_norm = intsum_norm + integrand_prefacs
          end do
        end do
      end if
      ff_intensity = (SLDsuperball - SLDmatrix)**2 * intsum_ff/intsum_norm
      print *, " Done."
    end subroutine formfactor

    subroutine sld(R, SLDcube, SLDmatrix, x, sld_array)
      double precision, intent(in) :: R, SLDcube, SLDmatrix

      integer, parameter :: Nx=4
      double precision, dimension(Nx), intent(out) :: x
      double precision, dimension(Nx), intent(out) :: sld_array

      x(1) = 0
      x(2) = R
      x(3) = R
      x(4) = R + 100

      sld_array(1) = SLDcube
      sld_array(2) = SLDcube
      sld_array(3) = SLDmatrix
      sld_array(4) = SLDmatrix
    end subroutine sld
  endmodule superball_new