module cube
use math
implicit none
contains
  double precision function p_cube(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np

    double precision :: qa_2, qxa_2, qya_2, qza_2
    double precision :: six, siy, siz
    double precision :: a, SLDcube, SLDmatrix, phi, theta
    double precision :: costheta, sintheta, cosphi, sinphi

    a = p(1)
    SLDcube = p(2)
    SLDmatrix = p(3)
    phi = p(4)
    theta = p(5)

    cosphi = cos(phi)
    sinphi = sin(phi)
    costheta = cos(theta)
    sintheta = sin(theta)

    qa_2 = q*a/2d0
    qxa_2 = qa_2*sintheta*cosphi
    qya_2 = qa_2*sintheta*sinphi
    qza_2 = qa_2*costheta
    if (qxa_2 == 0d0) then
      six = 1d0
    else
      six = sin(qxa_2)/qxa_2
    end if

    if (qya_2 == 0d0) then
      siy = 1d0
    else
      siy = sin(qya_2)/qya_2
    end if

    if (qza_2 == 0d0) then
      siz = 1d0
    else
      siz = sin(qza_2)/qza_2
    end if
    p_cube = a**3 *(SLDcube - SLDmatrix)*six*siy*siz
  end function p_cube

  subroutine oriented_formfactor(&
    q, a, SLDcube, SLDmatrix, phi, theta, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: a, SLDcube, SLDmatrix
    double precision, intent(in) :: phi, theta
    integer, intent(in) :: Nq

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq

    p = (/a, SLDcube, SLDmatrix, phi, theta/)
    !$omp parallel
    !$omp do
    do iq=1, Nq
      ff_intensity(iq) = abs(p_cube(q(iq), p, Np))**2
    end do
    !$omp end do
    !$omp end parallel
  end subroutine oriented_formfactor

  subroutine formfactor(&
    q, a, SLDcube, SLDmatrix, siga,&
    x_herm, w_herm, x_leg, w_leg, N_gh_order, N_gl_order, Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: a, SLDcube, SLDmatrix
    double precision, intent(in) :: siga
    double precision, dimension(N_gh_order), intent(in) :: x_herm, w_herm
    double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
    integer, intent(in) :: N_gh_order, N_gl_order, Nq

    double precision, dimension(Nq), intent(out) :: ff_intensity
    double precision :: current_a, current_phi, current_th, current_sinth
    double precision :: w_a_i, w_phi_i, w_theta_i
    double precision :: intsum_norm
    double precision, dimension(Nq) :: hff, intsum_ff

    integer :: ia, ith, iphi
    double precision :: integrand_prefacs

    intsum_ff = 0d0
    intsum_norm = 0d0
    if (siga > tolerance) then
      do ia=1, N_gh_order
        current_a = a*dexp(sq2 * x_herm(ia) * siga)
        w_a_i = w_herm(ia)
        do ith=1, N_gl_order
          current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
          current_sinth = sin(current_th)
          w_theta_i = pi/4d0 * w_leg(ith)

          do iphi=1, N_gl_order
            current_phi = pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
            w_phi_i =  pi/4d0 * w_leg(iphi)
            call oriented_formfactor(q, current_a,&
              SLDcube, SLDmatrix, current_phi, current_th, Nq, hff)

            integrand_prefacs = w_a_i * w_phi_i * w_theta_i * current_sinth

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
          call oriented_formfactor(q, a,&
            SLDcube, SLDmatrix, current_phi, current_th, Nq, hff)

          integrand_prefacs = w_phi_i * w_theta_i * current_sinth

          intsum_ff = intsum_ff + integrand_prefacs * hff
          intsum_norm = intsum_norm + integrand_prefacs
        end do
      end do
    end if
    ff_intensity = intsum_ff/intsum_norm
  end subroutine formfactor

  subroutine magnetic_oriented_formfactor(&
    q, a, SLDcube, SLDmatrix, phi, theta, &
    mag_SLDcube, mag_SLDmatrix, &
    xi, sin2alpha, plus_or_minus, Nq, ff_out)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: a, SLDcube, SLDmatrix, phi, theta
    double precision, intent(in) :: mag_SLDcube, mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    integer, intent(in) :: Nq

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_out

    integer :: iq

    p = (/a, SLDcube, SLDmatrix, phi, theta/)
    !$omp parallel
    !$omp do
    do iq=1, Nq
      ff_out(iq) = call_magnetic(q(iq), p, Np)
    end do
    !$omp end do
    !$omp end parallel

    contains
      double precision function call_magnetic(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np

        double precision, dimension(Np) :: p_mag
        p_mag = (/p(1), mag_SLDcube, mag_SLDmatrix, phi, theta/)
        call magnetic_scattering(&
          q, xi, sin2alpha, plus_or_minus, p, p_mag, p_cube, p_cube,&
          Nq, Np, Np, call_magnetic)
      end function call_magnetic
  end subroutine magnetic_oriented_formfactor

  subroutine magnetic_formfactor(&
    q, a, SLDcube, SLDmatrix, siga, mag_SLDcube, mag_SLDmatrix,&
    xi, sin2alpha, plus_or_minus,&
    x_herm, w_herm, x_leg, w_leg, N_gh_order, N_gl_order, &
    Nq, ff_out)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: a, SLDcube, SLDmatrix
    double precision, intent(in) :: siga
    double precision, intent(in) :: mag_SLDcube, mag_SLDmatrix
    double precision, intent(in) :: xi, sin2alpha, plus_or_minus
    double precision, dimension(N_gh_order), intent(in) :: x_herm, w_herm
    double precision, dimension(N_gl_order), intent(in) :: x_leg, w_leg
    integer, intent(in) :: Nq, N_gh_order, N_gl_order
    double precision, dimension(Nq), intent(out) :: ff_out

    integer, parameter :: Np=5
    double precision, dimension(Np) :: p

    double precision :: current_a, current_phi, current_th, current_sinth
    double precision :: w_a_i, w_phi_i, w_theta_i

    double precision :: intsum_norm
    double precision, dimension(Nq) :: hff, intsum_ff

    integer :: ia, ith, iphi
    double precision :: integrand_prefacs

    intsum_ff = 0d0
    intsum_norm = 0d0
    if (siga > tolerance) then
      do ia=1, N_gh_order
        current_a = a*dexp(sq2 * x_herm(ia) * siga)
        w_a_i = w_herm(ia)
        do ith=1, N_gl_order
          current_th = pi/4d0 * (x_leg(ith) + 1d0) ! -1 to 1 for 0 to pi/2
          current_sinth = sin(current_th)
          w_theta_i = pi/4d0 * w_leg(ith)

          do iphi=1, N_gl_order
            current_phi = pi/4d0 * (x_leg(iphi) + 1d0) ! -1 to 1 for 0 to pi/2
            w_phi_i =  pi/4d0 * w_leg(iphi)
            call magnetic_oriented_formfactor(q, current_a,&
              SLDcube, SLDmatrix, current_phi, current_th,&
              mag_SLDcube, mag_SLDmatrix,&
              xi, sin2alpha, plus_or_minus, Nq, hff)

            integrand_prefacs = w_a_i * w_phi_i * w_theta_i * current_sinth

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
          call magnetic_oriented_formfactor(q, a,&
            SLDcube, SLDmatrix, current_phi, current_th,&
            mag_SLDcube, mag_SLDmatrix,&
            xi, sin2alpha, plus_or_minus, Nq, hff)

          integrand_prefacs = w_phi_i * w_theta_i * current_sinth

          intsum_ff = intsum_ff + integrand_prefacs * hff
          intsum_norm = intsum_norm + integrand_prefacs
        end do
      end do
    end if
    ff_out = intsum_ff/intsum_norm
  end subroutine magnetic_formfactor

  subroutine sld(a, SLDcube, SLDmatrix, x, sld_array)
    double precision, intent(in) :: a, SLDcube, SLDmatrix

    integer, parameter :: Nx=4
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array

    x(1) = 0
    x(2) = a
    x(3) = a
    x(4) = a + 100

    sld_array(1) = SLDcube
    sld_array(2) = SLDcube
    sld_array(3) = SLDmatrix
    sld_array(4) = SLDmatrix
  end subroutine sld
end module cube
