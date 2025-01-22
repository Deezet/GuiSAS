module ellipsoid_cs
use math
implicit none
contains
    double precision function p_ellipsoid(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
            
        double precision :: R_z, R_r, alpha, SLDellipsoid, SLDmatrix
        double precision :: R_theta
        R_z = p(1)
        R_r = p(2)
        alpha = p(3)
        SLDellipsoid = p(4)
        SLDmatrix = p(5)
            
        R_theta = (R_r**2*(cos(alpha*deg))**2 + R_z**2*(sin(alpha*deg))**2)**0.5
        
        p_ellipsoid = 4d0*pi*R_r**2*R_z*(SLDellipsoid - SLDmatrix)*(sin(q*R_theta)-q*R_theta*cos(q*R_theta))/(q*R_theta)**3
    end function p_ellipsoid

    double precision function p_ellipsoid_cs(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        
        double precision :: R_z, R_r, d_s, alpha, SLDcore, SLDshell, SLDmatrix    !R_z = Radius in z direction, R_r = radius in xy direction, d_s = thickness of the shell
        double precision :: ff_amp_core, ff_amp_shell

        R_z = p(1)
        R_r = p(2)
        d_s = p(3)
        alpha = p(4)
        SLDcore = p(5)
        SLDshell = p(6)
        SLDmatrix = p(7)

        ff_amp_shell = p_ellipsoid(q, (/R_z+d_s, R_r+d_s, alpha, SLDshell, SLDmatrix/), 5)
        ff_amp_core = p_ellipsoid(q, (/R_z, R_r, alpha, SLDcore, SLDshell/), 5)
        p_ellipsoid_cs = ff_amp_shell + ff_amp_core
    end function p_ellipsoid_cs

  double precision function ff_ellipsoid_cs(q, p, Np)
    double precision, intent(in) :: q
    double precision, dimension(Np), intent(in) :: p
    integer, intent(in) :: Np
        
    ff_ellipsoid_cs = abs(p_ellipsoid_cs(q, p, Np))**2
  end function ff_ellipsoid_cs

  subroutine formfactor(q, R_z, R_r, d_s, alpha, SLDcore, SLDshell, SLDmatrix, &
                    sigR_r, sigd_s, sigalpha,&
                    Nq, ff_intensity)
    double precision, dimension(Nq), intent(in) :: q
    double precision, intent(in) :: R_z, R_r, d_s, alpha, SLDshell, SLDcore, SLDmatrix
    double precision, intent(in) :: sigR_r, sigd_s, sigalpha
    integer, intent(in) :: Nq
        
    integer, parameter :: Np=7
    double precision, dimension(Np) :: p
    double precision, dimension(Nq), intent(out) :: ff_intensity
    integer :: iq
    double precision :: alphamin, alphamax, R_rmin, R_rmax, d_smin, d_smax
    double precision :: alpha_probability
    p = (/R_z, R_r, d_s, alpha, SLDcore, SLDshell, SLDmatrix/)
    call get_cutoff_lognormal(R_r, sigR_r, R_rmin, R_rmax) 
    call get_cutoff_lognormal(d_s, sigd_s, d_smin, d_smax) 
    call get_cutoff_gaussian(alpha, sigalpha, alphamin, alphamax)

    !$omp parallel
    !$omp do
    do iq=1, Nq
        call integrate_three_size_distributions(q(iq), p, Np, &
                        2, R_rmin, R_rmax, sigR_r, &
                        3, d_smin, d_smax,sigd_s, &
                        4, alphamin, alphamax, sigalpha, &
                        ff_ellipsoid_cs, lognormal, lognormal, sindeg_gaussian,&
                        ff_intensity(iq))
    end do
    !$omp end do
    !$omp end parallel
        
    call integrate_prob_distribution(alpha, sigalpha, alphamin, alphamax,&
                    sindeg_gaussian, alpha_probability)
    ff_intensity = ff_intensity / alpha_probability
  end subroutine formfactor
    

  subroutine sld(R_z, d_s, SLDcore, SLDshell, SLDmatrix, x, sld_array)
    !R_z, R_r, d_s, alpha, SLDcore, SLDshell, SLDmatrix
    double precision, intent(in) :: R_z, d_s, SLDcore, SLDshell, SLDmatrix
   
    integer, parameter :: Nx=6
    double precision, dimension(Nx), intent(out) :: x
    double precision, dimension(Nx), intent(out) :: sld_array
    
      x(1) = 0
      x(2) = R_z
      x(3) = R_z
      x(4) = R_z + d_s
      x(5) = R_z + d_s
      x(6) = R_z + d_s + 100
    
      sld_array(1) = SLDcore
      sld_array(2) = SLDcore
      sld_array(3) = SLDshell
      sld_array(4) = SLDshell
      sld_array(5) = SLDmatrix
      sld_array(6) = SLDmatrix
  end subroutine sld
end module ellipsoid_cs
