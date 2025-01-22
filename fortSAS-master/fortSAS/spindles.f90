module spindle
use math
implicit none
contains
    subroutine amplitude_spindle(q, p, Np, p_spindle, r_long_axis, sin_angle)
        ! calculates form factor amplitude
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        double precision, intent(out) :: p_spindle
        double precision, dimension(3), intent(out) :: r_long_axis
        double precision, intent(out) :: sin_angle
        
        double precision :: RL, R0, phi, theta, psi, beta, SLDspindle, SLDmatrix
        double precision :: q_xy_angle
        
        double precision :: qx, qy, acos_arg
        double precision :: V, ratio, alpha, r, qr, K
        RL = p(1)
        R0 = p(2)
        phi = p(3)
        theta = p(4)
        beta = p(5)
        psi = p(6)
        SLDspindle = p(7)
        SLDmatrix = p(8)
        q_xy_angle = p(9)
        
        ratio =  RL/R0
        if (q /= 0d0) then
            ! get orientation of long axis for given phi, theta, psi, beta angles:
            call calc_vec_rotated_long(phi, theta, psi, beta, r_long_axis, sin_angle)

            ! from given q and angle value get position of pixel on detector: 
            qx = q * cos(q_xy_angle)
            qy = q * sin(q_xy_angle)
            
            ! what is angle between q vector and long axis? 
            ! acos results in NaN if argument is numerically slightly larger than 1
            acos_arg = (qx*r_long_axis(1) + qy*r_long_axis(2))/q
            if (acos_arg > 1d0) then
                alpha = 0.
            else if (acos_arg <-1d0) then
                alpha = pi
            else
                alpha = acos(acos_arg) ! numerically unstable
            end if
            
            ! effective radius for ellipsoidal model:
            ! see for example: http://gisaxs.com/index.php/Form_Factor:Ellipsoid_of_revolution
            r = R0*sqrt( 1d0 + (ratio**2-1d0)*cos(alpha)**2 )
            
            qr = q*r
            K = 3* (sin(qr) - qr*cos(qr)) / qr**3
            V = 4d0/3d0*pi*R0**2*RL
            p_spindle = V * (SLDspindle - SLDmatrix) * K 
            ! volume gets squared for formfactor here (opposed to cpp script)
        else
            p_spindle = V * (SLDspindle - SLDmatrix)
        endif
    end subroutine amplitude_spindle
    
    double precision function p_spindle(q, p, Np)
        ! function for purpose of obtaining amplitude in default way
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np

        double precision, dimension(3) :: r_long_axis
        double precision :: sin_angle

        call amplitude_spindle(q, p, Np, p_spindle, r_long_axis, sin_angle)
    end function p_spindle

    double precision function ff_spindle(q, p, Np)
        ! calculates form factor intensity
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        double precision :: p_spindle
        double precision, dimension(3) :: r_long_axis
        double precision :: sin_angle, sin_psi, psi
        
        call amplitude_spindle(q, p, Np, p_spindle, r_long_axis, sin_angle)
        ! get angle between z axis and x, y plane, for integration
        ! norm of r_long_axis is 1 by definition
        ! sin_angle =&
        !         sqrt(r_long_axis(1)**2 + r_long_axis(2)**2) !/ norm2(r_long_axis)
        ff_spindle = abs(p_spindle)**2
    end function ff_spindle
    
    double precision function get_sin_angle(q, p, Np)
        ! for the calculation of angle space for scaling
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        double precision :: phi, theta, psi, beta
        double precision, dimension(3) :: r_long_axis

        phi = p(3)
        theta = p(4)
        beta = p(5)
        psi = p(6)

        call calc_vec_rotated_long(phi, theta, psi, beta, r_long_axis, get_sin_angle)
        
        ! norm of r_long_axis is 1 by definition
        ! get_sin_angle = &
                !sqrt(r_long_axis(1)**2 + r_long_axis(2)**2) !/ norm2(r_long_axis)
    end function get_sin_angle

    double precision function ff_spindle_orientation_average(q, p, Np)
        double precision, intent(in) :: q
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        
        call integrate_size_distribution(q, p, Np, &
                            6, 0d0, pi, 1d0, &
                            ff_spindle, one, &
                            ff_spindle_orientation_average)

    end function ff_spindle_orientation_average
    
    
    subroutine calc_vec_rotated_long(phi, theta, psi, beta, r_long_axis, sin_angle)
        ! From phi, theta, psi, beta calculate r_long_axis
        double precision, intent(in) :: phi, theta, psi, beta
        double precision, dimension(3), intent(out) :: r_long_axis
        double precision, intent(out) :: sin_angle
        
        double precision :: cos_theta, sin_theta
        double precision :: cos_phi, sin_phi
        double precision :: cos_beta, sin_beta
        double precision :: cos_psi, sin_psi

        double precision, dimension(3) :: r_long, r_easy
        double precision, dimension(3) :: r_easy_cross_r_long
        double precision, dimension(3) :: r_long_rotated

        cos_phi = cos(phi)
        sin_phi = sin(phi)
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        cos_psi = cos(psi)
        sin_psi = sin(psi)
        cos_beta = cos(beta)
        sin_beta = sin(beta)

        ! r_easy: B field along (1,0,0), orientation of easy axis
        !                                precessing in a cone around B
        !r_easy = (/cos_theta,&
        !           sin_theta*cos_phi,&
        !           sin_theta*sin_phi/)
        
        ! r_long: long axis of r_easy, vector lies in plane perpendicular to r_easy
        r_long = (/-sin_theta,&
                    cos_theta*cos_phi,&
                    cos_theta*sin_phi/)

        ! to obtain all possible vectors of r_long, integrate psi
        ! use for this purpose Rodrigues rotation formula, need cross product:
        r_easy_cross_r_long = (/0d0,&
                               -sin_phi,&
                                cos_phi/)
        ! rotation formula:
        r_long_rotated = r_long*cos_psi + r_easy_cross_r_long*sin_psi
        !last part rotation formula not necessary as r_easy, r_long are perpendicular by definition
        !see: r_easy*dot_product(r_easy, r_long)*(1d0-cos(psi)) 

        sin_angle  = sqrt(r_long_rotated(2)**2 + r_long_rotated(3)**2)
        
        !sin_angle = sin_psi
        !print *, "sin_angle", sin_angle
        !print *, "sin_psi", sin_psi
        ! rotate frame of reference by beta around y, to switch magnetic field direction
        r_long_axis = (/cos_beta*r_long_rotated(1) - sin_beta*r_long_rotated(3),&
                       r_long_rotated(2),&
                       sin_beta*r_long_rotated(1) + cos_beta*r_long_rotated(3)/)
    end subroutine calc_vec_rotated_long
    

    
    subroutine formfactor(q, RL, R0, beta, SLDspindle, SLDmatrix, &
                        phi_arr, p_phi, q_xy_angle,&
                        sigRL, theta_arr, p_theta, Ntheta, &
                        Nq, Nphi, ff_intensity)
        double precision, dimension(Nq), intent(in) :: q
        double precision, intent(in) :: RL, R0, beta,  SLDspindle, SLDmatrix, q_xy_angle
        double precision, dimension(Nphi), intent(in) :: phi_arr, p_phi
        double precision, dimension(Ntheta), intent(in) :: theta_arr, p_theta
        double precision, intent(in) :: sigRL
        integer, intent(in) :: Nq, Nphi, Ntheta
        
        integer, parameter :: Np=9
        double precision, dimension(Np) :: p, parallel_p
        double precision, dimension(Nq), intent(out) :: ff_intensity
        
        integer :: iq, iphi
        integer :: counter
        double precision :: beta_deg 
        double precision :: qval
        double precision :: angle_integration
        double precision :: phimin, phimax, phi_probability, RLmin, RLmax
        double precision :: theta, psi
        double precision :: angle_integration_sum, ff_phi_sum, p_phi_sum
        double precision :: dphi, dtheta
        double precision :: ff_help_trapezoid, ff_help_trapezoid_sum
        double precision :: ff_outer, ff_outer_sum, Fq
        double precision :: integrand_i, integrand_i_plus_1
        integer :: itheta
        double precision :: p_theta_sum, phi_diff
        
        ! user enters angles in degree: translate to rad:
        beta_deg = beta*deg
        
        dphi = phi_arr(2)-phi_arr(1)
        dtheta = theta_arr(2)-theta_arr(1)
        
        p = (/RL, R0, 0d0, 0d0, beta_deg, 0d0,  SLDspindle, SLDmatrix, q_xy_angle*deg/)
        counter = 0
        
        call get_cutoff_lognormal(R0, sigRL, RLmin, RLmax)  
        
        !$omp parallel private(qval, iphi, itheta) &
        !$omp& private(parallel_p, ff_help_trapezoid, ff_help_trapezoid_sum, Fq) &
        !$omp& private(integrand_i, integrand_i_plus_1) &
        !$omp& private(ff_outer, ff_outer_sum, p_theta_sum) &
        !$omp& private(p_phi_sum,phi_diff) shared(counter)
         !$omp do          
        do iq=1, Nq
            parallel_p = p
            qval = q(iq)
            ff_help_trapezoid = 0d0
            p_theta_sum = 0d0
            integrand_i = 0d0
            integrand_i_plus_1 = 0d0
            phi_diff = 0d0
            ff_help_trapezoid_sum = 0d0
            ff_outer = 0d0
            ff_outer_sum = 0d0
            
            do itheta=1, Ntheta
            ff_help_trapezoid_sum = 0d0 
            p_phi_sum = 0d0
            parallel_p(3) = phi_arr(1)
            parallel_p(4) = theta_arr(itheta) 
                !!TRAPEZOIDAL
            call integrate_size_distribution(qval, parallel_p, Np, &
                        2, RLmin, RLmax, sigRL, &
                        ff_spindle_orientation_average, lognormal,&
                        ff_help_trapezoid) 
            integrand_i = ff_help_trapezoid * p_phi(1)
            
            do iphi=1, Nphi-1
                parallel_p(3) = phi_arr(iphi+1)
                call integrate_size_distribution(qval, parallel_p, Np, &
                        2, RLmin, RLmax, sigRL, &
                        ff_spindle_orientation_average, lognormal,&
                        ff_help_trapezoid)
                
                integrand_i_plus_1 = ff_help_trapezoid * p_phi(iphi+1)
                
                phi_diff = phi_arr(iphi+1) - phi_arr(iphi)
                
                ff_help_trapezoid_sum = ff_help_trapezoid_sum +&
                phi_diff*(integrand_i_plus_1 + integrand_i)/2
                
                p_phi_sum = p_phi_sum +&
                phi_diff * (p_phi(iphi+1) + p_phi(iphi))/2
                
                integrand_i = integrand_i_plus_1
                
            end do
            ff_outer = (ff_help_trapezoid_sum/p_phi_sum) * p_theta(itheta)
            ff_outer_sum = ff_outer_sum + ff_outer
            p_theta_sum = p_theta_sum + p_theta(itheta)
            
            end do
            Fq = ff_outer_sum *dtheta /p_theta_sum    
            ff_intensity(iq) = Fq
            counter = counter + 1
        write(*, '(1a1,I3, A2)', advance="no") char(13), int(100d0*counter/(Nq)), ' %'            
        end do
         !$omp end do
        !$omp end parallel
        
        ! update progress in line:

    
        !! divide by integration over angle space
        !do iphi=1, Nphi
        !    p(3) = phi_arr(iphi)
        !    call integrate_two_size_distributions(0d0, p, Np, &
        !        4, 0d0, 2d0*pi, 1d0, &
        !        5, 0d0, pi, 1d0, &
        !        get_sin_angle, one, one,&
        !        angle_integration)
        !    angle_integration_sum = angle_integration_sum + angle_integration
        !    
        !    
        !end do
        !
        !ff_intensity = ff_intensity / angle_integration_sum 
        
        
    end subroutine formfactor

    subroutine formfactor_azimuth(angles, q, RL, R0, beta, SLDspindle, SLDmatrix, &
                        phi_arr, p_phi, sigRL, theta_arr, p_theta, Ntheta, Nangles, Nphi, &
                        ff_intensity)

        double precision, intent(in), dimension(Nangles) :: angles ! in deg
        double precision, intent(in) :: q
        double precision, intent(in) :: RL, R0, beta,  SLDspindle, SLDmatrix
        double precision, intent(in) :: sigRL
        integer, intent(in) :: Nangles, Nphi, Ntheta
        double precision, dimension(Nphi), intent(in) :: phi_arr, p_phi
        double precision, dimension(Ntheta), intent(in) :: theta_arr, p_theta
        
        integer, parameter :: Np=9
        double precision, dimension(Np) :: p, parallel_p
        double precision, dimension(Nangles), intent(out) :: ff_intensity
        
        integer :: iangle, iphi
        integer :: counter
        double precision :: phi_deg, beta_deg, sigphi_deg
        double precision :: qval
        double precision :: phimin, phimax, phi_probability, RLmin, RLmax
        double precision :: theta, phi
        double precision :: angle_integration_sum, angle_integration
        double precision :: dphi, p_phi_sum, dtheta
        double precision :: ff_help_trapezoid, ff_help_trapezoid_sum
        double precision :: ff_outer, ff_outer_sum, Fq
        double precision :: integrand_i, integrand_i_plus_1
        integer :: itheta
        double precision :: p_theta_sum, phi_diff
        
        ! user enters angles in degree: translate to rad:
        beta_deg = beta*deg
        
        dphi = phi_arr(2)-phi_arr(1)
        dtheta = theta_arr(2)-theta_arr(1)

        p = (/RL, R0, 0d0, 0d0, beta_deg, 0d0,  SLDspindle, SLDmatrix, 0d0/)
        counter = 0
        
        call get_cutoff_lognormal(R0, sigRL, RLmin, RLmax)  

        !$omp parallel private(qval, iphi, itheta) &
        !$omp& private(parallel_p, ff_help_trapezoid, ff_help_trapezoid_sum, Fq) &
        !$omp& private(integrand_i, integrand_i_plus_1) &
        !$omp& private(ff_outer, ff_outer_sum, p_theta_sum) &
        !$omp& private(p_phi_sum,phi_diff) shared(counter)
         !$omp do          
        do iangle=1, Nangles
            parallel_p = p
            parallel_p(9) = angles(iangle) * deg
            ff_help_trapezoid = 0d0
            p_theta_sum = 0d0
            integrand_i = 0d0
            integrand_i_plus_1 = 0d0
            phi_diff = 0d0
            ff_help_trapezoid_sum = 0d0
            ff_outer = 0d0
            ff_outer_sum = 0d0
            
            do itheta=1, Ntheta
            ff_help_trapezoid_sum = 0d0 
            p_phi_sum = 0d0
            parallel_p(3) = phi_arr(1)
            parallel_p(4) = theta_arr(itheta) 
                !!TRAPEZOIDAL
            call integrate_size_distribution(q, parallel_p, Np, &
                        2, RLmin, RLmax, sigRL, &
                        ff_spindle_orientation_average, lognormal,&
                        ff_help_trapezoid) 
            integrand_i = ff_help_trapezoid * p_phi(1)
            
            do iphi=1, Nphi-1
                parallel_p(3) = phi_arr(iphi+1)
                call integrate_size_distribution(q, parallel_p, Np, &
                        2, RLmin, RLmax, sigRL, &
                        ff_spindle_orientation_average, lognormal,&
                        ff_help_trapezoid)
                
                integrand_i_plus_1 = ff_help_trapezoid * p_phi(iphi+1)
                
                phi_diff = phi_arr(iphi+1) - phi_arr(iphi)
                
                ff_help_trapezoid_sum = ff_help_trapezoid_sum +&
                phi_diff*(integrand_i_plus_1 + integrand_i)/2
                
                p_phi_sum = p_phi_sum +&
                phi_diff * (p_phi(iphi+1) + p_phi(iphi))/2
                
                integrand_i = integrand_i_plus_1
                
            end do
            ff_outer = (ff_help_trapezoid_sum/p_phi_sum) * p_theta(itheta)
            ff_outer_sum = ff_outer_sum + ff_outer
            p_theta_sum = p_theta_sum + p_theta(itheta)
            
            end do
            Fq = ff_outer_sum *dtheta /p_theta_sum    
            ff_intensity(iangle) = Fq
            counter = counter + 1
            ! update progress in line:
            write(*, '(1a1,I3, A2)', advance="no") char(13), int(100d0*counter/(Nangles)), ' %'
            
        end do
         !$omp end do
        !$omp end parallel
        
        
    
        
        
        
        
        
        !! divide by integration over angle space
        !do iphi=1, Nphi
        !    p(3) = phi_arr(iphi)
        !    call integrate_two_size_distributions(0d0, p, Np, &
        !                4, 0d0, 2d0*pi, 1d0, &
        !                5, 0d0, pi, 1d0, &
        !                get_sin_angle, one, one,&
        !                angle_integration)
        !    angle_integration_sum = angle_integration_sum + angle_integration
        !    
        !    
        !end do
        !ff_intensity = ff_intensity / angle_integration_sum 
    end subroutine formfactor_azimuth
    
    
    
    subroutine get_2dimage_precalculated_theta_distribution(&
                    qy, qz, RL, R0, beta, SLDspindle, SLDmatrix, &
                    phi_arr, p_phi, sigRL, theta_arr, p_theta, Nphi, Ntheta, Nqy, Nqz, ff_intensity)
                    
                        !subroutine get_2dimage(qx, qy, RL, R0, phi, beta, SLDspindle, SLDmatrix, &
                        !sigRL, sigphi,&
                        !Nqx, Nqy, ff_intensity)
                    
        double precision, dimension(Nqy), intent(in) :: qy
        double precision, dimension(Nqz), intent(in) :: qz
        double precision, intent(in) :: beta, SLDspindle, SLDmatrix
        double precision, dimension(Nphi), intent(in) :: phi_arr, p_phi
        integer, intent(in) :: Nphi, Ntheta
        double precision, intent(in) :: sigRL
        integer, intent(in) :: Nqy, Nqz
        double precision, dimension(Ntheta), intent(in) :: theta_arr, p_theta
        double precision :: beta_deg, phi_deg, R0, RL
        
        integer, parameter :: Np=9
        integer :: counter
        double precision, dimension(Np) :: p, parallel_p
        double precision, dimension(Nqy, Nqz), intent(out) :: ff_intensity
        
        integer :: iqy, iqz
        double precision :: angle_integration, angle_integration_sum
        double precision :: theta_deg, sigth_deg
        double precision :: qyval, qzval, qval, qy2, qz2
        double precision :: RLmin, RLmax, thmin, thmax, phi_probability
        double precision :: integrand_i, integrand_i_plus_1

        ! define arrays for phi and their probability

        ! set this correctly
        !integer, parameter :: Nphi=199 
        ! set this correctly

        integer :: iphi, itheta
        double precision :: ff_help_trapezoid, ff_help_trapezoid_sum, p_phi_sum
        double precision :: sin_sum, phi_sum, prob, ff_help_ave
        double precision :: ref_ff, phi_diff , Fq, V
        double precision :: p_phi_help, p_phi_help_plus_1
        double precision :: dphi, dtheta
        double precision :: ff_outer, ff_outer_sum, p_theta_sum
        
        ! user enters angles in degree: translate to rad:
        
        beta_deg = beta*deg
        
        dphi = phi_arr(2)-phi_arr(1)
        dtheta = theta_arr(2)-theta_arr(1)
        
        
        call get_cutoff_lognormal(R0, sigRL, RLmin, RLmax)
        !call get_cutoff_gaussian(phi_deg, sigphi_deg, phimin, phimax)

        p = (/RL, R0, 0d0, 0d0, beta_deg, 0d0,  SLDspindle, SLDmatrix, 0d0/)
        
        
        counter = 0
        
        !$omp parallel private(qyval, qzval, qval, qy2, iqy, iqz, iphi, itheta) &
        !$omp& private(parallel_p, ff_help_trapezoid, ff_help_trapezoid_sum, Fq) &
        !$omp& private(integrand_i, integrand_i_plus_1 ) &
        !$omp& private(ff_outer, ff_outer_sum, p_theta_sum) &
        !$omp& private(p_phi_sum,phi_diff) shared(counter)
         !$omp do         
         do iqy=1, Nqy/2
            qyval = qy(iqy)
            qy2 = qyval**2
            parallel_p = p
            do iqz=1, Nqz
                qzval = qz(iqz)
                qval = sqrt(qy2 + qzval**2)
                parallel_p(9) = atan2(qzval, qyval)

                ff_help_trapezoid = 0d0
                p_theta_sum = 0d0
                integrand_i = 0d0
                integrand_i_plus_1 = 0d0
                phi_diff = 0d0
                ff_help_trapezoid_sum = 0d0
                ff_outer = 0d0
                ff_outer_sum = 0d0
                
                
                do itheta=1, Ntheta
                ff_help_trapezoid_sum = 0d0 
                p_phi_sum = 0d0
                parallel_p(3) = phi_arr(1)
                parallel_p(4) = theta_arr(itheta) 
                 !!TRAPEZOIDAL
                call integrate_size_distribution(qval, parallel_p, Np, &
                            2, RLmin, RLmax, sigRL, &
                            ff_spindle_orientation_average, lognormal,&
                            ff_help_trapezoid) 
                integrand_i = ff_help_trapezoid * p_phi(1)
                
                do iphi=1, Nphi-1
                    
                    parallel_p(3) = phi_arr(iphi+1)
                    call integrate_size_distribution(qval, parallel_p, Np, &
                            2, RLmin, RLmax, sigRL, &
                            ff_spindle_orientation_average, lognormal,&
                            ff_help_trapezoid)
                    
                    integrand_i_plus_1 = ff_help_trapezoid * p_phi(iphi+1)
                    
                    phi_diff = phi_arr(iphi+1) - phi_arr(iphi)
                    
                    ff_help_trapezoid_sum = ff_help_trapezoid_sum +&
                    phi_diff*(integrand_i_plus_1 + integrand_i)/2
                    
                    p_phi_sum = p_phi_sum +&
                    phi_diff * (p_phi(iphi+1) + p_phi(iphi))/2
                    
                    integrand_i = integrand_i_plus_1
                    
                    
                end do
                
                ff_outer = (ff_help_trapezoid_sum/p_phi_sum) * p_theta(itheta)
                ff_outer_sum = ff_outer_sum + ff_outer
                p_theta_sum = p_theta_sum + p_theta(itheta)
                
                end do
                
                Fq = ff_outer_sum *dtheta /p_theta_sum    
                ff_intensity(iqy, iqz) = Fq
                ff_intensity(Nqy - iqy, iqz) = Fq
                
            end do
            counter = counter + 1
            ! update progress in line:
            write(*, '(1a1,I3, A2)', advance="no") char(13), int(100d0*counter/(Nqy/2)), ' %'
            
        end do
         !$omp end do
        !$omp end parallel
        
        print *, '' ! make new line for next output
        
        ! divide by integration over angle space
        !do iphi=1, Nphi
        !p(3) = phi_arr(iphi)
        !call integrate_two_size_distributions(0d0, p, Np, &
        !            4, 0d0, 2d0*pi, 1d0, &
        !            5, 0d0, pi, 1d0, &
        !            get_sin_angle, one, one,&
        !            angle_integration)
        !angle_integration_sum = angle_integration_sum + angle_integration
        !
        !end do
        !ff_intensity = ff_intensity / angle_integration
        !print *, '' ! make new line for next output

    end subroutine get_2dimage_precalculated_theta_distribution   
    
    
    
    
end module spindle
