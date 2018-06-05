function [dvelx_by_dt, dvely_by_dt, diss_x, diss_y,PF_sum,curv_a] = velocity_der(particles, int_fluid, h,...
    kernel, d, alpha, c_0, g, viscosity, slip_cond, rho_0, p_0, Xi, gamma,...
    my, eqn_of_state, int_boundary)

    % initialisation
    dvelx_by_dt(int_fluid) = 0;
    dvely_by_dt(int_fluid) = 0;
    dx = h;
    PF_sum(int_fluid) = 0;
    epsilon = 0.01; % parameter to avoid zero denominator   
    diss_x(int_fluid) = 0;
    diss_y(int_fluid) = 0;
    global testcase
    % cutoff radius
    if kernel == 1
        r_c = 3 * h;
    elseif kernel == 2
        r_c = 2 * h;
    end   
     
    % range-search
    Y = particles(int_fluid,1:2);   %List of fluid particle positions
    NS = particles(:,1:2);          %List of all particle positions
    idx = rangesearch(NS,Y,r_c);    %List of all particles within range of each fluid particle
    idx_all = rangesearch(NS,NS,r_c);
    
    parfor a = int_fluid % loop over all fluid particles
        
        % calculate interaction of one fluid particle with all particles within range r_c
        for b = idx{a}(2:end)

            %distance between particles
            drx = particles(a,1) - particles(b,1);
            dry = particles(a,2) - particles(b,2);
            rad = sqrt(drx^2 + dry^2);
               
            % nondimensional distance
            q = rad / h;
            % derivative of Kernel
            der_W = kernel_der(kernel, d, h, q) / h;
                        
            % velocity difference between particles
            dvx = particles(a,6) - particles(b,6);
            dvy = particles(a,7) - particles(b,7);
            
            % pressure of particles
            p_a = particles(a,8);
            p_b = particles(b,8);
            
            % particle type
            flag_a = particles(a,3);
            flag_b = particles(b,3);
            
            % densities, mass
            rho_a = particles(a,5);
            m_a   = particles(a,4);
            if flag_b == 1 % if b is a boundary particle
                % calculate density
                rho_b = boundary_rho(particles, a, rho_0, p_0, Xi, gamma,eqn_of_state);
                m_b = rho_0(flag_a) * dx*dx;
                  
            else % if b is a fluid particle
                rho_b = particles(b,5);
                m_b = particles(b,4);
            end
            rho_ab = 0.5 * (rho_a + rho_b);     
            
% RG        % Surface Tension and Virial Pressure
            %   1 = Cosine Pairwise Force
            %   3 = PF 3
            %   4 = Kernel Weighted Attraction
            %   5 = Cohesion Force
            [PF, Virial_Pressure,adhesion_F, der_color_field_i] = PairwiseForce(5, particles, h, rad, a, b, idx, rho_0, p_0, r_c, rho_b, m_b, Xi, gamma, eqn_of_state, int_boundary,idx_all);
            PF_sum(a) =  PF_sum(a) + PF;
            curv_a(a) = der_color_field_i;
            if der_color_field_i>0
                VP=0;
            else
                VP=1;
            end
            virial_radius = max(particles(:,1))-3*h;
            % momentum equation (eq7) Pressure Gradient part
            p_ab = (rho_b * p_a + rho_a * p_b) / (rho_a + rho_b);   %(eq8)
            pressure_fact = - 1/m_a * ((m_a/rho_a)^2 + (m_b/rho_b)^2) * (p_ab) * der_W + (1/(4*(2*virial_radius)^2))*PF*rad*VP;  %(eq7)
            
            % Monaghan formulation
            % pressure_fact = - m_b * (p_a / rho_a^2 + p_b / rho_b^2) * der_W ;
            
            % acceleration due to pressure gradient (eq7) + PF
            dvelx_by_dt(a) = dvelx_by_dt(a) + pressure_fact * drx / rad + 1/m_a * PF * drx/rad;
            dvely_by_dt(a) = dvely_by_dt(a) + pressure_fact * dry / rad + 1/m_a * PF * dry/rad;  
                                   
            % if free slip condition aplies only consider particles which are not fluid particles                
            if slip_cond == 1 || flag_b ~= 1
                
                % artificial viscosity (eq11)
                if viscosity == 1 || viscosity == 3 || viscosity == 4
                    
                    if flag_b ~= 1 % so if particle b is a fluid
                        alpha_ab = 0.5 * (alpha(flag_a) + alpha(flag_b));
                        c_ab = 0.5 * (c_0(flag_a) + c_0(flag_b));
                    else
                        alpha_ab = alpha(flag_a);
                        c_ab = c_0(flag_a);
                    end
                    
                    if (drx * dvx + dry * dvy) < 0
                        visc_art_fact = m_b * alpha_ab * h * c_ab * (dvx * drx + dvy * dry)...
                            /(rho_ab * (rad^2 + epsilon * h^2)) * der_W;
                    else
                        visc_art_fact = 0;               
                    end                   
                    dvelx_by_dt(a) = dvelx_by_dt(a) + visc_art_fact * drx / rad;
                    dvely_by_dt(a) = dvely_by_dt(a) + visc_art_fact * dry / rad;
                                        
                % laminar viscosity (eq 10)              
                elseif viscosity == 2
                    
                    my_a = my(flag_a);                    
                    if flag_b ~= 1
                        my_b = my(flag_b);
                    else
                        my_b = my(flag_a);                        
                    end
                    my_ab = 2 * my_a * my_b / (my_a + my_b); % (eq 9)
                    visc_lam_fact = 1/m_a * my_ab * ((m_a/rho_a)^2 + (m_b/rho_b)^2) * der_W;
                    
                    dvelx_by_dt(a) = dvelx_by_dt(a) + visc_lam_fact * dvx / rad;
                    dvely_by_dt(a) = dvely_by_dt(a) + visc_lam_fact * dvy / rad;
                end   
            end
        end  
        % Gravity
        if testcase ~=6
            dvelx_by_dt(a) = dvelx_by_dt(a) + g(1);
            dvely_by_dt(a) = dvely_by_dt(a) + g(2);
        end
    end
end