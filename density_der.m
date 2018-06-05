function [drho_by_dt, dissipation] = density_der(particles, int_fluid, h, kernel, d, rho_0, p_0, Xi, gamma, dx, vel_wall_cond, eqn_of_state, delta_sph)

% initialisation
drho_by_dt(int_fluid) = 0;
dissipation(int_fluid) = 0;
%Antuono terms initialisation
der_rho_L(d,int_fluid) = 0;
diss_L(d,d,int_fluid) = 0;

% cutoff radius
if kernel == 1
    r_c = 3 * h;
elseif kernel == 2
    r_c = 2 * h;
end

% range-search
Y = particles(int_fluid,1:2);
NS = particles(:,1:2);
[idx] = rangesearch(NS,Y,r_c);

%Antuono delta-SPH
if (delta_sph == 2)
    [diss_L] = calc_corr(particles, idx, int_fluid, h, kernel, d);
    [der_rho_L] = calc_der_rho_corr(particles, idx, int_fluid, h, kernel, d, diss_L);
end

parfor a = int_fluid % loop over all fluid particles
    
    % calculate interaction of one fluid particle with all other particles
    for b = idx{a}(2:end)

        flag_b = particles(b,3);
        flag_a = particles(a,3);
        
        %distance between particles
        drx = particles(a,1) - particles(b,1);
        dry = particles(a,2) - particles(b,2);
        rad = sqrt(drx^2 + dry^2);
        
        % nondimensional distance
        q = rad / h;
        %kernel and derivative values
%         Wab = kernel_fct(kernel, d, h, q);  %only for mixed correction (unsused)
        DWab = kernel_der(kernel, d, h, q) / h;
        
        % kernel derivative with respect to x_a
        Fab = [(drx/rad);(dry/rad)]*DWab;
        
        % densities, mass
        rho_a = particles(a,5);
        if particles(b,3) == 1 %if b boundary particle
            % calculate density of boundary particles, depending on interacting fluid particle
            rho_b = boundary_rho(particles, a, rho_0, p_0, Xi, gamma,eqn_of_state);
            m_b = rho_0(flag_a) * dx*dx;
            % velocity difference between particles
            if vel_wall_cond == 1
                dvx = particles(a,6) - particles(b,6);
                dvy = particles(a,7) - particles(b,7);
            elseif vel_wall_cond == 2
                dvx = particles(a,6) - 0;
                dvy = particles(a,7) - 0;
            end
        else % straightforward if fluid particle
            rho_b = particles(b,5);
            m_b = particles(b,4);
            % velocity difference between particles
            dvx = particles(a,6) - particles(b,6);
            dvy = particles(a,7) - particles(b,7);
        end
        
                
        % continuity equation (eq6)        
        drho_by_dt(a) = drho_by_dt(a) + rho_a * m_b / rho_b * [dvx;dvy]'*Fab;
                       
        % delta-SPH model (general)
        if delta_sph ~= 0 && flag_b ~= 1            
            if (delta_sph == 2) %Antuono                
                psi_ab = (rho_b-rho_a) - 0.5 * (der_rho_L(:,b) + der_rho_L(:,a))' * [-drx;-dry];
                Rab = [-drx;-dry]/(rad^2);
            elseif (delta_sph == 1) % Moltani                
                psi_ab = (rho_b-rho_a);
                Rab = [-drx;-dry]/(rad^2);
            end
            
            %Dissipative term            
            dissipation(a) = dissipation(a) + 2 * (m_b/rho_b * psi_ab * Rab'*Fab);
        end        
    end
end
