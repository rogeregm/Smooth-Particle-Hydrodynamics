function [pressure,vel_x,vel_y] = boundary_update(particles, int_fluid, int_boundary, kernel, d, h, g, slip_cond, PF_sum)

    if slip_cond == 1
        % prescribed wall velocity
        v_wall(int_boundary,1:2) = particles(int_boundary,(4:5));        
    end
    % velocity components
    vel_x = particles(:,6);
    vel_y = particles(:,7);
    
    % initialisation pressure
    pressure(1:length(int_boundary)) = 0;
    
    % cutoff radius
    if kernel == 1
        r_c = 3 * h;
    elseif kernel == 2
        r_c = 2 * h;
    end
    
    % range-search
    Y = particles(int_boundary,1:2);
    NS = particles(int_fluid,1:2);
    [idx] = rangesearch(NS,Y,r_c);

%     boundary_samp_dens = Boundary_Sampling_density(particles, d, h, r_c, int_boundary);
    
    for n = 1:length(int_boundary)% over all boundary particles
        w = int_boundary(n); 
        
        % initialisation wall pressure
        sum_pW = 0;
        sum_rhorW = 0;
        sum_W = 0;
        sum_PFW = 0;
        % initialisation wall velocity
        sum_vWx = 0;
        sum_vWy = 0;        
        
        for m = 1:length(idx{n})% interaction with all fluid particles in viscinity
            f = idx{n}(m); 

            %distance between particles
            drx = particles(w,1) - particles(f,1);
            dry = particles(w,2) - particles(f,2);
            rad = sqrt(drx^2 + dry^2);
            
            % nondimensional distance
            q = rad / h;
            % Kernel function
            W = kernel_fct(kernel, d, h, q);

            % build up pressure of wall particle
            sum_pW = sum_pW + particles(f,8) * W;
            sum_rhorW = sum_rhorW + particles(f,5) * W *(drx * g(1) + dry * g(2)); 
            sum_W = sum_W + W;
%             sum_PFW = sum_PFW + particles(f,5) * PF_sum(f) * W * (drx + dry);
            
            if slip_cond == 1
                % building up the SPH average for wall velocity
                sum_vWx = sum_vWx + vel_x(f) * W;
                sum_vWy = sum_vWy + vel_y(f) * W;
            end
        end

        if sum_W == 0
            pressure(n) = 0;
        else
            % combining terms to get pressure of wall particle (eq 27)
            pressure(n) = ( sum_pW + sum_rhorW -sum_PFW) / sum_W ; 
            if slip_cond == 1
                % calculate wall velocity (eq 22 & 23) 
                %v_wall(1:2) = particles(w,(4:5));
                vel_x(w) = 2 * v_wall(w,1) - sum_vWx / sum_W;
                vel_y(w) = 2 * v_wall(w,2) - sum_vWy / sum_W;
            end
        end        
    end
end