function [p] = pressure(p_0,rho_0,gamma,Xi,particles, int_fluid, eqn_of_state)
    
    p = particles(1:length(int_fluid),8);
      
    for flag = 2:size(p_0,2) % for the different fluid particles
        % make a vector of integers for each fluid phase
        ct = 0;
        int = 0;
        for m = 1:size(particles(:,3));
            if particles(m,3) == flag
                ct = ct + 1;
                int(ct) = m;
            end
        end
        % end making vector
        
        rho = particles(int,5);
        if eqn_of_state == 1
        % pressure for weakly compressible fluid
        p(int) = p_0(flag) * ( (rho/rho_0(flag)).^gamma(flag) - 1 ) + Xi(flag);
        % pressure for ideal gas law
        elseif eqn_of_state == 2 % ideal gas
            R = 287;
            %T = particles(int,9);
            p(int) = rho.^1.4;
        end
        
    end  
end