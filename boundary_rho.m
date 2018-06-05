function [density] = boundary_rho(particles, int, rho_0, p_0, Xi, gamma,eqn_of_state)

    % flag depending on fluid phase
    flag = particles(int,3);
    
    if eqn_of_state == 1
        % calculate density from pressure for wall particles
        density = rho_0(flag) * ( (particles(int,8) - Xi(flag)) / p_0(flag) + 1) .^(1/gamma(flag));
    elseif eqn_of_state == 2
        % ideal gas law
        R = 287;
        density = particles(int,8).^1/1.4;
    end
    
end