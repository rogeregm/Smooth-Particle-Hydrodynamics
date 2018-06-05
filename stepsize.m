function [dt] = stepsize(c_0, particles, h, alpha, d, g, viscosity, my)

    % CFL condition (eq 19)
    c_max = max(c_0);
    vel = sqrt(particles(:,6).^2 + particles(:,7).^2); % velocity magnitude of each particle
    v_max = max(vel);
    dt(1) = 0.25 * h / (c_max + (v_max));
    
    % viscous condition (eq 20)
    
    if viscosity == 1  || viscosity == 3 || viscosity == 4 % artificial
        my = 0.5/(d+2) * max(alpha) * h * c_max;
    elseif viscosity == 2       % laminar
        my = max(my);
    end
    dt(2) = 0.125 * h^2 / my;
    
    % body force condition (eq 21)
    dt(3) = 0.25 * sqrt(h/norm(g));

    % global timestep
    dt = min(dt);

end