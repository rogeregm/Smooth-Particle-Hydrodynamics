function [particles, int] = initialisation(particles, phase, flag, rho_0, dx, d, vel)
    
    int = size(particles,1) + (1 : size(phase,1));
    particles(int,1:9) = 0;
    
    particles(int,1:2)= phase;                 % coordinates of phase particles
    particles(int, 3) = flag;                  % integer flag, 1 for boundary
    particles(int, 5) = rho_0;                 % density
    particles(int, 4) = particles(int,5) * dx^d;% mass of particles is fixed
    particles(int, 6) = vel(1);                % initial velocity x-component
    particles(int, 7) = vel(2);                %                  y-component
    %particles(int,8) = 0;                     % pressure (will be calculated from density)    
    %particles(int,9) = 0;                     % temperature
    
end