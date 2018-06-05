 function [d, dx, h, kernel, time_nondim, viscosity, delta_sph, slip_cond,...
    scheme, testcase, dt_save,particles, rho_0,gamma,c_0,p_0,Xi,my,alpha,...
    a_wall, int_fluid, int_boundary,int_bound_remove,t_remove, t_damp_fact, G,...
    Height, vel_wall_cond, eqn_of_state, delta] = SPH_input()

    % dimensions
    d = 2; 
    % distance between particles
    dx = 0.05;
    % smoothing length
    h = dx;
    % define Kernel
    kernel = 1; % 1 quintic spline, 2 wendland
    % nondimensional running time (t/t_ref)
    time_nondim = 5;
    % viscosity model
    viscosity = 1; % 1 artificial, 2 laminar
    % viscosity coefficiant value
    alpha = 0.02; 
    % delta-SPH model
    delta_sph = 2; % 0 no delta-SPH, 1 Moltani, 2 Antuono
    % delta-SPH value
    delta = 0.1;
    % boundary condition
    slip_cond = 1; % 1 no slip, 2 free_slip
    % wall velocity condition
    vel_wall_cond = 1; % 1 vel_wall from boundary_update, 2 vel_wall in continuity = 0
    % scheme
    scheme = 1; % 1 velocity verlet, 2 - to be included at later time
    % equation of state
    eqn_of_state = 1; % 1 weakly compressible, 2 ideal gas (not tested)
    % saving time increment
    dt_save = 0.01;
    % testcase 
    global testcase
    testcase = 2; % 1 multifluid, 2 hydrostatic, 3 dambreak, 4 ESA Case, 5 SF Flat interface, 6 Free Bubble,7 Meniscus,8 Drop, 9 Crown
    
    if testcase == 1 % multifluid
        G = 1 * [ 0, -9.81];
        Height = 0.1;
        v_max = sqrt(abs(G(2)) * Height);        
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary] = multifluid(kernel, dx, d, v_max, alpha);
        t_damp_fact = 1;
        t_remove = 0;
        int_bound_remove = [];
        
    elseif testcase == 2 % hydrostatic        
        G = 1 * [ 0, -9.81];
        Height = 1.0;
        v_max = sqrt(abs(G(2)) * Height);
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary] = hydrostat(kernel, dx, d, v_max, alpha);
        % hydrostatic initialisation
%         particles(int_fluid,8) = rho_0(2) .* abs(G(2)) .* (Height - particles(int_fluid,2));
%         particles(int_fluid,5) = rho_0(2) * ( (particles(int_fluid,8) - Xi(2)) / p_0(2) + 1) .^(1/gamma(2));
        % initialisation time parameters
        t_damp_fact = 1;        
        t_remove = 0;
        int_bound_remove = []; 
        
    elseif testcase == 3   % dambreak        
        G = 1 * [ 0, -1.0];
        Height = 1;
        v_max = sqrt(abs(G(2)) * Height); 
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary, int_bound_remove] = dambreak(kernel, dx, d, v_max, alpha);
        t_damp_fact = 0;
        t_remove = 0;
        %initial pressure, hydrostatic
        particles(int_fluid,8) = rho_0(2) .* abs(G(2)) .* (Height - particles(int_fluid,2));
        particles(int_fluid,5) = rho_0(2) * ( (particles(int_fluid,8) - Xi(2)) / p_0(2) + 1) .^(1/gamma(2));
    
    elseif testcase == 4 % ESA Case
        dx = 0.0015;
        G = 1 * [ 0, -1.0];
        Height = 0.03;
        v_max = sqrt(abs(G(2)) * Height);
        % create particle matrix and fluid properties
        [particles, rho_0,gamma,c_0,p_0,Xi,my,alpha, a_wall, int_fluid, int_boundary] = ESA_Case(kernel, dx, d, v_max, alpha);
        t_damp_fact = 1;        
        t_remove = 0;
        int_bound_remove = [];
        
    elseif testcase == 5 % Surface Tension Flat interface
        G = 1 * [ 0, -9.81];
        Height = 0.1;
        v_max = sqrt(abs(G(2)) * Height);        
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary] = SF_Flat_interface(kernel, dx, d, v_max, alpha);
        t_damp_fact = 1;
        t_remove = 0;
        int_bound_remove = [];
        
    elseif testcase == 6 % Surface Tension Free Bubble
        G = 1* [ 0, -1.0];
        Height = 0.4;
        v_max = sqrt(abs(G(2)) * Height);        
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary] = Free_Bubble(kernel, dx, d, v_max, alpha);
        t_damp_fact = 0;
        t_remove = 0;
        int_bound_remove = [];
        
    elseif testcase == 7 % Meniscus      
        G = 1 * [ 0, -9.81];
        Height = 1.0;
        v_max = sqrt(abs(G(2)) * Height);
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary] = Meniscus(kernel, dx, d, v_max, alpha);
        t_damp_fact = 1;        
        t_remove = 0;
        int_bound_remove = [];
        
    elseif testcase == 8 % Drop      
        G = 1 * [ 0, -9.81];
        Height = 0.1;
        v_max = sqrt(abs(G(2)) * Height);
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary] = Drop(kernel, dx, d, v_max, alpha);
        t_damp_fact = 1;        
        t_remove = 0;
        int_bound_remove = []; 
    elseif testcase == 9 % Crown      
        G = 1 * [ 0, -9.81];
        Height = 0.4;
        v_max = sqrt(abs(G(2)) * Height);
        % create particle matrix and fluid properties
        [particles, rho_0, gamma, c_0, p_0, Xi, my, alpha, a_wall, int_fluid, int_boundary] = Crown(kernel, dx, d, v_max, alpha);
        t_damp_fact = 1;        
        t_remove = 0;
        int_bound_remove = []; 
    end
end
