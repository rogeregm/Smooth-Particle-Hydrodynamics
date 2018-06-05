function [particles, rho_0,gamma,c_0,p_0,Xi,my,alpha, a_wall, int_fluid, int_boundary, int_bound_remove] = dambreak(kernel, dx, d, v_max, alpha)

    % rectangular domain
    % first fluid phase
    % specify coordinates of edges for fluid
    f_lowleft =  [-0.0  0.0];
    f_lowright = [ 2.0  0.0];
    f_upleft =   [-0.0  1.0];
    f_upright =  [ 2.0  1.0];
    % properties of fluid (will be stored)
    flag = 2;
    rho_0(flag) = 1;                                    % density    
    gamma(flag) = 7;                                    % pressure exponent
               % artificial speed of sound, c_0 = 10 * v_max = 10 * sqrt(g*H)
    c_0(flag) = 10*v_max;
    p_0(flag) = rho_0(flag) * c_0(flag)^2 / gamma(flag);% reference pressure 
    Xi(flag) = 0.0 * p_0(flag);                         % background pressure
    my(flag) = 0.01;                                    % viscosity
    alpha(flag) = alpha;                                % artificial visc factor
    % initial velocities
    vel = [0,0];
    % create particle matrix
    fluid = create_fluid(dx, f_lowleft, f_lowright, f_upleft, f_upright);    
    particles = [];
    [particles, int_f1] = initialisation(particles, fluid, flag, rho_0(flag), dx, d, vel);

    %% second fluid phase
    multi_ph = 2;
    if multi_ph == 1
    f_lowleft =  [-0.2  0.0];
    f_lowright = [ 0.0  0.0];
    f_upleft =   [-0.2  1.2];
    f_upright =  [ 0.0  1.2];
    % properties of fluid (will be stored)
    flag = 3;
    rho_0(flag) = 0.25;                                    % density
    gamma(flag) = 7;                                    % pressure exponent
               % artificial speed of sound, c_0 = 10 * v_max = 10 * sqrt(g*H)
    c_0(flag) = 10*v_max;
    p_0(flag) = rho_0(flag) * c_0(flag)^2 / gamma(flag);% reference pressure 
    Xi(flag) = 0.05 * p_0(flag);                         % background pressure
    my(flag) = 0.01;                                    % viscosity
    alpha(flag) = 0.24;                                 % artificial visc factor
    % initial velocities
    vel = [0,0];
    % create particle matrix
    fluid = create_fluid(dx, f_lowleft, f_lowright, f_upleft, f_upright);
    [particles, int_f2] = initialisation(particles, fluid, flag, rho_0(flag), dx, d, vel);
    
    % integer of fluid particles in matrix
    int_fluid = 1:max(int_f2);
    else
        int_fluid = 1:max(int_f1);
    end
    
    
    

    %% specify coordinates of edges for boundary
    width = dx*round(5.366/dx);
    b_lowleft =  [-0.0   0.0];
    b_lowright = [ width  0.0];
    b_upleft =   [-0.0   3.0];
    b_upright =  [ width  3.0];
    % properties
    v_wall = [0,0]; % prescribed wall velocity
    a_wall = [0,0]; % wall acceleration
    flag = 1;       % for boundary
    %set artificial viscosity of boundary to 0
    alpha(flag) = 0;
    % create boundary matrix
    boundary = create_boundary(kernel, dx, b_lowleft, b_lowright, b_upleft, b_upright);
    [particles, int_boundary] = initialisation(particles, boundary, flag, 0, dx, d, [0,0]);
    particles(int_boundary,4) = v_wall(1);
    particles(int_boundary,5) = v_wall(2);
    
    %% create removable wall for hydrostatic pressure build up
    b_lowleft =  [ 2.0  0.0];
    b_lowright = [ 2.0+3*dx  0.0];
    b_upleft =   [ 2.0  1.2];
    b_upright =  [ 2.0+3*dx  1.2];
    bound_remove = create_fluid(dx, b_lowleft, b_lowright, b_upleft, b_upright);
    [particles, int_bound_remove] = initialisation(particles, bound_remove, flag, 0, dx, d, [0,0]);
   
    int_boundary = min(int_boundary):max(int_bound_remove);    
    
    %% plot initial positions of particles
    figure(1)
    hold on    
    plot(particles(int_f1,1), particles(int_f1,2), '.')
    if multi_ph == 1
    plot(particles(int_f2,1), particles(int_f2,2), 'g.')
    end
    plot(particles(int_boundary,1), particles(int_boundary,2), 'r.')
    plot([b_upleft(1,1), b_lowleft(1,1), b_lowright(1,1), b_upright(1,1)],...
        [b_upleft(1,2), b_lowleft(1,2), b_lowright(1,2), b_upright(1,2)], 'r', 'linewidth', 2)
    axis('equal')
   
   
end


