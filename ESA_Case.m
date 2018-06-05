function [particles, rho_0,gamma,c_0,p_0,Xi,my,alpha, a_wall, int_fluid, int_boundary] = ESA_Case(kernel, dx, d, v_max, alpha)

    % rectangular domain
    % first fluid phase
    % specify coordinates of edges for fluid
    width = dx*round(0.295/dx); % in m
    height = dx*round(0.03/dx);
    f_lowleft =  [-width  0.0];
    f_lowright = [ width  0.0];
    f_upleft =   [-width  height];
    f_upright =  [ width  height];
    % properties of fluid (will be stored)
    flag = 2;
    rho_0(flag) = 1;                                    % density  (water)  
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
    int_fluid = 1:max(int_f1); 
    

    %% specify coordinates of edges for boundary
    height = dx*round(0.25/dx);
    b_lowleft =  [-width  0.0];
    b_lowright = [ width  0.0];
    b_upleft =   [-width  height];
    b_upright =  [ width  height];
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
    
    %% plot initial positions of particles
    figure(1)
    hold on    
    plot(particles(int_f1,1), particles(int_f1,2), '.')
    plot(particles(int_boundary,1), particles(int_boundary,2), 'r.')
    plot([b_upleft(1,1), b_lowleft(1,1), b_lowright(1,1), b_upright(1,1)],...
        [b_upleft(1,2), b_lowleft(1,2), b_lowright(1,2), b_upright(1,2)], 'r', 'linewidth', 2)
    axis('equal')
%    
   
end


