% SPH solver
% Jan Eichstaedt, CID 00883791
% Imperial College London, Department of Aeronautics
% Begin: 29.05.2014
% Current Version: 26.06.2014
%
% Roger Gonzalez,CID:00839673
% Imperial College London, Department of Aeronautics
% Surface tension and testcases 5-9
% 01.05.2017 - 15.09.2017

%% 
close all
clear
clc

%% Initialisation

% Start new calculation?
cont = input('Start new calculation? \n yes - 1, no - 2 , Drop Equilibrium - 3 \n');
% cont =1;
global testcase
if cont == 1  
    % setup variables and domain etc...
    [d, dx, h, kernel, time_nondim, viscosity, delta_sph, slip_cond,...
    scheme, testcase, dt_save,particles, rho_0,gamma,c_0,p_0,Xi,my,alpha,...
    a_wall, int_fluid, int_boundary, int_bound_remove, t_remove,t_damp_fact,...
    G, Height, vel_wall_cond, eqn_of_state, delta] = SPH_input();
    PF_sum(int_fluid) = 0;
    % reference values  
    H = max(particles(int_fluid,2)) - min(particles(int_fluid,2)) + dx;% height of fluid
    if abs(H - Height) >= 1e-06
        disp('wrong height specified');
%         break
    end
    v_ref = sqrt(abs(G(2)) * H );
    t_ref = H / v_ref; 
    p_ref = rho_0 * abs(G(2)) * H;

    % time initialisation
    t = 0;      % time
    n_dt = 0;   % number of current timestep
    t_damp = t_damp_fact * t_ref; % time over which gravitiy is introduced
    t_end = time_nondim * t_ref;
    
    % saving initialisation
    n_save = ceil(time_nondim/dt_save);
    if n_save >= 1001
        disp('Too many saving steps specified!');
%         break
    end
    %save_pos = zeros(size(particles,1),size(particles,2),n_save);
    save_pos_t = zeros(1,n_save);
    n_save = 1;
    % initial particles values
    particles_ini = particles;
   

elseif cont == 2
    % load previous workspace
    try
        load save_var_all;
    catch
        disp('Cannot load file save_var_all.')
%         break
    end
    %t_end = 8 * t_ref; %manually adjust t_end if different from SPH_input.m
    
elseif cont == 3 
    try
        load save_DROP_EQUIL;
    catch
        disp('Cannot load file save_DROP_EQUIL.');
    end
    t=0;
end

% sloshing
if testcase == 4
    freq = 0.4;
    disp = 0.02;
    om = 2 * pi() * freq;
    amp = disp * om * om;
end

alpha
delta
tic
%% Simulate Particle Movement

while t <= t_end % loop over timestep
    % transient gravity force
    if t_damp_fact == 0 || t >= t_damp
        g = G;        
    else
        g = 0.5 * (sin((-0.5 + t/t_damp) * pi) + 1 ) * G;
    end
    
    % sloshing for ESA Case
    if testcase == 4 && t >= t_damp        
        g(1) = amp * sin(om * t);
    end
    
    % remove temporary boundaries, for dambreak case
    if (t-t_remove) >= 0 && testcase == 3
        int_boundary = min(int_boundary):min(int_bound_remove)-1;        
        particles = particles([int_fluid int_boundary],1:8);
    end
    
    % display progress
    prog = floor(100*t/t_end);    
    fprintf('progress = %d%%, time %f, timestep no %d\n',prog, t, n_dt); 
    
    % determine time step size (eq19,20,21)
    dt = stepsize(c_0, particles, h, alpha, d, g, viscosity, my);
    t = t + dt;
    n_dt = n_dt + 1;

    % time stepping schemes    
    if scheme == 1 % velocity- Verlet
        % (eq 14)
        if n_dt == 1
            if testcase ~= 6
            [particles(int_boundary,8),particles(:,6),particles(:,7)] = boundary_update(particles, int_fluid, int_boundary, kernel, d, h, g, slip_cond,PF_sum);
            end
            [dvelx_by_dt, dvely_by_dt,PF_sum,curv_a] = velocity_der(particles, int_fluid, h, kernel, d, alpha, c_0, g, viscosity, slip_cond, rho_0, p_0, Xi, gamma, my, eqn_of_state, int_boundary);
        end
        % Velocity Update
        particles(int_fluid,6) = particles(int_fluid,6) + dt/2 * dvelx_by_dt';
        particles(int_fluid,7) = particles(int_fluid,7) + dt/2 * dvely_by_dt';
        
        % (eq 15) Position Update
        particles(int_fluid,1) = particles(int_fluid,1) + dt/2 * particles(int_fluid,6);
        particles(int_fluid,2) = particles(int_fluid,2) + dt/2 * particles(int_fluid,7);
        
        % (eq 16) % Density and Pressure update
        if testcase ~=6
            [particles(int_boundary,8),particles(:,6),particles(:,7)] = boundary_update(particles, int_fluid, int_boundary, kernel, d, h, g, slip_cond,PF_sum);
        end
        [drho_by_dt, dissipation] = density_der(particles, int_fluid, h, kernel, d, rho_0, p_0, Xi, gamma, dx, vel_wall_cond, eqn_of_state, delta_sph);
        particles(int_fluid,5) = particles(int_fluid,5) + dt * drho_by_dt' + dt * delta * h * c_0(2) * dissipation';
        particles(int_fluid,8) = pressure(p_0,rho_0,gamma,Xi,particles, int_fluid, eqn_of_state);
        
        % (eq 17) Position Update
        particles(int_fluid,1) = particles(int_fluid,1) + dt/2 * particles(int_fluid,6);
        particles(int_fluid,2) = particles(int_fluid,2) + dt/2 * particles(int_fluid,7);
        
        % (eq 18) Velocity update 
        if testcase ~= 6
            [particles(int_boundary,8),particles(:,6),particles(:,7)] = boundary_update(particles, int_fluid, int_boundary, kernel, d, h, g, slip_cond,PF_sum);
        end
        [dvelx_by_dt, dvely_by_dt, diss_x, diss_y,PF_sum,curv_a] = velocity_der(particles, int_fluid, h, kernel, d, alpha, c_0, g, viscosity, slip_cond, rho_0, p_0, Xi, gamma, my, eqn_of_state, int_boundary);
        particles(int_fluid,6) = particles(int_fluid,6) + dt/2 * dvelx_by_dt';
        particles(int_fluid,7) = particles(int_fluid,7) + dt/2 * dvely_by_dt';
    end
        
    % saving particle data in time saving increments 'dt_save'
    if t >= n_save*dt_save*t_ref 
        % save time
        save_pos_t(n_save) = t;
        % create file name and content (only fluid particles)
        namestr = ['save_pos',num2str(n_save,'%3.3d')];
        eval([namestr,' = particles(int_fluid,1:8);']);
        % create directory
        dirname='particle_results';
        dirstatus=exist(dirname,'dir');
        if(dirstatus==0)
            mkdir(dirname)
        end
        % save fluid particle matrix as .mat-file
        save([dirname,'/',namestr],namestr);        
        % save fluid particle matrix as .vtu-file
        eval(['save_vtu(',namestr,',n_save, dirname)']);        
        % save workspace
        save ('save_var_all')
        % delete saved matrix from workspace
        clear(namestr)
        n_save = n_save + 1;
    end
end

t_sim = toc;

%% final savings
% reconstruct particle saving matrix
save_pos = zeros(size(particles,1),size(particles,2),n_save-1);
for i = 1:n_save-1
    namestr = ['save_pos',num2str(i,'%3.3d')];
    eval(['load ',dirname,'/',namestr,'.mat']);    
    eval(['save_pos(int_fluid,1:8,i) = ',namestr,';']);
    if testcase~=6
        save_pos(int_boundary,1:2,i) = particles(int_boundary,1:2);
    end
    eval(['clear ',namestr]);
end
if testcase ~= 6
    % save boundary particle matrix as .vtu-file
    namestr = 'boundary';
    eval([namestr,' = particles(int_boundary,1:8);']);
    eval(['save_vtu(',namestr,',0, dirname)']);
    % save workspace
    str = ['save_',datestr(clock,30)];
    save(str)
end
    fprintf('progress = 100%%, Simulation completeda after %d sec \n',t_sim);

