%% Time evolution plots

%close all

% pressure at fixed particle(1) or fixed position(2) ?
type = 2;
width = 0.75;
% position
if testcase == 2
    pos = [1/2*h 1/2*h];
elseif testcase == 3
    pos = [5.366 0.2];
    %pos = [ 5*h 5*h]
    load dambreak_ex_data
end

%% time evolution of pressure of one particle
if type == 1
    
    if testcase == 2    
        idx = rangesearch(particles(:,1:2), pos, h);
    elseif testcase == 3
        idx = rangesearch(particles(int_boundary,1:2), pos, h);
    end
    pos = particles(idx{1}(1),1:2)
    n=size(save_pos_t,2);
    pres(1:n,1) = save_pos(idx{1}(1),8,1:n);    % pressure
    yy(1:n) = save_pos(idx{1}(1),2,1:n);        % y-position
    time = save_pos_t;                          % time
    % transient gravity force
    idx = find(abs(time - t_damp) < 1e-3);
    g(1:idx(1)) = 0.5 * (sin((-0.5 + time(1:idx(1))/t_damp) * pi) + 1 ) * G(2);
    g(idx(1)+1 : n) = G(2);

    if testcase == 2
        % analytical values
        yy_ana = H;
        pres_ana = (rho_0(2) .* abs(g) .* (yy_ana-yy))';
    elseif testcase == 3
        %experimental values
        load dambreak_ex_data
        time_ex = dambreak_exp_data(:,1);
        pres_ex = dambreak_exp_data(:,2);
    end
%% time evolution of pressure at one position  
elseif type == 2      
    nn = size(save_pos_t,2);
    for t = 1:nn;

        sum_rhoW = 0;
        sum_pW = 0;
        sum_W = 0;
        if testcase == 2
            idx = rangesearch(save_pos(:,1:2,t), pos, 3*h);
        elseif testcase == 3
            idx = rangesearch(save_pos(int_fluid,1:2,t), pos, 2*3*h);
        end
        for m = idx{1}

            %distance between particles
            drx = pos(1) - save_pos(m,1,t);
            dry = pos(2) - save_pos(m,2,t);
            rad = sqrt(drx^2 + dry^2);

            % nondimensional distance
            q = rad /(2* h);
            W = kernel_fct(kernel, d, h, q);
            %sum_rhoW = sum_rhoW + save_pos(m,5,t) * W;
            sum_pW = sum_pW + save_pos(m,8,t) * W;
            sum_W = sum_W + W;        
        end
        press1(t) = sum_pW / sum_W;
        %rho2(t) = sum_rhoW / sum_W;
    end
    %press2 = p_0(2) * ( (rho2/rho_0(2)).^gamma(2) - 1 ) + Xi(2);

    time = save_pos_t;  % time
    % transient gravity force
    idx = find(abs(time - t_damp) < 1e-3);
    g(1:idx(1)) = 0.5 * (sin((-0.5 + time(1:idx(1))/t_damp) * pi) + 1 ) * G(2);
    g(idx(1)+1 : nn) = G(2);

    pres = press1';

    if testcase == 2
        % analytical values
        pres_ana = (rho_0(2) .* abs(g) .* (H-pos(2)));
    elseif testcase == 3
        %experimental values    
        time_ex = dambreak_exp_data(:,1);
        pres_ex = dambreak_exp_data(:,2);
    end
end

%{
%plot(time, press1-press2)

press(1,:) = press1;
press(2,:) = press2;
press(3,:) = press3;
dt = 5;
nt = 870/dt;
for n = 1:nt
    press_int(1:3,n) = mean(press(1:3,(n-1)*dt+1:dt*n),2);
    time_int(n) = time( (n-1) * dt + round(dt/2));
end

figure
hold on
plot((time_int-t_remove)/t_ref, press_int/p_ref(2))
plot((time_int-t_remove)/t_ref, press_int(3,:)/p_ref(2),'linewidth',2,'color', 'r')
xlabel('non-dim. time')
ylabel('non-dim. pressure')
grid on
legend
axis([1.8 7.2 -0.5 1.5])


%}
%% plotting

figure
set(gcf, 'Position', (540*width + 20)*[0.5, 0.5, 1.0, 0.75]); 
hold on
grid on
box on
p_sim = pres /p_ref(2);
xlabel('t / t_{ref}')
ylabel('p / p_{ref}')
if testcase == 3
    time1 = (time-t_remove)/t_ref;
    %interpolation
    dt = 10;
    nt = 870/dt;
    for n = 1:nt
        press_int(n) = mean(p_sim((n-1)*dt+1:dt*n));
        press_int(isnan(press_int)) = 0;
        time_int(n) = time1( (n-1) * dt + round(dt/2));
    end
    plot(time1, p_sim,'color',[1 1 1] * 0.5,'linestyle','-','linewidth',1)
    plot(time_int, press_int,'color','k','Marker','o','linewidth',1.5,'linestyle','none')
    p_ex = pres_ex/p_ref(2);
    plot(time_ex, p_ex,'linestyle','none','color','k','Marker','+','linewidth',1.5)
    legend( 'Simulation','Simulation, smoothed','Experiment','Location','Best')    
    axis([2 7.4 -0.5 1.5])    
elseif testcase == 2
    time1 = (time)/t_ref;
    plot(time1, p_sim,'color','k','linestyle','-','linewidth',1.5)
    p_ana = pres_ana/p_ref(2);
    plot(time1, p_ana,'linewidth',2,'color',[1 1 1] * 0.5,'linestyle','--')
    legend( 'Simulation','Analytical','Location','Best')
    axis([0 8 0.9775 0.995])
end

