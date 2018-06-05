%% Energy evolution for postprocessing

%% initialisation
clear particles
%close all
clc

time = save_pos_t;
n = length(time);

% results = dir('particle_results/*.mat');

E_kin = zeros(n,1);
E_pot = zeros(n,1);
E_pV  = zeros(n,1);
E_tot1 = zeros(n,1);
E_tot2 = zeros(n,1);

% particles = []
%% loop over timesteps
for t = 2:n
    
        if n<10
        filename=['save_pos00',num2str(n)];
    elseif n<100
        filename=['save_pos0',num2str(n)];
    else
        filename=['save_pos',num2str(n)];
    end 

    % load fluid
    clear particles
    load([filename,'.mat'])
    particles=eval(filename);
    
    % particle data of current timestep
%     particles(:,:) = save_pos(int_fluid,:,t);

    
    % kinetic Energy
    mass = particles(:,4);
    vel_sq = (particles(:,6).^2 + particles(:,7).^2); % v^2
    kin = mass .* vel_sq / 2;
    kin(isnan(kin)) = 0;
    E_kin(t) = sum(kin);
        
    if (time(t) - t_damp) <= 1e-3
        % transient gravity force
        g = 0.5 * (sin((-0.5 + time(t)/t_damp) * pi) + 1 ) * abs(G(2));
        
        % potential Energy
        pot = mass .* g .* particles(:,2);
        pot(isnan(pot)) = 0;
        E_pot(t) = sum(pot);
        E_pot_ref = E_pot(t);

        % pressure volume work
        volume = particles(:,4)./ particles(:,5);
        rho = particles(:,5);
        volume_0 = dx * dx;
        pV = mass .* c_0(2).^2./(7.*(7-1)).*((rho./rho_0(2)).^(7-1) + (7-1).*rho_0(2)./rho - 7);
        %pV = volume .* particles(:,8);
        %pV =  p_0(2) ./ (1-gamma(2)) .*  volume_0.^(gamma(2)) .* ( volume.^(1-gamma(2)) - volume_0.^(1-gamma(2))) + (Xi(2) - p_0(2)) .* (volume - volume_0); 
        %pV(isnan(pV)) = 0;
        E_pV(t) = sum(pV);
        E_pV_ref = E_pV(t);        
        
    elseif  time(t) >= t_damp
        g = abs(G(2)); 
        
        % potential Energy
        pot = mass .* g .* particles(:,2);
        pot(isnan(pot)) = 0;
        E_pot(t) = sum(pot) - E_pot_ref;
        
        % pressure volume work
        volume = particles(:,4)./ particles(:,5);
        rho = particles(:,5);
        volume_0 = dx * dx;
        pV = mass .* c_0(2).^2./(7.*(7-1)).*((rho./rho_0(2)).^(7-1) + (7-1).*rho_0(2)./rho - 7);
        %pV = volume .* particles(:,8);
        %pV =  p_0(2) ./ (1-gamma(2)) .*  volume_0.^(gamma(2)) .* ( volume.^(1-gamma(2)) - volume_0.^(1-gamma(2))) + (Xi(2) - p_0(2)) .* (volume - volume_0); 
        pV(isnan(pV)) = 0;
        E_pV(t) = sum(pV) - E_pV_ref;
    end
    
    % total Energy
    E_tot1(t) = E_kin(t) + E_pot(t);
    E_tot2(t) = E_kin(t) + E_pot(t) + E_pV(t);    
    
end

%% plotting
if testcase == 3
    time1 = (time-t_remove)/t_ref;
    range = 0.55;
    t1 = 0;
    E_ref = rho_0(2) * abs(G(2)) * 2 * H^2 *(5.366-2) / (2 * 5.366);
else %testcase == 2
    time1 = time/t_ref;
    range = 1e-2;
    t1 = 1;
    E_ref = 0.5 * rho_0(2) * H^3;
    %range = max(max(E_kin(200:end)),max(E_pot(200:end)));
end
%{
%subplot(2,1,1)
hold on
grid on
plot(time/t_ref,E_tot1,'k')
plot(time/t_ref,E_kin ,'g')
plot(time/t_ref,E_pot,'b')
legend('total(mechanical)','kinetic','potential','Location','EastOutside')
axis([0 max(time/t_ref) -range range])

subplot(2,1,2)
%}

figure
fz = 9;
width = 1;
lw = 1.5;
hold on
grid on
box on
%set(gcf, 'Position', (540*width + 20)*[0.2, 0.2, 1.0, 0.30]);
plot(time1(101:end),E_pot(101:end)/E_ref,'ob','linestyle','-','linewidth',lw);
plot(time1(101:end),E_pV(101:end)/E_ref ,'color','r','linestyle','-','linewidth',lw);
plot(time1(101:end),E_kin(101:end)/E_ref ,'og','linestyle','-','linewidth',lw);
plot(time1(101:end),E_tot2(101:end)/E_ref,'k','linewidth',1.2*lw);
%legend('E_{ext} -  E_{ext}^{(0)}','E_{int} -   E_{int}^{(0)}','E_{tot} -  E_{tot}^{(0)}','10 E_{kin}','Location','EastOutside')
%legend('(E_{tot} -  E_{tot}^{(0)}) / E_{ref}','(E_{kin} -  E_{kin}^{(0)}) / E_{ref}','(E_{ext} -  E_{ext}^{(0)}) / E_{ref}','(E_{int} -   E_{int}^{(0)}) / E_{ref}','Location','EastOutside')
%legend('E_{ext} -  E_{ext}^{(0)}','E_{int} -   E_{int}^{(0)}','Location','EastOutside')
%legend('E_{tot} -  E_{tot}^{(0)}','10 * E_{kin} -  E_{kin}^{(0)}','Location','EastOutside')
%axis([0 7.5 min(E_pot)/E_ref max(E_kin)/E_ref])
xlabel('t / t_{ref}', 'fontsize',fz)
ylabel('Energy', 'fontsize',fz)
set(gca, 'fontsize',fz) 

% logarithmic plot
%{
fz = 9;
width = 1;
lw = 1.5;
figure
semilogy(time1(121:end),-E_tot1(121:end),'k','linewidth',1.2*lw)
hold on
semilogy(time1(121:end),E_kin(121:end) ,'color',[0 0.7 0],'linestyle','-','linewidth',lw);
semilogy(time1(121:end),-E_pot(121:end),'color','b','linestyle','-','linewidth',lw);
semilogy(time1(121:end),-E_pV(121:end) ,'color','r','linestyle','-','linewidth',lw);
legend('-(E_{tot} -  E_{tot}^{(0)}) / E_{ref}','(E_{kin} -  E_{kin}^{(0)}) / E_{ref}','-(E_{ext} -  E_{ext}^{(0)}) / E_{ref}','-(E_{int} -   E_{int}^{(0)}) / E_{ref}','Location','EastOutside')
axis([t1 7.5 -range range])
xlabel('t / t_{ref}', 'fontsize',fz)
ylabel('Energy', 'fontsize',fz)
set(gca, 'fontsize',fz) 
set(gca,'YDir','reverse')
set(gcf, 'Position', (540*width + 20)*[0.2, 0.2, 1.0, 0.40]);
grid on
box on

%}