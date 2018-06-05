%% Function to calculate the energy of each fluid in a multi-fluid system

function [KE_alpha, KE_beta, PE_alpha, PE_beta] = Post_multifluid_energy()


%% initialisation
clear particles
%close all
clc

time = save_pos_t;
n = length(time);


E_kin_alpha = zeros(n,1);
E_kin_beta  = zeros(n,1);
E_pot_alpha = zeros(n,1);
E_pot_beta  = zeros(n,1);
E_pV_alpha  = zeros(n,1);
E_pV_beta   = zeros(n,1);
E_tot1_alpha = zeros(n,1);
E_tot2_alpha = zeros(n,1);
E_tot1_beta = zeros(n,1);
E_tot2_beta = zeros(n,1);

% particles = []
%% loop over timesteps
for t = 2:n
    
    % particle data of current timestep
    particles(:,:) = save_pos(int_fluid,:,t);
    
    % kinetic Energy
    mass = particles(:,4);
    vel_sq = (particles(:,6).^2 + particles(:,7).^2); % v^2
    kin = mass .* vel_sq / 2;
    kin(isnan(kin)) = 0;
    
    % summing the KE of all the particles of each fluid
    % E_kin_alpha is going to rewrite the value at t for every iteration in
    % this loop
    i=1;
    for i = 1:int_fluid
        if particles(i,3) == 2
            E_kin_alpha(t) = E_kin_alpha(t) + kin(i);
        elseif particles(i,3) == 3
            E_kin_beta(t) = E_kin_beta(t) + kin(i);
        end
    end
        
    if (time(t) - t_damp) <= 1e-3
        % transient gravity force
        g = 0.5 * (sin((-0.5 + time(t)/t_damp) * pi) + 1 ) * abs(G(2));
        
        % potential Energy
        pot = mass .* g .* particles(:,2);
        pot(isnan(pot)) = 0;
        % Summing PE of all particles of each fluid at this timestep
        i=1;
        for i = 1:int_fluid
            if particles(i,3) == 2
                E_pot_alpha(t) = E_pot_alpha(t) + pot(i);
            elseif particles(i,3) == 3
                E_pot_beta(t) = E_pot_beta(t) + pot(i);
            end
        end
%         E_pot(t) = sum(pot);
%         E_pot_ref = E_pot(t);

        % pressure volume work
        volume = particles(:,4)./ particles(:,5);
        rho = particles(:,5);
        volume_0 = dx * dx;
        pV = mass .* c_0(2).^2./(7.*(7-1)).*((rho./rho_0(2)).^(7-1) + (7-1).*rho_0(2)./rho - 7);
        %pV = volume .* particles(:,8);
        %pV =  p_0(2) ./ (1-gamma(2)) .*  volume_0.^(gamma(2)) .* ( volume.^(1-gamma(2)) - volume_0.^(1-gamma(2))) + (Xi(2) - p_0(2)) .* (volume - volume_0); 
        pV(isnan(pV)) = 0;
        i=1;
        for i = 1:int_fluid
            if particles(i,3) == 2
                E_pV_alpha(t) = E_pV_alpha(t) + pV(i);
            elseif particles(i,3) == 3
                E_pV_beta(t) = E_pV_beta(t) + pV(i);
            end
        end
%         E_pV(t) = sum(pV);
%         E_pV_ref = E_pV(t);        
        
    elseif  time(t) >= t_damp
        g = abs(G(2)); 
        
        % potential Energy
        pot = mass .* g .* particles(:,2);
        pot(isnan(pot)) = 0;
        % Summing PE of all particles of each fluid at this timestep
        j=0;
        for j = 1:int_fluid
            if particles(j,3) == 2
                E_pot_alpha(t) = E_pot_alpha(t) + pot(j);
            elseif particles(j,3) == 3
                E_pot_beta(t) = E_pot_beta(t) + pot(j);
            end
        end
%         E_pot(t) = sum(pot) - E_pot_ref;
        
        % pressure volume work
        volume = particles(:,4)./ particles(:,5);
        rho = particles(:,5);
        volume_0 = dx * dx;
        pV = mass .* c_0(2).^2./(7.*(7-1)).*((rho./rho_0(2)).^(7-1) + (7-1).*rho_0(2)./rho - 7);
        %pV = volume .* particles(:,8);
        %pV =  p_0(2) ./ (1-gamma(2)) .*  volume_0.^(gamma(2)) .* ( volume.^(1-gamma(2)) - volume_0.^(1-gamma(2))) + (Xi(2) - p_0(2)) .* (volume - volume_0); 
        pV(isnan(pV)) = 0;
        % Summing PE of all particles of each fluid at this timestep
        i=1;
        for i = 1:int_fluid
            if particles(i,3) == 2
                E_pV_alpha(t) = E_pV_alpha(t) + pV(i);
            elseif particles(i,3) == 3
                E_pV_beta(t) = E_pV_beta(t) + pV(i);
            end
        end
%         E_pV(t) = sum(pV) - E_pV_ref;
    end
    
    % total Energy
    E_tot1_alpha(t) = E_kin_alpha(t) + E_pot_alpha(t);
    E_tot2_alpha(t) = E_kin_alpha(t) + E_pot_alpha(t) + E_pV_alpha(t);
    E_tot1_beta(t)  = E_kin_beta(t) + E_pot_beta(t);
    E_tot2_beta(t)  = E_kin_beta(t) + E_pot_beta(t) + E_pV_beta(t);  
    
end

end

