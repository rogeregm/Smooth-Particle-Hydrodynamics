%% Surface Tension Force model
%  Roger Gonzalez
%  30/05/17

function [PF, Virial_Pressure, adhesion_F, der_color_field_i] = PairwiseForce(ST_model, particles, h, dist, a, b, domain,  rho_0, p_0, r_c, rho_b, m_b, Xi, gamma, eqn_of_state, int_boundary,idx_all)
% Have to individually import mass and density for the boundary particle


%% Cosine Pairwise Force
if ST_model == 1
    phase_flag = abs(particles(a,3) - particles(b,3));
    %Calibration constants for surface tension
    if length(domain{a}) < 27 || phase_flag ~= 0
        S_ab = (16^2)*0.00001;    %different fluid phases or surface particle
        VP_switch = 0;
    elseif phase_flag == 0 
        S_ab = (16^2)*0.1;    %same fluid phase or bulk particle
        VP_switch = 1; % To ensure virial pressure only adds when all particles are in same phase
    end

    if dist < h
        PF = -100*S_ab*cos(1.5*pi*dist / (3*h));
        Virial_Pressure = 0;
    else
        PF = 0;
        xi = (h^3) *(-8 +9*pi^2) / (27*pi^2);
        Virial_Pressure = - xi*(particles(a,5)^2)*S_ab * VP_switch;
    end
    %Specific energy due to the pairwise force model
%     PF_spec_E = (8/81) * ((h^4)/(pi^4)) * (9/4*pi^3 - 6*pi -4) * particles(a,5) * particles(b,5) * S_ab;
%     assignin(ws, 'PF_spec_E', PF_spec_E);
end

%% PF-3
if ST_model ==3
    phase_flag = abs(particles(a,3) - particles(b,3));
    if length(domain{a}) < 27 || phase_flag ~= 0 
        S_ab = (16^2)*2;    %same fluid phase
        VP_switch = 1; % To ensure virial pressure only adds when all particles are in same phase
    else
        S_ab = (16^2)*0.00001;    %different fluid phases
        VP_switch = 0;
    end

    if dist <= r_c
        eps=r_c/3.5;
        eps_0 = eps/2;
        A = (eps/eps_0)^3;
        psi_eps = exp(-(dist^2)/(2*(eps^2)));
        psi_eps_0 = exp(-(dist^2)/(2*(eps_0^2)));
        PF = S_ab*dist*(-A*psi_eps_0 + psi_eps);
        Virial_Pressure = 0;
    else
        PF = 0;
        xi = (r_c^3) *(-8 +9*pi^2) / (27*pi^2);
        Virial_Pressure = - xi*(particles(a,5)^2)*S_ab * VP_switch;
    end
    
end

%% Kernel Weighted Attraction (Becker & Teschner)
%   - Upon further consideration this model is pointless
%   - It is not useful with the kernel that we are using and less accurate
%     than Cohesion Force (Akinci)

if ST_model == 4
   
    W = kernel_fct(1, 2, h, dist/h);
    PF = -( particles(b,4) / particles(a,4) ) * W;
    
    Virial_Pressure = 0;
    PF_spec_E = 0;
    
    if dist < 0.5*r_c
        PF=0;
    end
end

%% Cohesion Force model

if ST_model == 5
    
%     NS = particles(:,1:2);          %List of all particles
%     idx_all = rangesearch(NS,NS,r_c);    %List of all particles within range of each fluid particle
    
    cohesion_spline = 0;
    adhesion_spline = 0;
    adhesion_F = 0;
    
    %Cohesion Force
    if dist <= r_c && dist > 0.5*r_c
        cohesion_spline = ((r_c-dist)^3)*(dist^3);        
    elseif dist >= 0 && dist <= 0.5*r_c
        cohesion_spline = 2*((r_c-dist)^3)*(dist^3) - (r_c^6)/64;
    end   
    cohesion_spline = (32/(pi*r_c^9))*cohesion_spline;
    cohesion_F  = -particles(a,4) * particles(b,4) * cohesion_spline;
    
    %Surface Area Minimization
    %derivative of the color field for a fluid particle: particle in
    %question
    der_color_field_i = der_color_field(particles, a, a, r_c, h, rho_0, p_0, Xi, gamma, eqn_of_state, idx_all);
    
    % derivative of the color field of second particle
    der_color_field_j = der_color_field(particles, a, b, r_c, h, rho_0, p_0, Xi, gamma, eqn_of_state, idx_all);
    
    surf_area_min_F = -particles(a,4) * (der_color_field_i - der_color_field_j);
    
    if isnan(surf_area_min_F) || isnan(der_color_field_i) || isnan(der_color_field_i)
        fprintf('mass of fluid particle: %d \n ',particles(a,4));
        fprintf('derivative of color field of a: %d \n', der_color_field_i);
        fprintf('derivative of color field of b: %d \n', der_color_field_j);
        error('Surface Area Minimization force is NaN')
    end
    
    %Adhesion Force
    if particles(b,3) == 1
        sum_boundary_vol = Boundary_Sampling_density(particles, 2, h, r_c, int_boundary, 1);
        temp_index = b - int_boundary(1)+1;
        BSD = sum_boundary_vol(temp_index);
       if dist > 0.5*r_c && dist <= r_c
           adhesion_spline = (0.007/(r_c^3.25)) * (-(4*dist^2)/r_c + 6*dist - 2*r_c)^(1/4);
       end
    adhesion_F = - particles(a,4) * rho_0(particles(a,3)) * BSD * adhesion_spline;
    end
    
    % Correction Factor to amplify forces of the particles with
    % neighbourhood deficiency
    corr_factor=(2*rho_0(particles(a,3))) / ( particles(a,5) + particles(b,5) );
    contact_angle = 120;
    gamma = 1-0.75*cosd(contact_angle);
    beta = 1+abs(0.5*cosd(contact_angle));
%     gamma=0.1;
%     beta=1.5;
    PF = corr_factor * gamma * (cohesion_F) + beta * adhesion_F;

% %   virial Pressure
%     Cosine Force Virial Pressure
%     eps = (h^3 *(9*(pi^2)-8))/(27*(pi^2));
%     Virial_Pressure = -eps * (particles(a,5)^2) * gamma

%     %Derived Virial Pressure Term
%     if dist <= r_c && dist > 0.5*r_c
%         eps = 1/1260;        
%     elseif dist >= 0 && dist <= 0.5*r_c
%         eps = 61/645120;
%     end  
%     Virial_Pressure = - pi * particles(a,5)^2 * gamma * particles(a,4) * particles(b,4) * r_c^3 * eps;
end
    Virial_Pressure=0;

end