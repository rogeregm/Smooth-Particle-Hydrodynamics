%% Calculating the gradient of the smoothed color field
%part of the Surface area minimization in Akinci Surface tension model
%
%   Roger Gonzalez
%   04/07/2017

function [boundary_der_color_field] = der_color_field(particles, a, b, r_c, h, rho_0, p_0, Xi, gamma, eqn_of_state, range)
%
% a - fluid particle wrt which boundary values will be calculated
% b - particle in question
% r_c - kernel sampling range
% h - initial distance between particles = dx


boundary_der_color_field = 0;

    for k = range{b}(2:end)
        
        if particles(k,3) == 1
            rho_b = boundary_rho(particles, a, rho_0, p_0, Xi, gamma,eqn_of_state);
            m_b = rho_0(particles(a,3)) * h*h;
        else
            rho_b = particles(k,5);
            m_b = particles(k,4);
        end
        %distance between particles
        drx = particles(b,1) - particles(k,1);
        dry = particles(b,2) - particles(k,2);
        rad = sqrt(drx^2 + dry^2);
        W_der =  kernel_der(1, 2, h, rad)/h;

        boundary_der_color_field = boundary_der_color_field + (m_b/rho_b) *W_der;
    end

    % Serves to make the term scale independent
    boundary_der_color_field = r_c * boundary_der_color_field;
    
    %Warning Signs
    if isnan(boundary_der_color_field)
            fprintf('Mass \n',m_b);
            fprintf('desity \n',rho_b);
            error('derivative of the color field is NaN')
    end
end