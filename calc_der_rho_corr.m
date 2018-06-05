function [der_rho_L] = calc_der_rho_corr(particles, idx, int_fluid, h, kernel, d, diss_L)

der_rho_L(d,int_fluid) = 0;

parfor a = int_fluid % loop over all fluid particles
    
    % calculate interaction of one fluid particle with all other particles
    for b = idx{a}(2:end)
        
        flag_b = particles(b,3);
        if flag_b ~= 1 % if b is a fluid particle
            %distance between particles
            drx = particles(a,1) - particles(b,1);
            dry = particles(a,2) - particles(b,2);
            rad = sqrt(drx^2 + dry^2);
            
            % nondimensional distance
            q = rad / h;
            
            %kernel and derivative values
            DWab = kernel_der(kernel, d, h, q) / h;
            
            %kernel derivative with respect to x_a
            Fab = [(drx/rad);(dry/rad)]*DWab;
            
            %corrected kernel
            CFab = diss_L(:,:,a)\Fab;
            
            %define particle properties
            Vb = particles(b,4) / particles(b,5);
          
            der_rho_L(:,a) = der_rho_L(:,a) + (CFab * (particles(b,5)-particles(a,5)) * Vb);
            
        end %if fluid particle
    end %end loop over b
end %end loop over a