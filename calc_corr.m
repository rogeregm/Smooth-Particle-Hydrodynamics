function [diss_L] = calc_corr(particles, idx, int_fluid, h, kernel, d)

%initialisation
diss_L(2,2,int_fluid) = 0;

%loop to calculate correction matrix L
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
            Fab = [drx;dry]/rad * DWab;
            
            Vb = particles(b,4) / particles(b,5);                      
            
            %correction matrix
            diss_L(:,:,a) = diss_L(:,:,a) + Vb * [-drx;-dry]*Fab';
            
        end %end if b fluid particle
    end %end loop over b
end %end loop over a

%remove ill-conditiond matrices
parfor a = int_fluid    
    if rcond(diss_L(:,:,a)) < 1.0e-14
%         a
%         bad_cond_mat = diss_L(:,:,a)
        diss_L(:,:,a) = eye(d);
    end
end