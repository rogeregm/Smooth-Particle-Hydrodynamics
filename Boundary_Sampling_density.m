function [sum_boundary_vol] = Boundary_Sampling_density(particles, d, h, r_c, int_boundary,kernel)
%To calculate the boundary particle density

Y = particles(int_boundary,1:2);
range = rangesearch(Y,Y,r_c);
sum_boundary_vol =zeros(1,length(int_boundary));

for n = 1:length(int_boundary)
    Bi = int_boundary(n);
    for k = range{n}(2:end)
        Bj = int_boundary(k);
        %distance between particles
        drx = particles(Bi,1) - particles(Bj,1);
        dry = particles(Bi,2) - particles(Bj,2);
        rad = sqrt(drx^2 + dry^2);
            
        % nondimensional distance
        q = rad / h;
 
        W = kernel_fct(kernel, d, h, q);
        sum_boundary_vol(n) = sum_boundary_vol(n) + W;
    end
end
sum_boundary_vol = 1./sum_boundary_vol ;

end