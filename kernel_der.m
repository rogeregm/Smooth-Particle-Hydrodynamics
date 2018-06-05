function [der_W] = kernel_der(kernel, d, h, q)

    if kernel == 1
        % Quintic Spline
        % cutoff radius r_c = 3 * h;

        % normalisation parameter
        if d == 2
            %alpha_d = -5 * 7 / (478 * pi * h^2); % for 2D
            alpha_d = -0.0233072092393989 / (h * h);
        elseif d == 3
            alpha_d = -5 / (120 * pi * h^3); % for 3D
        elseif d == 1
            alpha_d = -5 / (120 * h);
        end

        % derivative of Weighting function
        if q < 3 && q >= 2
            der_W = alpha_d * (3-q)^4;
        elseif q < 2 && q >= 1
            der_W = alpha_d * ((3-q)^4 - 6 * (2-q)^4);
        elseif q < 1
            der_W = alpha_d * ((3-q)^4 - 6 * (2-q)^4 + 15 * (1-q)^4);
        elseif q >= 3
            der_W = 0;
        end


    elseif kernel == 2
        % Quintic or Wendland Kernel
        % cutoff radius r_c = 2 * h;

        % normalisation parameter alpha_d
        if d == 2
            alpha_d = 7 / (4 * pi * h^2); % for 2D
        elseif d == 3
            alpha_d = 21 / (16 * pi * h^3); % for 3D
        end

        % Derivative of Weighting function
        if abs(q) >= 2
            der_W = 0;
        else
            der_W = alpha_d * (1-q/2)^3 * (-5*q);
        end

    end

end