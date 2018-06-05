function [fluid] = create_fluid(dx,f_lowleft, f_lowright, f_upleft, f_upright)
    
    % create fluid particle vector
    no_rows = round(abs(f_upleft(1,2) - f_lowleft(1,2)) / dx) ;
    no_cols = round(abs(f_upleft(1,1) - f_upright(1,1)) / dx) ;
    for n = 1 : no_rows
        % x-coordinate
        fluid(1+(n-1)*no_cols : n*no_cols,1) = f_lowleft(1,1)+dx/2:dx:f_lowright(1,1)-dx/2;
        % y-coordinate
        fluid(1+(n-1)*no_cols : n*no_cols,2) = f_lowleft(1,2) + (n-1/2) * dx;
    end

end