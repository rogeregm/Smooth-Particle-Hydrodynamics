function [boundary] = create_boundary(kernel, dx, b_lowleft, b_lowright, b_upleft, b_upright)

    % number of boundary particles for kernel support
    if kernel == 2
        k = 2; % wendland
    elseif kernel == 1
        k = 3; % quintic spline
    elseif kernel == 3
        k = 5; % for free bubble case
    end
    global testcase
    % create boundary particles vector
    no_rows = round(abs(b_upleft(1,2) - b_lowleft(1,2)) / dx) ;
    no_cols = round(abs(b_upleft(1,1) - b_upright(1,1)) / dx + 2*k);
    for n = 1 : no_rows
        % x-coordinate
        boundary(1+(n-1)*2*k : 2*n*k,1) = [b_upleft(1,1)-(k*dx-dx/2):dx:b_upleft(1,1)-dx/2, b_upright(1,1)+dx/2:dx:b_upright(1,1)+(k*dx-dx/2)];
        % y-coordinate    
        boundary(1+(n-1)*2*k : 2*n*k,2) = b_upleft(1,2) - (n-1/2) * dx;    
    end
    for n = 1 : k
        % x-coordinate
        boundary(2*no_rows*k+1+(n-1)*no_cols : 2*no_rows*k+n*no_cols,1) = b_lowleft(1,1)-(k*dx-dx/2):dx:b_lowright(1,1)+(k*dx-dx/2);
        % y-coordinate
        boundary(2*no_rows*k+1+(n-1)*no_cols : 2*no_rows*k+n*no_cols,2) = b_lowleft(1,2) - (n-1/2) * dx;
    end
    
    if  testcase == 6
        for n= 1:k
            % x-coordinate
            boundary(2*no_rows*k+1+(k+(n-1))*no_cols : 2*no_rows*k+(k+n)*no_cols,1) = b_lowleft(1,1)-(k*dx-dx/2):dx:b_lowright(1,1)+(k*dx-dx/2);
            % y-coordinate
            boundary(2*no_rows*k+1+(k+(n-1))*no_cols : 2*no_rows*k+(k+n)*no_cols,2) = b_upleft(1,2) + (n-1/2) * dx;
        end
    end
end