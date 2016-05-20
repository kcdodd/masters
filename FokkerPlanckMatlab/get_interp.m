function y = get_interp(y_i, x_min, dx, x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    n = max(size(y_i, 1),size(y_i, 2));
    
    % assume all elements same interval
    
    a = floor((x-x_min)/dx);
    nx = a+1;
    
    if (nx > n || nx < 1)
        disp(['value not found: ' num2str(x)]);
    elseif (nx == n)
        y = y_i(n);
    else
        if (nx == 1 || nx == n-1)
            
            % do linear interp for first and last segments
            y = (y_i(nx+1)-y_i(nx))*(x-(a*dx+x_min))/dx + y_i(nx);
        else
            % do cubic interp for all other segments
            alpha = (x-(a*dx+x_min))/dx;
        
            a3 = (3*(y_i(nx) - y_i(nx+1))+y_i(nx+2)-y_i(nx-1))/6;
            a2 = (y_i(nx+1)+y_i(nx-1)-2*y_i(nx))/2;
            a1 = (6*y_i(nx+1) - 3*y_i(nx) - y_i(nx+2) - 2*y_i(nx-1))/6;
            
            y = a3*alpha^3 + a2*alpha^2 + a1*alpha + y_i(nx);
        end
        
        
    end
    
end
