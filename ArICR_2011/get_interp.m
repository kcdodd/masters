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
        y = (y_i(nx+1)-y_i(nx))*(x-(a*dx+x_min))/dx + y_i(nx);
    end
    
end
