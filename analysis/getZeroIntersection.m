function [x, y, intersects] = getZeroIntersection(x1a, y1a, x2a, y2a, x1b, y1b, x2b, y2b)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    numerator = x2b*(x2a*(y1a - y1b) + x1a*(y1b - y2a)) + x1b*(x1a*(y2a - y2b) + x2a*(-y1a + y2b));
    denominator = (-(x1b - x2b))*(y1a - y2a) + (x1a - x2a)*(y1b - y2b);
    
    x = numerator/denominator;
    y = y1a + (y2a-y1a)*(x-x1a)/(x2a-x1a);
    
    
    
    if (x > 0 && x < 1)
        intersects = true;
        
    else
        intersects = false;
    end
    
end

