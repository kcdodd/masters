function [x_out, y_out] = findZeroIntersection(x_in, y_in, f1, f2)
%FINDZEROINTERSECTION Finds points where both functions are simultaneously
%zero. 
%   Uses bi-linear interpolation.

    %make sure size of array is correct

    nx1 = size(f1, 2);
    nx2 = size(f2, 2);
    nx = size(x_in, 2);

    ny1 = size(f1, 1);
    ny2 = size(f2, 1);
    ny = size(y_in, 2);

    if (~(nx1 == nx && nx2 == nx && ny1 == ny && ny2 == ny))
        disp('array sizes must match');
        return;
    end

    % go cell by cell
    for ix = 1:(nx-1)
        for iy = 1:(ny-1)
            iszero = 0;
            % first see if there are any zeros in f1 in this cell.
            % find intercepts along each side of the cell. If all
            % intercepts are outisde cell boundaries then there are no
            % zeros inside cell.
            
            [haszero_a, zeropts_a] = getZeroPts(f1(iy,ix), f1(iy,ix+1), f1(iy+1,ix), f1(iy+1,ix+1));
            
            if (haszero_a)
                %f1 has zeros in cell somewhere. Now see if they intersect
                %the zeros of f2.
                
                [haszero_b, zeropts_b] = getZeroPts(f2(iy,ix), f2(iy,ix+1), f2(iy+1,ix), f2(iy+1,ix+1));
                
                if (haszero_b)
                    [zero_x, zero_y, intersects] = getZeroIntersection(zeropts_a(1,1), zeropts_a(1,2), zeropts_a(2,1), zeropts_a(2,2), zeropts_b(1,1), zeropts_b(1,2), zeropts_b(2,1), zeropts_b(2,2));
                    
                    if (intersects)
                        x_out = zero_x*(x_in(ix+1)-x_in(ix)) + x_in(ix);
                        y_out = zero_y*(y_in(iy+1)-y_in(iy)) + y_in(iy);
                        
                        return;
                    end

                end
                
            end
            
        end
    end

end

