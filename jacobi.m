%TODO general boundary intervals
function [u] = jacobi(x_int, y_int, u_bound, N, tol)
% Jacobi iteration method for elliptic PDEs (e.g. 2D laplace).
% u_bound = [u_left, u_top, u_right, u_bottom]
    del_x = abs(diff(x_int)) / N(1);
    del_y = abs(diff(y_int)) / N(2);
    beta = del_x / del_y;

    u0 = mean(u_bound) / 4;
    
    I = N(1) + 1;
    J = N(2) + 1;
    u(1:I, 1:J) = u0;
    %TODO: need to take into account intervals for bounds
    u(:, 1) = u_bound(1); % left
    u(I, :) = u_bound(2); % top
    u(:, J) = u_bound(3); % right
    u(1, :) = u_bound(4); % bottom
    u_new = u;
    k = 0;
    err = tol + 1; % To begin iterations.
    while(err > tol)
        err = 0;
        for i=2:I-1
            for j=2:J-1
                u_new(i,j) = (1/(2*(1+beta^2))) * ...
                (u(i+1,j) + u(i-1,j) + (beta^2)*(u(i,j+1) + u(i,j-1)));
                err_new = abs(u_new(i,j) - u(i,j));
                if (err < err_new)
                   err = err_new; 
                end
            end
        end
        u = u_new; % TODO: Don't need to copy entire matrix.
        k = k + 1;
    end
end

