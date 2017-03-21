function [u, n, err] = ftcs(alpha, del_t, x_int, u_bound, u0, N, tol)
    % FTCS iteration method for 1D unsteady heat equation.
    % Note, u_bound = [u_left, u_right]
    del_x = abs(diff(x_int)) / N;
    d = (alpha*del_t) / (del_x)^2;

    I = N + 1;
    
    % Initial conditions.
    u(1:I,1) = u0;
    
    % Boundary conditions.
    u(1) = u_bound(1);
    u(I) = u_bound(2);
    
    u_new = u;
    n = 0;
    curr_err = tol + 1; % To begin iterations.
    err = [];
    while(curr_err > tol)
        % i=1 and i=I are the boundary.
        for i=2:I-1
            u_new(i) = u(i) + d*(u(i+1) - 2*u(i) + u(i-1));
        end
        curr_err = max(abs(u_new - u));
        
        u = u_new;
        n = n + 1; % Counting time intervals.
        err(n) = curr_err; % Storing a history of errors.
    end
end