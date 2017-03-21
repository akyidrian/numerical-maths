function [u, n, err] = cn(alpha, del_t, x_int, u_bound, u0, N, tol)
    % C-N iteration method for 1D unsteady heat equation.
    % Note, u_bound = [u_left, u_right]
    del_x = abs(diff(x_int)) / N;
    d = (alpha*del_t) / (del_x)^2;

    I = N + 1;
    
    % Initial conditions.
    u(1:I,1) = u0;
    
    % Boundary Conditions.
    u(1) = u_bound(1);
    u(I) = u_bound(2);
    
    % Tridiagonal matrix corresponding to nth time step vector.
    n_tridiag = full(gallery('tridiag',I,d,2*(1-d),d));
    n_tridiag(1,1:I) = zeros(1,I);
    n_tridiag(1,1) = 1; % Ignoring boundary condition.
    n_tridiag(I,1:I) = zeros(1,I);
    n_tridiag(I,I) = 1; % Ignoring boundary condition.
    
    % Tridiagonal matrix corresponding to n+1 time step vector.
    np1_tridiag = full(gallery('tridiag',I,-d,2*(1+d),-d));
    np1_tridiag(1,1:I) = zeros(1,I);
    np1_tridiag(1,1) = 1;
    np1_tridiag(I,1:I) = zeros(1,I);
    np1_tridiag(I,I) = 1;
    
    n = 0;
    curr_err = tol + 1;
    err = [];
    while(curr_err > tol)
        u_new = np1_tridiag\n_tridiag*u;
        curr_err = max(abs(u_new - u));
        u = u_new;
        n = n + 1;
        err(n) = curr_err; % Storing a history of errors.
    end
end
