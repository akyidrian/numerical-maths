function [u, k, del_x, del_y, err] = lsor(w, x_int, y_int, u_bound, u0, N, tol)
% LSOR iteration method for solving 2D Laplace equation.
% Note, boundary condition vector u_bound has the following format:
% u_bound = [u_left, u_top, u_right, u_bottom].
    del_x = abs(diff(x_int)) / N(1);
    del_y = abs(diff(y_int)) / N(2);
    beta = del_x / del_y;
      
    I = N(1) + 1;
    J = N(2) + 1;
    u(1:I, 1:J) = u0;
    
    % Boundary conditions.
    u(1, :) = u_bound(1); % left BC
    u(I, :) = u_bound(3); % right BC
    u(:, J) = u_bound(2); % top BC
    u(:, 1) = u_bound(4); % bottom BC

    % Creating tridiagonal matrix corresponding to the k+1 iteration vector.
    kp1_tridiag = full(gallery('tridiag', I, w, -2*(1+beta^2), w));
    kp1_tridiag(1,1:I) = zeros(1,I);
    kp1_tridiag(1,1) = 1; % Ignoring boundary condition.
    kp1_tridiag(I,1:I) = zeros(1,I);
    kp1_tridiag(I,I) = 1; % Ignoring boundary condition.
    
    u_new = u;
    k = 0;
    curr_err = tol + 1; % To begin iterations.
    err = [];
    while (curr_err > tol)
        
        % j=1 and j=J are the boundary.
        for j=2:J-1
            u_kph_temp = kp1_tridiag\((-2*(1-w)*(1+beta^2)*u(:,j)) ...
                - ((w*beta^2)*(u(:,j+1) + u_new(:,j-1))));
            u_new(2:I-1,j) = u_kph_temp(2:I-1);
        end
        
        curr_err = max(max(abs(u_new - u)));
        
        u = u_new;
        k = k + 1;
        err(k) = curr_err; % Storing a history of errors.
    end
end