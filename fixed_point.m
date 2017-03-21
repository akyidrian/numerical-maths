function [xvect, xdif, nit] = fixed_point(x0, tol, nmax, fun)
% Runs Fixed Point method.
    nit = 0;
    xvect = x0;
    xdif = [];
    err = tol + 1;
    while (nit<nmax) && (err>tol)
        nit = nit + 1;
        x=xvect(nit);
        x = fun(x);
        xvect = [xvect;x];

        err = abs(x - x0); xdif=[xdif;err];
        x0 = x;
    end
end