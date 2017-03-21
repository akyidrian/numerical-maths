function [xvect,xdif,fx,nit]=bisect(a,b,tol,nmax,fun)
%BISECT Bisection method
% [XVECT,XDIF,FX,NIT]=BISECTION(A,B,TOL,NMAX,FUN) tries to find a zero
% of the continuous function FUN in the interval [A,B] using the bisection 
% method. FUN accepts real scalar input x and returns a real scalar value.
% XVECT is the vector of iterates, XDIF the vector of the differences between
% consecutive iterates, FX the residual. TOL specifies the tolerance of the
% method.
err=tol+1; 
nit=0; 
xvect=[]; fx=[]; xdif=[];
while nit<nmax & err>tol
    nit=nit+1; 
    c=(a+b)/2; x=c; fc=fun(x); xvect=[xvect;x]; 
    fx=[fx;fc]; x=a; 
    if fc*fun(x)>0
        a=c; 
    else 
        b=c; 
    end 
    err=0.5*abs(b-a); xdif=[xdif;err];
end
return