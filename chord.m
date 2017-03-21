function [xvect,xdif,fx,nit]=chord(a,b,x0,tol,nmax,fun)
%CHORD Chord method
% [XVECT,XDIF,FX,NIT]=CHORD(A,B,TOL,NMAX,FUN) tries to find a zero
% of the continuous function FUN in the interval [A,B] using the chord 
% method. FUN accepts real scalar input x and returns a real scalar value.
% XVECT is the vector of iterates, XDIF the vector of the differences between
% consecutive iterates, FX the residual. TOL specifies the tolerance of the
% method.
x=a; fa=fun(x); 
x=b; fb=fun(x); 
r=(fb-fa)/(b-a);
err=tol+1; nit=0; xvect=x0; x=x0; fx=fun(x); xdif=[];
while nit<nmax & err>tol
    nit=nit+1; 
    x=xvect(nit); 
    xn=x-fx(nit)/r; 
    err=abs(xn-x); 
    xdif=[xdif; err]; 
    x=xn; 
    xvect=[xvect;x]; 
    fx=[fx;fun(x)];
end
return