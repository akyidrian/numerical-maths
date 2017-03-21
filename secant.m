function [xvect,xdif,fx,nit]=secant(xm1,x0,tol,nmax,fun)
%SECANT Secant method
% [XVECT,XDIF,FX,NIT]=SECANT(XM1,X0,TOL,NMAX,FUN) tries to find a zero
% of the continuous function FUN  using the secant method. FUN accepts 
% real scalar input x and returns a real scalar value. XVECT is the 
% vector of iterates, XDIF the vector of the differences between
% consecutive iterates, FX the residual. TOL specifies the tolerance of the
% method.
x=xm1; fxm1=fun(x); 
xvect=[x]; fx=[fxm1]; 
x=x0; fx0=fun(x);  
xvect=[xvect;x]; fx=[fx;fx0]; 
err=tol+1; nit=0; xdif=[];
while nit<nmax & err>tol
    nit=nit+1; 
    x=x0-fx0*(x0-xm1)/(fx0-fxm1); 
    xvect=[xvect;x]; 
    fnew=fun(x); 
    fx=[fx;fnew]; 
    err=abs(x0-x); 
    xdif=[xdif;err];
    xm1=x0; fxm1=fx0; 
    x0=x; fx0=fnew;
end
return