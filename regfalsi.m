function [xvect,xdif,fx,nit]=regfalsi(xm1,x0,tol,nmax,fun)
%REGFALSI Regula Falsi method
% [XVECT,XDIF,FX,NIT]=REGFALSI(XM1,X0,TOL,NMAX,FUN) tries to find a zero
% of the continuous function FUN in the interval [XM1,X0] using the Regula Falsi
% method. FUN accepts real scalar input x and returns a real scalar value.
% XVECT is the vector of iterates, XDIF the vector of the differences between
% consecutive iterates, FX the residual. TOL specifies the tolerance of the
% method.
nit=0; 
x=xm1; f=fun(x); fx=[f]; 
x=x0; f=fun(x); fx=[fx, f]; 
xvect=[xm1,x0]; xdif=[]; f=tol+1; kprime=1;
while nit<nmax & abs(f)>tol
    nit=nit+1; 
    dim=length(xvect); 
    x=xvect(dim); 
    fxk=fun(x); 
    xk=x; i=dim; 
    while i>=kprime
        i=i-1; x=xvect(i); fxkpr=fun(x);
        if fxkpr*fxk<0 
            xkpr=x; kprime=i; break; 
        end 
    end
    x=xk-fxk*(xk-xkpr)/(fxk-fxkpr); 
    xvect=[xvect, x]; f=fun(x); 
    fx=[fx, f]; err=abs(x-xkpr); xdif=[xdif, err];
end
return