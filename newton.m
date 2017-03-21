function [xvect,xdif,fx,nit]=newton(x0,tol,nmax,fun,dfun)
%NEWTON Newton method
% [XVECT,XDIF,FX,NIT]=NEWTON(XM1,X0,TOL,NMAX,FUN,DFUN) tries to find a zero
% of the continuous function FUN using the Newton method starting from
% the initial guess X0. FUN and DFUN accept
% real scalar input x and returns a real scalar value. XVECT is the 
% vector of iterates, XDIF the vector of the differences between
% consecutive iterates, FX the residual. TOL specifies the tolerance of the
% method.
err=tol+1; nit=0; xvect=x0; x=x0; fx=fun(x); xdif=[];
while nit<nmax & err>tol
    nit=nit+1; 
    x=xvect(nit); 
    dfx=dfun(x);
    if dfx==0
        err=tol*1.e-10; %WTF?
        fprintf('Stop for vanishing dfun\n');
%         return
    else
        xn=x-fx(nit)/dfx; err=abs(xn-x); xdif=[xdif; err];
        x=xn; xvect=[xvect;x]; fx=[fx;fun(x)];
    end
end
return