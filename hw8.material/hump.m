function [u,du]=hump(x)

%
% function [u,du]=hump(x)
%
% x: independent variable 
% u: function value at x
% du: du/dx at x
%

lambda=1; a=0.5; b=1;
uL=1; uR=1;

% u(x)
u=a*(cos(2*pi*(x./lambda))+1)+b;
u(x>0.5*lambda)=uR; u(x<-0.5*lambda)=uL;

% u'(x)
du=-(2*a*pi*sin((2*pi*x)/lambda))/lambda;
du(x>0.5*lambda)=0; du(x<-0.5*lambda)=0;
