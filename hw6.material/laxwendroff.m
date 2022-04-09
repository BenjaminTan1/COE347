function [mG,pG,G]=laxwendroff(x,s)
%
% lax-friedrichs 
%

G=1-s^2+s^2*cos(x)-i*s*sin(x);
Ge=exp(-i*s*x);
mG=abs(G);
pG=angle(G)./angle(Ge);
