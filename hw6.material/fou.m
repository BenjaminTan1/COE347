function [mG,pG,G]=fou(x,s)
%
% lax-friedrichs 
%

G=1-s+s*cos(x)-i*s*sin(x);
Ge=exp(-i*s*x);
mG=abs(G);
pG=angle(G)./angle(Ge);
