function u=step(x)

% u(x)
u=zeros(size(x));
u(x>0)=1; u(x<0)=2;
