function [xm,U_]=advanceconservative(uh,fh,uL,uR,L,sigma,N,T,method)

%
% function [xm,U_]=advanceconservative(uh,fh,uL,uR,L,sigma,N,T,method)
%
% --inputs--
% uh: handle to initial conditions
% fh: handle to f(u) and df/du, i.e. [f,df]=feval(fh,u)
% uL: left state
% UR: right state
% L:  domain size so that x \in [-L,L]
% sigma: Courant number
% N: number of cells
% T: final time
% method: one of 'first-order-upwind', 'law-wendroff',
%                'richtmyer', 'maccormack'
%
% --outputs--
% xm: location of midpoints
% U_: numerical solution
%

% initial conditions
[xm,U_]=integrated(uh,-L,L,N);
h=xm(2)-xm(1);

if ( ~iscolumn(xm) ) xm=xm'; end;
if ( ~iscolumn(U_) ) U_=U_'; end;

%
%     
%    F(1)  F(2)   F(3)       F(i-1)  F(i)   F(i+1)            F(N)   F(N+1)
%    |------|------|--- //  ---|-------|------|------|--- // ---|-----|
%       1      2                  i-1      i     i+1               N
%   -L                                                                L   
%

t=0; l=0;
if method == "lax-wendroff"
    while t<T

        % find max(abs(lambda))
        for i=1:N
            [dummy,l]=feval(fh,U_(i)); abslambda_(i)=abs(l);
        end

        % pick the smallest between the time step
        % and the final time
        k=min([h*sigma/max(abslambda_),(T-t)]);

        % fluxes 
        F_=fluxevaluate(fh,k/h,U_,uL,uR,method);

        % apply them
        U_ = U_ - (1/2)*(k/h)*F_;
        t = t+k; l=l+1;

        %fprintf('%5.5d, t=%e, h=%e\n',l,t,k);
    end
else
    while t<T

        % find max(abs(lambda))
        for i=1:N
            [dummy,l]=feval(fh,U_(i)); abslambda_(i)=abs(l);
        end

        % pick the smallest between the time step
        % and the final time
        k=min([h*sigma/max(abslambda_),(T-t)]);

        % fluxes 
        F_=fluxevaluate(fh,k/h,U_,uL,uR,method);

        % apply them
        dF_ = F_(2:end)-F_(1:end-1);
        U_ = U_ - (k/h)*dF_;
        t = t+k; l=l+1;

        %fprintf('%5.5d, t=%e, h=%e\n',l,t,k);
    end
end
