function out=wave_solve(c,L,n,sigma,T,M,u0,method)

%
% --inputs--
% c:        advective speed
% L:        domain size [0,L]
% n:        number of interior grid points
% sigma:    Courant number
% T:        final time
% M:        number of solutions recorded between [0,T]
% u0:       function that prescribes the initial conditions u0(x)
% method:   integration method, one of:
%           'forward-upwind' 
%           'implicit-central'
%           'beam-warming'
%           'lax-wendroff' 
%         ...
%
% --outputs--
% out.h   grid spacing
% out.k   time step size
% out.l   number of time steps taken
% out.x:  spatial locations so that out.x(1)=0 and out.x(end)=L
% out.TT: out.TT(1)=0 and out.TT(end)=T with
%         length(out.TT)=M+2;
% out.U:  numerical solution as matrix
%         out.U(:,j) is the numerical solution at time out.TT(j)
%         with j=1,\dots,M+2
%         size(out.U,1)=length(out.x)
%         size(out.U,2)=length(out.TT)
%

% set output to empty
out=[];

% work on grid
h=L/(n+1); % grid spacing recovered from the number of interior points
out.h = h; % store it

out.x=[0:h:L]; % actual grid array, including x=0 and x=L
N=length(out.x); % number of overall points

% time outputs
out.TT=linspace(0,T,M+2);

% build the matrix for the updates
switch lower(method)

case 'forward-upwind'

  if ( c<0 ) error('please specify a positive advective speed'); end;

  A = -diag(ones(N,1),0) - ...
      -diag(ones(N-1,1),-1);
  A(1,n+1)=1; % periodic boundary on U(1)=U_0

case 'implicit-central'
  % nothing to do

case 'beam-warming'
  % nothing to do

case 'lax-wendroff'
  % nothing to do

otherwise
  error('method is unknown');
end

% time step size recovered from Courant number
k=sigma*h/c;
out.k = k; % store it

% initial conditions
U_=u0(out.x)'; t=0; j=1;

% store initial conditions
out.U(:,j)=U_; j=j+1;

% integrate in time
l=0;
while t<out.TT(end)

   % pick the smallest between the time step
   % and the time step to get to the next out.TT(j)
   k_=min([k,(out.TT(j)-t)]);
   sigma_=k_*c/h;

   % fprintf('Time: %f; Sigma = %f; Time step = %f\n',t,sigma_,k_);
   
   % zero the update
   dU_=zeros(size(U_));

   switch lower(method)

   case {'forward-upwind'}

      dU_=sigma_*A*U_; % Euler fwd step

   case {'implicit-central'}

      % Instead of solving a system where Ax = b with A and x being known,
      % implicit methods instead solve with A and b being known. Therefore,
      % we construct the matrix A and its appropriate boundary conditions
      % with that in mind:

      % fill in...
      A = diag(ones(N,1),0) - ...
      (sigma_/2)*diag(ones(N-1,1),-1) - ...
      -(sigma_/2)*diag(ones(N-1,1),1);
      %A(1,n+1)=1; % periodic boundary on U(1)=U_0

      % Boundary Conditions
      A(1,N-1)=-(1/2)*sigma_;
      A(N,2)=(1/2)*sigma_;
      
      dU_ = A \ U_ - U_;

   case 'beam-warming'

      % fill in...
      % Creates a matrix that multiplies by U_ to get the next U_ in time.
      % When the boundaries are reached, utilize the value on the opposite
      % side of the matrix due to the periodic nature of the solution.

     temp = 3*diag(ones(N,1),0) + ...
      -4*diag(ones(N-1,1),-1) + ...
      diag(ones(N-2,1),-2);

     temp = (-sigma/2) * temp;

     temp = temp + (sigma^2/2) * (diag(ones(N,1),0) + ...
            -2*diag(ones(N-1,1),-1) + ...
            diag(ones(N-2,1),-2));

     A = zeros(N);
     A = temp;

     % Boundary Conditions
     A(1,1) = -(sigma/2)*3 + (sigma^2/2);
     A(1,N-2:N-1) = [-(sigma/2) + (sigma^2/2), (sigma/2)*4 + (sigma^2/2)*-2];
     A(2,1:2) = [(sigma/2)*4 + (sigma^2/2)*-2, -(sigma/2)*3 + (sigma^2/2)];
     A(2,N-1) = -(sigma/2) + (sigma^2/2);

     dU_ = A * U_;
   
   case 'lax-wendroff'
      
      % fill in...
      % Creates a matrix that multiplies by U_ to get the next U_ in time.
      % When the boundaries are reached, utilize the value on the opposite
      % side of the matrix due to the periodic nature of the solution.


     temp = diag(ones(N-1,1),1) + ...
            -diag(ones(N-1,1),-1);

     temp = (-sigma/2) * temp;

     temp = temp + (sigma^2/2) * (-2*diag(ones(N,1),0) + ...
            diag(ones(N-1,1),1) + ...
            diag(ones(N-1,1),-1));

     A = temp;

     % Boundary Conditions
     A(1,N-1) = (sigma/2) + (sigma^2/2);
     A(N,2) = (-sigma/2) + (sigma^2/2);
     dU_ = A * U_;

   otherwise
     error('method is unknown');
   end

   % update
   U_=U_ + dU_;

   % advance time to reflect update
   t=t+k_;
   l=l+1; % step counter

   % store
   if ( t==out.TT(j) )
     out.TT(j)=t; % adjust recorded time
     out.U(:,j)=U_; j=j+1;
   end

end

out.l=l; % number of steps

