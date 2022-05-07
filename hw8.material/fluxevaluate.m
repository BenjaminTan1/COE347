function F=fluxevaluate(fh,kh,U,uL,uR,method)

%
% function F=fluxevaluate(fh,kh,U,uL,uR,method)
%
% fh: function handle to the flux function
% kh: k/h factor
% U:  solution at time t = t_n
% uL: left boundary condition u(-L)=uL
% uR: right boundary condotion u(+L)=uR
% method: one of 'first-order-upwind', 'law-wendroff',
%                'richtmyer', 'maccormack'

N=length(U);
F=zeros(N+1,1);

%
%     
%    F(1)  F(2)   F(3)       F(i-1)  F(i)   F(i+1)            F(N)   F(N+1)
%    |------|------|--- //  ---|-------|------|------|--- // ---|-----|
%       1      2                  i-1      i     i+1               N
%   -L                                                                L   
%

switch lower(method)
case 'first-order-upwind'

  % F(i) = F^L_i, the left numerical flux
  for i=2:N
     if ( U(i)>0 ) F(i)=feval(fh,U(i-1)); end;
     if ( U(i)<0 ) F(i)=feval(fh,U(i)); end;
  end

  % BCs

  % left
  if ( U(1)>0 ) F(1)=feval(fh,uL); end
  if ( U(1)<0 ) F(1)=feval(fh,U(1)); end

  % right
  if ( U(N)>0 ) F(N+1)=feval(fh,U(N)); end
  if ( U(N)<0 ) F(N+1)=feval(fh,uR); end
  

case 'lax-wendroff'
   F = zeros(N,1);
   % code here
   for i = 2:N-1
        if (U(i)>0) 
            F1 = feval(fh,U(i+1))-feval(fh,U(i-1));
            F2 = -kh*((U(i+1)+U(i))/2)*(feval(fh,U(i+1))-feval(fh,U(i)));
            F3 = kh*((U(i-1)+U(i)/2)*(feval(fh,U(i))-feval(fh,U(i-1))));
            F(i) = F1 + F2 + F3;
        end
        if (U(i)<0)
            F1 = feval(fh,U(i+2))-feval(fh,U(i));
            F2 = -kh*((U(i+2)+U(i+1))/2)*(feval(fh,U(i+2))-feval(fh,U(i+1)));
            F3 = kh*((U(i-1)+U(i)/2)*(feval(fh,U(i+1))-feval(fh,U(i))));
            F(i) = F1 + F2 + F3;
        end
   end
   % BCs

   % left
   if ( U(1)>0 )
       F1 = feval(fh,U(2))-feval(fh,uL);
       F2 = kh*((U(2)+U(1))/2)*(feval(fh,U(2))-feval(fh,U(1)));
       F3 = -kh*((uL+U(1)/2)*(feval(fh,U(1))-feval(fh,uL)));
       F(1) = F1 + F2 + F3; 
   end
   if ( U(1)<0 )
       F1 = feval(fh,U(3))-feval(fh,U(1));
       F2 = kh*((U(3)+U(2))/2)*(feval(fh,U(3))-feval(fh,U(2)));
       F3 = -kh*((uL+U(2)/2)*(feval(fh,U(2))-feval(fh,U(1))));
       F(1) = F1 + F2 + F3; 
   end
   % right
   if ( U(N)>0 ) 
       F1 = feval(fh,U(N))-feval(fh,U(N-2));
       F2 = kh*((uR+U(N-1))/2)*(feval(fh,U(N))-feval(fh,U(N-1)));
       F3 = -kh*((U(N-2)+U(N-1)/2)*(feval(fh,U(N-1))-feval(fh,U(N-2))));
       F(N) = F1 + F2 + F3; 
   end
   if ( U(N)<0 )
       F1 = feval(fh,uR)-feval(fh,U(N-1));
       F2 = kh*((uR+U(N))/2)*(feval(fh,uR)-feval(fh,U(N)));
       F3 = -kh*((U(N-1)+U(N)/2)*(feval(fh,U(N))-feval(fh,U(N-1))));
       F(N) = F1 + F2 + F3; 
   end

case 'richtmyer'

   % code here
   for i = 2:N
        if (U(i)>0) 
            FMid = (U(i-1) + U(i))/2 - (kh/2) * (feval(fh,U(i)) - feval(fh,U(i-1)));
            F(i) = feval(fh, FMid); 
        end
        if (U(i)<0)
            FMid = (U(i) + U(i+1))/2 - (kh/2) * (feval(fh,U(i+1)) - feval(fh,U(i)));
            F(i) = feval(fh, FMid); 
        end
    end

    % BCs
    % left
    if ( U(1)>0 ) 
        FMid = U(1) - (kh/2) * (feval(fh, U(1)) - feval(fh, uL));
        F(1)=feval(fh,FMid); 
    end
    if ( U(1)<0 ) 
        FMid = U(1) - (kh/2) * (feval(fh, U(2)) - feval(fh, U(1)));
        F(1)=feval(fh,FMid); 
    end

    % right
    if ( U(N)>0 ) 
        FMid = U(N-1) - (kh/2) * (feval(fh, U(N)) - feval(fh, U(N-1)));
        F(N+1)=feval(fh,FMid); 
    end
    if ( U(N)<0 )
        FMid = U(N) - (kh/2) * (feval(fh, uR) - feval(fh, U(N)));
        F(N+1)=(1/2)*feval(fh,FMid);
    end

case 'maccormack'

    % code here
    for i = 2:N
        if (U(i)>0) 
            UProv = U(i) - kh * (feval(fh, U(i))-feval(fh, U(i-1))); 
            F(i) = feval(fh, UProv) + (1/kh)*(1/4)*UProv - (1/kh)*(1/4)*U(i); 
        end
        if (U(i)<0)
            UProv = U(i) - kh * (feval(fh,U(i+1)) - feval(fh,U(i)));
            F(i) = feval(fh, UProv) + (1/kh)*(1/4)*UProv - (1/kh)*(1/4)*U(i);
        end
    end

    % BCs
    % left
    if (U(i)>0) 
        UProv = uL - kh * (feval(fh, U(1))-feval(fh, uL)); 
        F(1) = feval(fh, UProv) + (1/kh)*(1/4)*UProv - (1/kh)*(1/4)*U(1); 
    end
    if (U(i)<0)
        UProv = U(1) - kh * (feval(fh,U(2)) - feval(fh,U(1)));
        F(1) = feval(fh, UProv) + (1/kh)*(1/4)*UProv - (1/kh)*(1/4)*U(1);
    end

    % right
    if (U(N)>0) 
        UProv = U(N) - kh * (feval(fh, U(N))-feval(fh, U(N-1))); 
        F(N+1) = feval(fh, UProv) + (1/kh)*(1/4)*UProv - (1/kh)*(1/4)*U(1); 
    end

    if (U(N)<0)
        UProv = uR - kh * (feval(fh,uR) - feval(fh,U(N)));
        F(N+1) = feval(fh, UProv) + (1/kh)*(1/4)*UProv - (1/kh)*(1/4)*U(1);
    end

otherwise
  error('method is unknown');

end


