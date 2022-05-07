function [xm,U,Q,W]=integrated(fh,A,B,N)

%
% function [xm,U,Q,W]=integrated(fh,A,B,N)
%
% returns the value of function with function handle fh
% integrated between A and B with N intervals
%
% A,B: extremes of the domain
% N: number of intervals, each of width h = (B-A) / N
% xm: mid point for each interval
% U: average in each
% Q: plot bars
% W: plot a staircase
%

h=(B-A)/N;
a=linspace(A,B-h,N)';
b=linspace(A+h,B,N)';
xm=(a+b)*0.5;

for i=1:N
  U(i)=quad(fh,a(i),b(i),1e-10)/h;
end
U=U';

% organize matrix
Q=zeros(N+1,N+1);
Q=[a;B];
for i=1:N
  Q(:,i+1)=NaN;
  Q(i:i+1,i+1)=U(i);
end

% new matrix
l=1;
for i=1:N
  W(l:l+1,1)=[a(i),b(i)]';
  W(l:l+1,2)=U(i);
  l=l+2;
end
