clear all
close all
clc

c=1;        % advective speed
L=2*pi;     % computational domain [0,L]
T=2*2*pi;   % end time
M=0;        % intermediate solutions

fexact='exact.dat';

sigma=[0.25, 0.5, 0.75, 1.0, 1.25]; % Courant number
n=25;       % number of interior points

% initial conditions
u0 = @(x) sin(x);  % anonymous function

%method='forward-upwind';
%method='implicit-central';
%method='beam-warming';
%method='lax-wendroff';

f = @(x) sin(x-T);

% solve
for i = 1:size(sigma,2)

    method1='lax-wendroff';
    outLax=wave_solve(c,L,n,sigma(i),T,M,u0,method1);
    method2='beam-warming';
    outBeam=wave_solve(c,L,n,sigma(i),T,M,u0,method2);
    method3='implicit-central';
    outImp=wave_solve(c,L,n,sigma(i),T,M,u0,method3);
    method4='forward-upwind';
    outFor=wave_solve(c,L,n,sigma(i),T,M,u0,method4);

    figure;
    plot(outLax.x,outLax.U(:,2));
    hold on;
    plot(outBeam.x,outBeam.U(:,2));
    hold on;
    plot(outImp.x,outImp.U(:,2));
    hold on;
    plot(outFor.x,outFor.U(:,2));
    hold on;
    fplot(f,[0,L]);
    
    tit = ['Plot of u(x) for Courant Number = ', string(sigma(i))];
    title(tit)
    xlabel('x');
    ylabel('u(x)');
    legend(method1,method2,method3,method4,'Exact');
    xlim([0,L]);

end

% dump
% fout=sprintf('%s_n%g_sigma%f.dat',method,n,sigma);
% dlmwrite(fout,[out.x',out.U],'delimiter',' ','precision','%e');
% dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');

% plot
% xx=linspace(0,L,1000);
% for i=1:size(out.U,2)
%   exact(:,i)=u0(xx-out.TT(i))';
%   plot(out.x,out.U(:,i),'ko-',...
%        xx,u0(xx-out.TT(i)),'r-');
%   axis([0,L,-1.1,1.1]);
%   xlabel('x');
%   ylabel('u(x) and numerical solution');
%   title(sprintf('Time is %f',out.TT(i)));
%   pause
% end

