clear all
close all
clc

c=1;        % advective speed
L=4*pi;     % computational domain [0,L]
T=2*2*pi;   % end time
M=0;        % intermediate solutions

fexact='exact.dat';

sigma = [0.25, 0.50, 0.75]; % Courant number
n = [25,50,100,200,400];       % number of interior points

%method='forward-upwind';
%method='implicit-central';
%method='beam-warming';
%method='lax-wendroff';

% initial conditions
u0 = @(x) sin(x);  % anonymous function
f = @(x) sin(x-T);

% solve
%out=wave_solve(c,L,n,sigma,T,M,u0,method);

method = ["forward-upwind", "lax-wendroff"];

% solve
e = zeros(size(n,2),size(method,2)*size(sigma,2));
for i = 1:size(sigma,2)
    for j = 1:size(method,2)
        figure;
        for k = 1:size(n,2)
            out=wave_solve(c,L,n(k),sigma(i),T,M,u0,method(j));

            e(k,size(method,2)*(i-1)+j) = errorCalc(out,f);

            plot(out.x,out.U(:,2));
            hold on;
            fplot(f,[0,L]);
            hold on;
        end
        legendCell = cellstr(num2str(n', 'n=%-d'));
        legendCell{end+1} = 'Exact';
        tit = ['Plot of u(x) for Courant Number = ', string(sigma(i)), ' using ', method(j)];
        title(tit);
        xlabel('x');
        ylabel('u(x)');
        legend(legendCell);
        xlim([0,L]);
    end
end

%% Plotting Error
% Creating Legend
sigma2 = zeros(1, 2*size(sigma,2));
leg = strings(1, 2*size(sigma,2));
for i = 1:size(sigma,2)
    sigma2(2*i-1) = sigma(i);
    sigma2(2*i) = sigma(i);
    leg(1,2*i-1) = ['Method = ' + method(1) + ', Sigma = ' + sigma2(2*i-1)];
    leg(1,2*i) = ['Method = ' + method(2) + ', Sigma = ' + sigma2(2*i-1)];
end

% Plotting
h = 1./(n+1);
figure;
for i = 1:size(e,2)
    loglog(h, e(:,i));
    hold on;
end
legendCell = cellstr(leg);
legend(legendCell);
title('Error for Different Methods at Different Meshes and Courant Numbers');
xlabel('h');
ylabel('error');

% Find C and alpha
eLog = log(e);
nLog = log(h);
p = linearize(nLog, eLog);

% Plot Test
figure;
for i = 1:size(p,2)
    pFunc = @(x) p(1,i)*x + p(2,i);
    fplot(pFunc,[nLog(end) nLog(1)]);
    hold on;
end
legend(legendCell);
title('Error for Different Methods at Different Meshes and Courant Numbers');
xlabel('log(h)');
ylabel('log(error)');

function [e] = errorCalc(out, f)
    sum = 0;
    for i = 1:size(out.x,2)
        sum = sum + (f(out.x(i)) - out.U(i,2))^2;
    end
    e = (1/size(out.x,2))*sum^(1/2);
end

function [pF] = linearize(n, e)
    p = zeros(2, size(e,2));
    for i = 1:size(e,2)
        p(:,i) = polyfit(n, e(:,i),1);
    end
    pF = p;
end

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

% dump
% fout=sprintf('%s_n%g_sigma%f.dat',method,n,sigma);
% dlmwrite(fout,[out.x',out.U],'delimiter',' ','precision','%e');
% dlmwrite(fexact,[xx',exact],'delimiter',' ','precision','%e');
% 
