%% Homework 6 Amplification Factor
% Constants
s = 0.5;
h = 0.1;
NCount = [25, 50, 100, 200, 400];
beta = (2*pi)./N;

% Forward-Upwind
GFou = 1-s+s*cos(h)-i*s*sin(h);
AmpFou = abs(GFou);
NFou = log(0.5)/log(AmpFou);
% Val = 5.542866375968506e+02

% Lax-Wendroff
GLax = 1-s^2+s^2*cos(h)-i*s*sin(h);
AmpLax = abs(GLax);
NLax = log(0.5)/log(AmpLax);
% Val = 2.962354606765640e+05

% Determine h Size
G = 0.999997660150584;
GFouF = @(x) abs(1-s+s*cos(x)-i*s*sin(x))-G;
hVal = -1*fzero(GFouF,0);
figure;
fplot(GFouF,[-pi, pi]);
xlabel('h');
ylabel('Amplification Factor')
% Val = 0.004326523156974