clear
clc
close all

%% 1

b = 2;
t = 8;
L = -5;
U = +5;

% numero minimo floating point
x_min = b^(L - 1);
% numero massimo floating point
x_max = b^U * (1 - b^(-t));

x_max/x_min

card = 2 * (b-1) * b^(t-1) * (U-L+1)

%% 2

n = 20;
A = gallery('lehmer', n);
b = ones(n, 1);

R = chol(A);
%   R' * y = b
%   R * x = y
y = R' \ b;
x = R \ y;

R(4, 5)

min(abs(y))

% - Cholesky: 1/3*n^3 + 2*n^2
%   Di cui
%   1/3*n^3 per risolvere R = chol(A)
%   n^2 per risolvere R'y=b
%   n^2 per risolvere Rx=y

nflops = 1/3*n^3 + 2*n^2

%% 3

% % Residuo
% res = b - A * x_it;
% % Residuo relativo/normalizzato
% res_rel = norm(res) / norm(b);
% res_norm = norm(b - A * x_it) / norm(b);
% % Stima (maggiorazione) errore relativo (quando deltaA ~ 0, wilkinson-higham)
% stima_err_rel = res_norm * cond(A);

b_norm = 10;
condA = 1e8;
res_rel = 1e-3 / b_norm;
stima_err_rel = res_rel * condA

%% 4

n = 10;
A = gallery('lehmer', n);
nmax = 100

tol_1 = 1e-1;
D_1 = qrbasic(A,tol_1,nmax);
K_A_1 = max(abs(D_1))/min(abs(D_1))

tol_2 = 1e-3;
D_2 = qrbasic(A,tol_2,nmax);
K_A_2 = max(abs(D_2))/min(abs(D_2))

K_A = max(abs(eig(A)))/min(abs(eig(A)))

%% 5

norm([1, sqrt(3)-1/2])

%% 6

alpha = 0;
f = @(x) exp(3 * x.^2) - 1;
df = @(x) 6*exp(3 * x.^2) * x;
d2f = @(x) 6*exp(3 * x^2) * (6*x.^2 + 1);

a = -1;
b = 1;
x = a:0.01:b;
plot(x, f(x))

toll = 1e-4;
% [xvect, it] = bisez(a,b,toll,f)

df(alpha)
d2f(alpha)

mol = 1;
nmax = 100;
x0 = 0.5;
[xvect, it] = newton(x0,nmax,toll,f,df,mol);
[p,c] = stimap(xvect);

% lim_k->inf (x_k+1 - alpha)/(x_k - alpha)^2 = 1/2 * d2phi_N(alpha) = 1/2 * d2f(alpha)/df(alpha)
%   mu = 1/2 * d2f(alpha)/df(alpha)

%% 7

% (x_10 - alpha)/(x_k - alpha)^p = 1/p! * d(p)phi(alpha)

%% 8

Phi = @(y1, y2) 0.5 * y1.^2 + (1 - y1/2) .* y2.^2;
GradPhi = @(y1, y2) [y1 - y2.^2; ...
                    -y1.*y2 + 2*y2];
HessPhi = @(y1, y2) [1, -y2; ...
                    -y2, 2*(1-y1/2)];
x_ex = [0, 0]';

x0 = [1/4, 1/2]';
toll = 1e-8;
nmax = 2;
[xvect,it] = newton_opt(GradPhi,HessPhi,x0,toll,nmax);
xvect(:, 3)