clear
clc
close all

%% T1

s = 0;
e = 3;
m = 21;
b = 2;
t = 5;

x = (-1)^s * m * b^(e - t)

%% T2

% divisione 5/x_n
% somma x_n + (5/x_n)
% divisione per 2

nflops = 50 * 3

%% T3

% fattorizzazione A 3 * (n - 1)
% risoluzione sistema lineare 2 * (n - 1) + 3 * (n - 1) + 1

n = 90;
jj = 10;

flops = 3 * (n - 1) + jj * (2 * (n - 1) + 3 * (n - 1) + 1)

%% T4

% [3/(2 * sqrt(2)) + gamma; 2]

%% T5

A = [4 -1 0 0 0
    1 4 0 0 0
    0 -1 3 0 0
    0 0 1 6 0
    0 0 0 9 -1];

sort(eig(A), 'descend')

% dato che gli autovalori lambda_2 e lambda_3 sono c.c. (stesso modulo) il
% metodo delle potenze inverse non ci converge
% lo shift deve stare più vicino a 3 che a 6 o -1 (e diverso da 3)
% 1 < s < 4, s != 3

%% T6

a = 0;
b = 5;
tol = 1e-3;

k_min = ceil(log2((b - a)/tol) - 1)

%% T7

phi = @(x) 0.25*(x.^4) + 0.5*(x.^2) - 3*x + 5;
dphi = @(x) x.^3 + x - 3;
d2phi = @(x) 3*x.^2 + 1;
x0 = 0;
nmax = 5;
tol = 1e-6;

[xvect,it] = newton(x0, nmax, tol, dphi, d2phi);
xvect(2)
xvect(3)
xvect(6)

%% T8

f = @(x) sin(pi/3 * x) * (x - 3).^2;
df = @(x) pi/3 * cos(pi/3 * x) * (x - 3).^2 + 2*(sin(pi/3 * x) * (x - 3));
df(3)
% ...
% molteplicità 3

x0 = 4;
nmax = 2;
tol = 1e-12;
[xvect,it] = newton(x0, nmax, tol, f, df, 3);
xvect(2)
xvect(3)

%% T9

phi = @(x) x - 1/3 * (1 - exp(1 - 3*x));
dphi = @(x) 1 - exp(1 - 3*x);
dphi(1/3)
d2phi = @(x) 3*exp(1 - 3*x);
d2phi(1/3)

%% T10

% | 6gamma | < 1 E 6gamma <= 0
% - 1/6 < gamma < 0