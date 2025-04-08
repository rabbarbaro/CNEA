clear
clc
close all

%% T1

b = 2;
t = 5;
L = -4;
U = 4;

% considerando i negativi, escludendo lo zero
card = 2 * (b-1) * b^(t-1) * (U-L+1)

s = 0;
m = bin2dec('11001');
e = -1;

x = (-1)^s * m * b^(e-t)

%% T2

n = 0;
S_n = 1;
while abs(exp(1) - S_n) > 1e-4 && n < 5
    n = n + 1;
    S_n = S_n + 1/factorial(n)
end
n

% con 5 iterazioni 2.7167
% 7 iterazioni per 1e-4

%% T3

n = 20;
% tridiagonale, uso thomas
% mi serve solo la diagonale di U (beta --> alpha)
% (n - 1)
% 2 * (n - 1)
% moltiplico tutti questi elementi
% n - 1 moltiplicazioni
nflops = 4 * (n - 1)

%% T4

% 1 V
% 2 F, devono esserlo tutti i minori principali
2/3 * 60^3
% 3 V
% 4 F
% 5 F (fill-in)

%% T5

A = [4 3 2
    3 4 1
    2 1 4
    0 1 2];
b = ones(4, 1);

[Q,R] = qr(A);
x = R \ (Q' * b);

Q(1, 2)
x(3)

%% T6

n = 100;
A = 5*diag(ones(n,1)) - 2*diag(ones(n-1,1),1) -2*diag(ones(n-1,1),-1);
s = 6;
tol = 1e-8;
nmax = 200;
x0 = ones(n, 1);

[lambda,x,iter] = invpowershift(A,s,tol,nmax,x0);

%% T7

% 1 F
% 2 V ?
% 3 V
% 4 F
% 5 V

%% T8

f = @(x) (x - 1) .* sin(x);
df = @(x) sin(x) + (x - 1) .* cos(x);
d2f = @(x) cos(x) + cos(x) - (x - 1) .* sin(x);
a = 0;
b = 1;
x = a:0.01:b;
plot(x, f(x), x, df(x), x, 0*x)
legend('f', 'df')

% graficamente molteplicità 1
[xvect,it] = newton(1,5,1e-12,df,d2f, 1)

f(xvect(5))

%% T9

f = @(x) sin(2 * pi * x);
g = @(x) pi * (x - 1/6) + sqrt(3)/2;
h = @(x) f(x) - g(x);
dh = @(x) pi * (- 1 + 2*cos(2 * pi * x));

x0 = 0;
tol = 1e-4;
nmax = 100;

[xvect,it] = newton(x0,nmax,tol,h,dh,1);
[p,c] = stimap(xvect);

xvect(end)
f(xvect(end))

[xvect,it] = newton(x0,nmax,tol,h,dh,2);
[p,c] = stimap(xvect);

%% T10

phi = @(x) eta * (1 - exp(2*x - 1)) + x;
dphi = -2 * eta * exp(2*x - 1) + 1;
alpha = 0.5;

% teoria: |phi'(alpha)| < 1
% |-2 * eta * exp(2*alpha - 1) + 1| < 1
% |-2 * eta * exp(2*1/2 - 1) + 1| < 1
% |-2 * eta + 1| < 1
% eta in (0, 1)

% teoria: phi'(alpha) > 0 (da mettere a sistema con la condizione di
% convergenza)
% -2 * eta + 1 > 0
% -2 * eta > -1
% eta < 1/2
% a sistema con eta in (0, 1)
% eta in (0, 1/2)