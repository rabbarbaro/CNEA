clear
clc
close all

%% T1

beta = 2;
t = 3;
L = -4;
U = 4;

x_max = beta^U * (1-beta^(-t))

m = bin2dec("110");
s = 0;
x = (-1)^s * m * beta^(U - t)

%% T2

x = 1;

for n = 1:10
    x = (x + 107/x) / 2;
end

% 1 somma + 1 divisione + 1 divisione
% 3 * 10 = 30

%% T3

% A SDP --> cholesky
n = 60;
n_sys = 20;

% la fattorizzazione Cholesky è O(1/3 * n^3) (una volta)
% i sistemi lineari sono O(2n^2)

nflops = 1/3 * n^3 + n_sys * 2 * n^2

%% T4

A = [4 6 2 1
    2 3 -1 1
    2 1 2 1
    1 1 1 1];
b = (1:4)';

[L, U, P] = lu(A);
y = fwsub (L, P * b);
x = bksub (U, y);
x(4)

P
% pivoting 2/3

%% T5

% A quadrata invertibile
% 1. F non per forza
% 2. F di B
% 3. V
% 4. V
% 5. F non sempre

%% T6

A = [3 1
    2 8];

D = qrbasic(A,1e-10,3)

%% T7

f = @(x) log(-x.^2 + 4*x - 3) * sin(pi * x);
df = @(x) 1/(-x.^2 + 4*x - 3) * (-2*x + 4) * sin(pi * x) + log(-x.^2 + 4*x - 3) * pi * cos(pi * x);
alpha = 2;

[xvect,it] = newton(1.5,100,1e-6,f,df,1);
[xvect,it] = newton(1.5,100,1e-6,f,df,2);
[xvect,it] = newton(1.5,100,1e-6,f,df,3);
[xvect,it] = newton(1.5,100,1e-6,f,df,4);

% 3

%% T8

f = @(x) 2 * log(x.^2 + 1) - 1;
x0 = 0;
x = 2;

% implemento metodo secanti
for k = 1:2
    q = (f(x) - f(x0)) / (x - x0);
    x0 = x;
    x = x - f(x) / q;
end

x

%% T9

phi = @(x) exp(-x.^2 / 2) - sin(x) + x;

x = 0:0.01:2;
plot(x, phi(x), x, x)

x0 = 0;
nmax = 100;
toll = 1e-8;
[succ, it] = ptofis(x0, phi, nmax, toll);

%% T10

phi = @(x) 1/2 * exp(2*x - 1) - x + 1/2;
dphi = @(x) exp(2*x - 1) - 1;
d2phi = @(x) 2 * exp(2*x - 1);

x = 0:0.01:1;
plot(x, phi(x), x, dphi(x), x, d2phi(x), x, x, x, 0*x)
legend('phi', 'dphi', 'd2phi', 'x')

x0 = 1;
nmax = 3;
toll = 1e-8;
[succ, it] = ptofis(x0, phi, nmax, toll)

% p = 2

