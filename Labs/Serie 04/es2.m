clear
clc
close all

%% 1

f = @(x) x.^3 - (2 + exp(1))*x.^2 + (2*exp(1) + 1)*x + (1 - exp(1)) - cosh(x - 1);
df = @(x) 3*x.^2 - 2*(2 + exp(1))*x + (2*exp(1) + 1) - sinh(x - 1);

x = 0.5:0.01:6.5;

plot(x, f(x), x, df(x))
grid on
hold on
plot(x, 0*x, '--')

d2f = @(x) 6*x - 2*(2 + exp(1)) - cosh(x - 1);
d2f(1)

% molteplicictà 2

%% 2

%% 3

mol = 2;
nmax = 1000;
toll = 1e-6;
x0 = .5;

[xvect_1, it_1] = newton(x0, nmax, toll, f, df)
[xvect_1_m, it_1_m] = newton(x0, nmax, toll, f, df, mol)

%%

x0 = 3;
[xvect_2, it_2] = newton(x0, nmax, toll, f, df)

x0 = 6;
[xvect_3, it_3] = newton(x0, nmax, toll, f, df)

%%

alpha1_ex = 1;

sol_ex_1 = alpha1_ex * ones(size(xvect_1));
err_1 = abs(xvect_1 - sol_ex_1);

sol_ex_1_m = alpha1_ex * ones(size(xvect_1_m));
err_1_m = abs(xvect_1_m - sol_ex_1_m);

it_1v = 0:it_1;
it_1_mv = 0:it_1_m;

figure
semilogy(it_1v, err_1)
hold on
semilogy(it_1_mv, err_1_m)
grid on
legend("Newton", "Newton modificato")