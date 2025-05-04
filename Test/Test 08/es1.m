clear
clc
close all

a = -2;
b = 2;
n = 5;

x_q = linspace(a, b, 1000);

f = @(x) 1 - exp(sin(x));

x_nod = linspace(a, b, n+1);
f_nod = f(x_nod);

P = polyfit(x_nod, f_nod, n);
P_dis = polyval(P, x_q);

err = abs(f(x_q) - P_dis);

[maxErr, idx] = max(err);
maxErr
x_q(idx)

plot(x_q, f(x_q), x_q, P_dis, x_q, err)