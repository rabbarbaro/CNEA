clear
clc
close all

n = 1;
a = -2;
b = 2;

f = @(x) exp(x) + exp(sin(2*pi*x));

M = 8;
h = (b-a)/M;

I = 0;
for x_k = a:h:b-h
    [x, w] = zplege(n, x_k, x_k + h);
    I = I + sum(w .* f(x));
end

I