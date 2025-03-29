clear
clc
close all

% genero le funzioni con il function handle
f = @(x) x .* sin(x) + 0.5 .^ sqrt(x);
g = @(x) x.^4 + log(x.^3 + 1);

x = 0:3;

f(x);
g(x);