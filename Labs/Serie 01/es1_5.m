clear
clc
close all

% creo le funzioni con le function handle
f = @(x) 2 + (x - 3) .* sin(5 .* (x - 3));
r1 = @(x) -x + 5;
r2 = @(x) x - 1;

x = 0:.01:6;

% plotto, posso usare hold on o fare tutto in una riga
plot(x, f(x));
hold on
plot(x, r1(x), '--', x, r2(x), '--')