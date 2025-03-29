clear
clc
close all

iter = 10;
S = 1e5;
% sarebbe x0
x(1) = 1e5;

% la 1 iterazione sta in x(2) e la 10 in x(11)
for n = 1:iter
    x(n+1) = 0.5 * (x(n) + S/x(n));
end