clear
clc
close all

iter = 10;
a = zeros(1, iter);
% inizializzo il primo elemento della successione
a(1) = 1;

% genero la successione come dal testo
for n = 2:iter
    a(n) = ((a(n-1))^2 + 2) / (2 * a(n-1));
end

% plotto
plot(a, 'x-')
hold on
plot(1:iter, sqrt(2).*ones(iter))