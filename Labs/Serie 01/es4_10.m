clear
clc
close all

% effettivamente si caga addosso per iterazioni >= 25 a causa di
% propagazione degli errori in aritmetica floating point 
iter = 20;

% inizializzo a, b, p, q
a_n(1) = sqrt(2);
b_n(1) = 2;
p_n(1) = 2 * a_n(1);
q_n(1) = 2 * b_n(1);

% itero secondo l'algoritmo
for n = 2:iter
    a_n(n) = sqrt(2) * sqrt(1 - sqrt(1 - 0.25 * a_n(n - 1) ^2));
    b_n(n) = a_n(n) / sqrt(1 - 0.25 * a_n(n - 1) ^2);
    p_n(n) = 2^n * a_n(n);
    q_n(n) = 2^n * b_n(n);
end

% plotto
plot(1:iter, p_n, 1:iter, q_n)
hold on
plot(1:iter, pi.*ones(iter), color="#EDB120")
legend('p_n', 'q_n', 'pi')