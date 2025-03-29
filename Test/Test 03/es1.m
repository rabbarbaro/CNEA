clear
clc
close all

n = 100;
A = diag(2 * ones(n, 1)) - diag(ones(n-1, 1), 1);
A(n, 1) = 1;
b = ones(n, 1);

A = sparse(A);

[L, U, P, Q] = lu(A);

isequal(L * U, P*A*Q)

y_star = L\(P*b);
x_star = U\y_star;

%% es 5 non ho sbatta di fare un nuovo file

clc
clear

n = 3;
A = [6 -2 -2
    -2 8 -4
    -2 -4 10];
b = ones(n, 1);
x0 = b;
toll = 1e-12;
nmax = 2;
alpha = 1;

P = [6 0 0
    -1 8 0
    -1 -2 10];

[x, k] = richardson(A, b, P, x0, toll, nmax, 1)

% A \ b