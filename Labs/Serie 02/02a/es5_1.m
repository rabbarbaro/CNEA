clear
clc
close all

% definiamo una matrice sparsa come da testo
n = 20;
A = diag(4*ones(1, n)) - diag(ones(1, n-1), -1) - diag(ones(1, n-1), 1);
A(1, :) = ones(1, n);
A(:, 1) = ones(n, 1);

% facendo la fattorizzazione LU standard (pivoting righe)
[L, U, P] = lu(A);

% definendo A come matrice sparsa posso fare pivoting totale
As = sparse(A);
[Ls, Us, Ps, Qs] = lu(As);

% stampo i pattern delle matrici con spy
tiledlayout(2, 3);
nexttile
spy(A)
title("A")
nexttile
spy(L)
title("L")
nexttile
spy(U)
title("U")
nexttile
% possiamo vedere gli effetti del pivoting totale
spy(Ps*As*Qs)
title("Ps*As*Qs")
nexttile
spy(Ls)
title("Ls")
nexttile
spy(Us)
title("Us")