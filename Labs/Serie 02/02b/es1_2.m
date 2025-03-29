clear
clc
close all

%% 1

n = 7;

% costruiamo A e b
A = diag(9 * ones(n, 1)) - diag(3 * ones(n-1, 1), 1) ...
    - diag(3 * ones(n-1, 1), -1) + diag(ones(n-2, 1), 2) ...
    + diag(ones(n-2, 1), -2);
b = [7 4 5 5 5 4 7]';

%% 2

% verifica dominanza diagonale stretta per righe

% avevamo già visto questo metodo per la verifica
% % inizializzo un vettore colonna vuoto
% r = zeros(n, 1);
% for ii = 1:n
%     % sommo gli elementi extradiagonali (con indice di colonna diverso da
%     % jj) sulle varie righe e salvo in un vettore
%     r(ii) = sum(abs(A(ii, [1:ii-1, ii+1:n])));
% end
% % if su un vettore è vero solo se sono verificate tutte le condizioni
% if diag(A) - r > 0
%     disp('A a dominanza diagonale per righe');
% else
%     disp('A non a dominanza diagonale per colonne');
% end

% questo metodo è più veloce (operazioni linearizzate)
A_diag = diag(abs(A));
A_outdiag = sum(abs(A), 2) - A_diag;
if(A_diag > A_outdiag)
    disp("A DDSR")
end

%% 3

% verifica SPD

% avevamo già visto questo metodo per la verifica (crit. Sylvester)
% % verifico se è simmetrica
% if isequal(A, A')
%     disp('A simmetrica');
%     nA = size(A, 1);
% 
%     % calcolo il determinante delle sottomatrici (criterio di Sylvester)
%     for ii = 1:nA
%         if det(A(1:ii, 1:ii)) > 0
%             % se siamo arrivati all'ultima iterazione allora le condizioni
%             % su tutti i determinanti sono verificate
%             if ii == nA
%                 disp('A DP');
%             end
%         else
%             error('A non DP');
%         end
%     end
% else
%     disp('A non simmetrica');
% end

% questo metodo è più veloce (calcolo autovalori in matlab ottimizzato)
if A == A'
    if eig(A) > 0
        disp("A SDP")
    end
end

%% 4

% creiamo la funzione che applica il metodo di Gauss-Seidel
% vedi file gs.m

%% 5

% definiamo i dati
x0 = zeros(n, 1);
toll = 1e-6;
nmax = 1e3;

% risolviamo con g-s
[x_gs, k_gs] = gs(A, b, x0, toll, nmax);

%% 6

% risolviamo con j
[x_j, k_j] = jacobi(A, b, x0, toll, nmax);

fprintf("Iterazioni G-S: %d\n", k_gs)
fprintf("Iterazioni J: %d\n", k_j)

% k_j > k_gs, il raggio spettrale della matrice di iterazione g-s sarà
% minore del raggio spettrale della matrice di iterazione j

% calcoliamo i raggi spettrali

% la matrice di iterazione g-s è B = I - T^-1 * A
P_gs = tril(A);
B_gs = eye(n) - P_gs \ A;       % inv(P_gs) * A
rho_B_gs = max(abs(eig(B_gs)));

% la matrice di iterazione j è B = I - D^-1 * A
P_j_inv = diag(1 ./ diag(A));
B_j = eye(n) - P_j_inv * A;
rho_B_j = max(abs(eig(B_j)));

fprintf("Raggio spettrale B G-S: %f\n", rho_B_gs)
fprintf("Raggio spettrale B J: %f\n", rho_B_j)

% entrambi sono minori di uno, convergono entrambi

% se invece di una tolleranza sul residuo normalizzato volessimo imporre
% una tolleranza sull'errore
% e_k = x_k - x_ex
% || e_k || <= |rho(B)|^k * || e_0 ||
% k_min = log(toll) / log(rho(B)) t.c. || e_k || / || e_0 || < toll

k_min_j = log(toll) / log(rho_B_j);
k_min_gs = log(toll) / log(rho_B_gs);

fprintf("Iterazioni minime per tolleranza su errore G-S: %f\n", k_min_j)
fprintf("Iterazioni minime per tolleranza su errore B J: %f\n", k_min_gs)