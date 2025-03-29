function det = detTriadiag(A)

% det = detTriadiag(A)
% Calcola il determinante di una matrice tridiagonale A (senza verifica)
% usando la fattorizzazione di Thomas di A
% IN
%   - A: matrice tridiagonale quadrata
% OUT
%   - det: determinante di A

%% inizializzazione

% dimensione del sistema
n = size(A, 1);

% sopradiagonale, sottodiagonale e diagonale principale di A
c = diag(A, 1);
e = diag(A, -1);
a = diag(A);

% prealloco i vettori dei parametri alpha e delta della fattorizzazione
alpha = zeros(n, 1);
delta = zeros(n-1, 1);

%% algoritmo fattorizzazione

% definisco il primo elemento di alpha
alpha(1) = A(1,1);

% calcolo alpha e delta per ogni iterazione
for ii = 2:n
    delta(ii-1) = e(ii-1) / alpha(ii-1);
    alpha(ii) = a(ii) - delta(ii-1)*c(ii-1);
end

%% calcolo determinante

% produttoria dei parametri alpha (diagonale di U)
det = prod(alpha);

end

