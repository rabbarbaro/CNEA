function [L, U, x] = thomas(A, b, isTriDiag)

% [L, U, x] = thomas(A, b, (opt) isTriDiag)
% Calcola la fattorizzazione di Thomas di A e risolve il sistema Ax = b
% IN
%   - A: matrice tridiagonale quadrata
%   - b: vettore dei termini noti
%   - (opt) isTriDiag: opzionale, verifica se la matrice è tridiagonale
% OUT
%   - L: matrice triangolare inferiore della fattorizzazione LU di A
%   - U: matrice triangolare superiore della fattorizzazione LU di A
%   - x: soluzione del sistema Ax = b 

%%  verifica degli input

n = length(b);

% verifica delle dimensioni
if size(A,1) ~= n || size(A,2) ~= n
    error('dimensioni incompatibili')
end

% opzionale perché computazionalmente molto intensa
if nargin == 3 && isTriDiag
    % verifica tridiagonalità
    for ii = 2:n
        if any(diag(A, ii)) || any(diag(A, -ii))
            error('A non tridiagonale')
        end
    end
end

%% inizializzazione

% dimensione del sistema
n = size(A, 1);

% sopradiagonale, sottodiagonale e diagonale principale di A
c = diag(A, 1);
e = diag(A, -1);
a = diag(A);

% prealloco L ed U inserendo le (sopra)diagonali costanti corrispondenti
L = diag(ones(n, 1));
U = diag(c, 1);

% prealloco i vettori dei parametri alpha e delta della fattorizzazione
alpha = zeros(n, 1);
delta = zeros(n-1, 1);

% prealloco i vettori soluzione
y = zeros(n, 1);
x = zeros(n, 1);

%% algoritmo fattorizzazione

% definisco il primo elemento di alpha
alpha(1) = A(1,1);

% calcolo alpha e delta per ogni iterazione
for ii = 2:n
    delta(ii-1) = e(ii-1) / alpha(ii-1);
    alpha(ii) = a(ii) - delta(ii-1)*c(ii-1);
end

% completo le matrici L ed U
L = L + diag(delta, -1);
U = U + diag(alpha);

%% algoritmo soluzione

% definisco il primo elemento di y
y(1) = b(1);

% calcolo y per ogni iterazione
for ii = 2:n
    y(ii) = b(ii) - delta(ii-1)*y(ii-1);
end

% definisco l'ultimo elemento di x
x(n) = y(n) / alpha(n);

% calcolo x per ogni iterazione
for ii = n-1:-1:1
    x(ii) = (y(ii) - c(ii)*x(ii+1))/alpha(ii);
end