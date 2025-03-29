function [x] = bksub(U, y)

% [x] = bksub(U, y)
% Sostituzione all'indietro per risoluzione di sistemi triangolari
% superiori
% IN
%   - U: matrice quadrata triangolare superiore
%   - y: vettore termini noti
% OUT
%   - x: soluzione del sistema Ux = y

%%  verifica degli input

n = length(y);

% verifica delle dimensioni
if size(U,1) ~= n || size(U,2) ~= n
    error('dimensioni incompatibili')
end
% verifica se L è triangolare superiore
if ~isequal(U, triu(U))
    error('L non triangolare superiore')
end
% verifica se la matrice è invertible: se è triangolare (verificato nella
% condizione precedente) allora il determinante è la produttoria degli
% elementi sulla diagonale
if prod(diag(U)) == 0
    error('matrice non invertibile')
end

%% inizializzazione

% prealloco il vettore soluzione e definisco l'ultimo elemento
x = zeros(n, 1);
x(n) = y(n)/U(n, n);

%% algoritmo

% notiamo che nell'algoritmo la sommatoria equivale al prodotto scalare
% fra il vettore riga ii-esima di U e il vettore colonna degli elementi del
% termine noto da ii+1 a n (gli ultimi n-ii+1 elementi)
for ii = n-1:-1:1
    x(ii) = (y(ii) - U(ii, ii+1:n) * x(ii+1:n)) / U(ii, ii);
end

end