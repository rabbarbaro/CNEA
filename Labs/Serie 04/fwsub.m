function [y] = fwsub(L, b)

% [y] = fwsub(L, b)
% Sostituzione in avanti per risoluzione di sistemi triangolari inferiori
% IN
%   - L: matrice quadrata triangolare inferiore
%   - b: vettore termini noti
% OUT
%   - y: soluzione del sistema Ly = b

%%  verifica degli input

n = length(b);

% verifica delle dimensioni
if size(L,1) ~= n || size(L,2) ~= n
    error('dimensioni incompatibili')
end

% verifica se L è triangolare inferiore
if ~isequal(L, tril(L))
    error('L non triangolare inferiore')
end

% verifica se la matrice è invertible: se è triangolare (verificato nella
% condizione precedente) allora il determinante è la produttoria degli
% elementi sulla diagonale, se uno degli elementi è nullo il determinante è
% nullo e dunque la matrice non è invertibile
if prod(diag(L)) == 0
    error('matrice non invertibile')
end

%% inizializzazione

% prealloco il vettore soluzione e definisco il primo elemento
y = zeros(n, 1);
y(1) = b(1)/L(1, 1);

%% algoritmo

% algoritmo "base"
% for ii = 2:n
%     % per ogni indice ii, faccio la sommatoria
%     s = 0;
%     for jj = 1:ii-1
%         s = s + L(ii, jj) * x(jj);
%     end
%     % definisco tutti gli elementi del vettore soluzione
%     x(ii) = (f(ii) - s) / L(ii, ii);
% end

% notiamo che nell'algoritmo la sommatoria equivale al prodotto scalare
% fra il vettore riga ii-esima di L e il vettore colonna dei primi ii-1
% elementi del termine noto
for ii = 2:n
    y(ii) = (b(ii) - L(ii, 1:ii-1) * y(1:ii-1)) / L(ii, ii);
end

end