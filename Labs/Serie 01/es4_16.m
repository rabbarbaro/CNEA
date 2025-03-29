clear
clc
close all

% definisco tutti i dati
a = 2;
b = 5;
k = 2;

% chiamo la funzione
res = fun_es4_16(a, b, k);

% verifico che sia corretto
% n = mod(k, 4) + 3;
% isequal(res, (a+b)^n);

function out = fun_es4_16(a, b, k)
    % calcolo n come indicato e prealloco il vettore v
    n = mod(k, 4) + 3;
    v = zeros(1, n);

    % calcolo un vettore con tutti i valori e poi sommo i suoi valori
    for jj = 0:n
        v(jj+1) = (factorial(n)/(factorial(jj) * factorial(n-jj)) ) ...
            * a^jj * b^(n - jj);
    end
    out = sum(v);

    % alternativamente:
    % sol=0;
    % for i=0:n
	%     sol=sol+factorial(n)/factorial(i)/factorial(n-i).*a^i.*b^(n-i);
    % end
end