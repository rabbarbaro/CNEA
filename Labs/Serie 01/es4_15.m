clear
clc
close all

% per numeri pari converge a 0.6, per numeri dispari oscilla fra due
% valori, di fatto conta solo se è pari o dispari
k = 6;
% calcolo r
r = 2.5 + mod(k,2);

% definisco lunghezza, inizializzo il primo valore
len = 100;
a = zeros(1, len);
a(1) = 0.5;

idx = 0;

% calcolo la successione
for n = 2:len
    a(n) = r * a(n-1) * (1 - a(n-1));
    % se trovo l'elemento ne salvo l'indice e blocco l'acquisizione
    if abs(a(n) - a(n-1)) < 1e-3 && ~idx
        idx = n;
    end
end

% plotto e printo
plot(1:len, a);
if idx
    fprintf ('Elemento %f trovato in posizione %d\n', a(idx), idx)
else
    disp('Elemento non trovato')
end