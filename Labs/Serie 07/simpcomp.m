function I = simpcomp(a, b, N, fun)

% I = simpcomp(a, b, N, f)
% Formula di Simpson composita: 
% IN:
%   - a, b: estremi di integrazione
%   - N: numero di sottointervalli (N=1 formula di integrazione semplice)
%   - f: funzione da integrare definita come anonymous
% OUT:
%   - I: integrale calcolato

% calcolo la lunghezza di ogni sottointervallo, il vettore dei nodi e il
% vettore dei punti medi
h = (b - a) / N;
x_k = a:h:b;
xk_m = (a + h/2):h:(b - h/2);
% calcolo l'integrale
I = h/6 * (sum(fun(x_k(1:end-1))) + 4 * sum(fun(xk_m)) + sum(fun(x_k(2:end))));

end