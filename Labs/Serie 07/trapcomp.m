function I = trapcomp(a, b, N, fun)

% I = trapcomp(a, b, N, f)
% Formula del trapezio composita:
% IN:
%   - a, b: estremi di integrazione
%   - N: numero di sottointervalli (N=1 formula di integrazione semplice)
%   - f: funzione da integrare definita come anonymous
% OUT:
%   - I: integrale calcolato

% calcolo la lunghezza di ogni sottointervallo e il vettore dei nodi
h = (b - a) / N;
x_k = a:h:b;
% calcolo l'integrale
I = h/2 * (fun(x_k(1)) + fun(x_k(end))) + h * sum(fun(x_k(2:end-1)));

end