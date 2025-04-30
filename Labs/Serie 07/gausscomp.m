function I = gausscomp(a, b, N, fun)

% I = gausscomp(a, b, N, f)
% Formula di Gauss composita a due nodi:
% IN:
%   - a, b: estremi di integrazione
%   - N: numero di sottointervalli (N=1 formula di integrazione semplice)
%   - f: funzione da integrare definita come anonymous
% OUT:
%   - I: integrale calcolato

% calcolo la lunghezza dei sottointervalli, il vettore dei nodi
% equispaziati (estremi sottointervalli) e il vettore dei punti medi
h = (b - a) / N;
x_k = a:h:b;
xm_k = (a + h/2):h:(b - h/2);

% nodi di quadratura per ogni sottointervallo
y0 = xm_k + h/2 * (-1/sqrt(3));
y1 = xm_k + h/2 * (1/sqrt(3));

% pesi
w0 = h/2;
w1 = h/2;

% calcolo l'integrale
I = w0 * sum(fun(y0)) + w1 * sum(fun(y1));

end