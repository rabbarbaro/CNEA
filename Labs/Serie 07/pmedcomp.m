function I = pmedcomp(a, b, N, fun)

% I = pmedcomp(a, b, N, f)
% Formula del punto medio composita: 
% IN:
%   - a, b: estremi di integrazione
%   - N: numero di sottointervalli (N=1 formula di integrazione semplice)
%   - f: funzione da integrare definita come anonymous
% OUT:
%   - I: integrale calcolato

% calcolo la lunghezza di ogni sottointervallo e il vettore dei punti medi
h = (b - a) / N;
xm_k = (a + h/2):h:(b - h/2);
% calcolo l'integrale
I = h * sum(fun(xm_k));

% alternativa (non necessaria poiché il calcolo è già ottimizzato sopra):
% I = 0;
% for x = xk
%     I = I + h*fun(x);
% end

end