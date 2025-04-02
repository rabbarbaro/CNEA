function [xvect, it] = biseznewton(a, b, nmax_b, nmax_n, toll, fun, dfun)

% [xvect, it] = biseznewton(a, b, nmax_b, nmax_n, toll, fun, dfun)
% Cerca lozero di f(x) nell'intervallo [a,b],utilizzando il metodo di
% bisezione per l'avvicinamento allo zero e poi il metodo di Newton
% IN
%   - a, b: estremi intervallo di partenza
%   - nmax_b: numero di iterazioni per il metodo di Bisezione
%   - nmax_n: numero massimo di iterazioni per il metodo di Newton
%   - toll_n: tolleranza sul test d'arresto il metodo di Newton
%   - fun: funzione inline/anonymus
%   - dfun: derivata di fun come funzione inline/anonymus
% OUT
%   - xvect: vettore (colonna) con tutte le iterate calcolate (l'ultima è
%       la soluzione)
%   - it: iterazioni totali (bisezione + newton) effettuate

xvect = [];
it = [];

disp('----------------- mi avvicino allo zero con bisezione ---------------')

% Metodo di Bisezione
toll_b = ( b - a ) / ( 2^( nmax_b + 1) );
[ xvect_b,  it_b ] = bisez( a, b, toll_b, fun );
it = it_b;
xvect = [ xvect_b ];

disp('----------------- lancio il metodo di Newton ---------------')


% Metodo di Newton
xv = xvect( end );
[ xvect_n, it_n ] = newton( xv, nmax_n, toll, fun, dfun );  
it = it + it_n;
xvect = [ xvect; xvect_n ]; 
