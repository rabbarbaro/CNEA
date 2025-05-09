function [t_h, u_h, iter_pf] = eulero_indietro_pto_fisso(f, t_max, y_0, h)

% [t_h, u_h, iter_pf] = eulero_indietro_pto_fisso(f, t_max, y_0, h)
%
% Risolve il problema di Cauchy (ODE primo ordine)
%       y' = f(t, y)
%       y(0) = y_0
% utilizzando il metodo di Eulero all'indietro (Eulero implicito)
%       u_(n+1) = u_n + h*f(t_(n+1), u_(n+1))
% e risolvendo l'quazione non lineare con il metodo delle iterazioni di
% punto fisso
% IN:
%   - f: anonymous function che descrive il problema di Cauchy che deve
%       ricevere in ingresso due argomenti: f = f(t, y)
%   - t_max: istante finale dell'intervallo temporale di soluzione, dato
%       l'istante iniziale t_0 = 0
%   - y_0: dato iniziale del problema di Cauchy
%   - h: ampiezza del passo di discretizzazione temporale
% OUT:
%   - t_h: vettore degli istanti in cui si calcola la soluzione discreta
%   - u_h: soluzione discreta calcolata nei nodi temporali t_h
%   - iter_pf: vettore contente il numero di iterazioni del metodo di punto
%       fisso ad ogni passo

% vettore degli istanti in cui risolvo la ODE
t_0 = 0;
t_h = t_0:h:t_max;

% ottengo il numero degli istanti e inizializzo il vettore delle soluzioni
% con il dato iniziale
Nh = length(t_h);
u_h = zeros(1, Nh);
u_h(1) = y_0;

% definisco i parametri per le iterazioni di punto fisso e inizializzo il
% vettore con il numero di iterazioni
nmax = 100;
tol = 1e-5;
iter_pf = zeros(1, Nh);

% ciclo iterativo che calcola u_(n+1) = u_n + h*f(t_(n+1), u_(n+1))
for ii = 1:Nh-1
    % definisco i parametri necessari per le iterazioni di punto fisso
    u_pf = u_h(ii);
    t_pf = t_h(ii+1);
    % definisco la funzione di iterazione di punto fisso
    phi = @(y) u_pf + h*f(t_pf, y);
    % calcolo u_(n+1) con le iterazioni di punto fisso (la guess iniziale
    % sarà u_n)
    [succ, it_pf] = ptofis(u_pf, phi, nmax, tol);

    % aggiorno u_(n+1) e salvo il numero delle iterazioni
    u_h(ii+1) = succ(end);
    iter_pf(ii+1) = it_pf;        
end

end

