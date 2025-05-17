function [t_h, u_h, vett_it_newton] = eulero_indietro_newton(f, df_dy, t_max, y_0, h)
    
% [t_h, u_h, vett_it_newton] = eulero_indietro_newton(f, df_du, t_max, valore_iniziale, delta_t)
%
% Risolve il problema di Cauchy (ODE primo ordine)
%       y' = f(t, y)
%       y(0) = y_0
% utilizzando il metodo di Eulero all'indietro (Eulero implicito)
%       u_(n+1) = u_n + h*f(t_(n+1), u_(n+1))
% e risolvendo l'equazione non lineare con il metodo di Newton
% IN:
%   - f: anonymous function che descrive il problema di Cauchy che deve
%       ricevere in ingresso due argomenti: f = f(t, y)
%   - df_dy: anonymous function della derivata parziale di f rispetto a y,
%       è funzione in t e y, df_dy = df_dy(t, y)
%   - t_max: istante finale dell'intervallo temporale di soluzione, dato
%       l'istante iniziale t_0 = 0
%   - y_0: dato iniziale del problema di Cauchy
%   - h: ampiezza del passo di discretizzazione temporale
% OUT:
%   - t_h: vettore degli istanti in cui si calcola la soluzione discreta
%   - u_h: soluzione discreta calcolata nei nodi temporali t_h
%   - vett_it_newton: vettore contente il numero di iterazioni del metodo
%       di newton ad ogni passo

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
vett_it_newton = zeros(1, Nh);

% ciclo iterativo che calcola u_(n+1) = u_n + h*f(t_(n+1), u_(n+1))
for ii = 1:Nh-1
    % definisco i parametri necessari per newton
    u_nwt = u_h(ii);
    t_nwt = t_h(ii+1);
    % definisco l'equazione da risolvere analiticamente e la sua derivata
    F = @(y) y - u_nwt - h * f(t_nwt, y);
    dF = @(y) 1 - h * df_dy(t_nwt, y);
    % calcolo u_(n+1) con newton (la guess iniziale sarà u_n)
    [u_newton, it_newton] = newton (u_nwt, nmax, tol, F, dF);
    
    % aggiorno u_(n+1) e salvo il numero delle iterazioni
    u_h(ii+1) = u_newton(end);
    vett_it_newton(ii+1) = it_newton;
end

