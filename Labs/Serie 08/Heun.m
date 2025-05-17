function [t_h, u_h] = Heun(f, t_max, y_0, h)

% [t_h, u_h] = Heun(f, t_max, y_0, h)
%
% Risolve il problema di Cauchy (ODE primo ordine)
%       y' = f(t, y)
%       y(0) = y_0
% utilizzando il metodo di Heun
%       u*_(n+1) = u_h(ii) + h*f(t_h(ii), u_h(ii))
%       u_(n+1) = u_n + h/2 * (f(t_n, u_n) + f(t_(n+1), u*_(n+1)))
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

% vettore degli istanti in cui risolvo la ODE
t_0 = 0;
t_h = t_0:h:t_max;

% ottengo il numero degli istanti e inizializzo il vettore delle soluzioni
% con il dato iniziale
Nh = length(t_h);
u_h = zeros(1, Nh);
u_h(1) = y_0;

% ciclo iterativo che calcola u*_(n+1), u_(n+1)
for ii = 1:Nh-1
    u_h(ii+1) = u_h(ii) +  h/2 * (f(t_h(ii), u_h(ii)) + f(t_h(ii+1), u_h(ii) + h*f(t_h(ii), u_h(ii))));
end

end