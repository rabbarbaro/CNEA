clear
clc
close all

% definisco la condizioone iniziale del problema di cauchy
y0 = 2;

% definisco l'intervallo di tempo
t0 = 0;
t_max = 10;

% definisco la funzione f del problema di cauchy e la soluzione esatta
f = @(t, y) (pi * (cos(pi*t))./(2 + sin(pi*t)) - 0.5) .* y;
y_ex = @(t) (2 + sin(pi*t)) .* exp(-t/2);

%% 1

% scrivo la funzione che implementa il metodo di Eulero in avanti per la
% risoluzione di un generico problema di Cauchy
% vedi eulero_avanti.m

%% 2

% risolvo il problema di cauchy con i due passi 0.05 e 0.01 con EA

h1 = 0.05;
[t_h1_EA, u_h1_EA] = eulero_avanti(f, t_max, y0, h1);

h2 = 0.01;
[t_h2_EA, u_h2_EA] = eulero_avanti(f, t_max, y0, h2);

% plotto la soluzione trovata
plot(t_h1_EA, u_h1_EA, t_h2_EA, u_h2_EA, t_h2_EA, y_ex(t_h2_EA), '--')
legend('Eulero avanti passo h1 = 0.05', 'Eulero avanti passo h2 = 0.01', 'Soluzione esatta', Location='northwest')

% per eulero in avanti abbiamo convergenza da sotto

%% 3

% scrivo la funzione che implementa il metodo di Eulero all'indietro per la
% risoluzione di un generico problema di Cauchy risolvendo l'quazione non
% lineare con il metodo delle iterazioni di punto fisso
% vedi eulero_indietro_pto_fisso.m

%% 4

% risolvo il problema di cauchy con i due passi 0.05 e 0.01 con EI (pto
% fisso)

[t_h1_EI, u_h1_EI, iter_pf_h1] = eulero_indietro_pto_fisso(f, t_max, y0, h1);

[t_h2_EI, u_h2_EI, iter_pf_h2] = eulero_indietro_pto_fisso(f, t_max, y0, h2);

% plotto la soluzione trovata
figure
plot(t_h1_EI, u_h1_EI, t_h2_EI, u_h2_EI, t_h2_EI, y_ex(t_h2_EA), '--')
legend('Eulero indietro passo h1 = 0.05', 'Eulero indietro passo h2 = 0.01', 'Soluzione esatta', Location='northwest')

% per eulero all'indietro abbiamo convergenza da sopra

%% 5

% plotto il numero di iterazioni eseguite dal metodo delle iterazioni di
% punto fisso per risolvere l'equazione non lineare in EI con h1 = 0.05
figure
plot(t_h1_EI, iter_pf_h1, '*-')

%% 6

% inizializzo il vettore coi passi e i vettori con gli errori
h_v = 0.04 ./ (2.^(0:4));
e_h_EA = [];
e_h_EI = [];

% per ogni passo
for h = h_v
    % risolvo il problema di cauchy con EA e EI (pto fisso)
    [t_hv_EA, u_hv_EA] = eulero_avanti(f, t_max, y0, h);
    [t_hv_EI, u_hv_EI, iter_pf_hv] = eulero_indietro_pto_fisso(f, t_max, y0, h);

    % calcolo il massimo degli errori e salvo nei vettori
    e_h_EA = [e_h_EA, max(abs(y_ex(t_hv_EA) - u_hv_EA))];
    e_h_EI = [e_h_EI, max(abs(y_ex(t_hv_EI) - u_hv_EI))];
end

% plotto in scala logaritmica gli errori
figure
loglog(h_v, e_h_EA, h_v, e_h_EI)
hold on
grid on
loglog(h_v, h_v, h_v, h_v.^2, LineStyle="--")
legend('Errore EA', 'Errore EI', 'h', 'h^2')

% troviamo il risultato teorico, hanno entrambi ordine 1