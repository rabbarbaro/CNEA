clear
clc
close all

% definisco i parametri per il problema di cauchy (problema modello)
y0 = 2;
lambda = -42;

% definisco l'intervallo di tempo
t0 = 0;
t_max = 1;

% definisco la funzione f del problema di cauchy e la soluzione esatta
f = @(t, y) lambda * y;
y_ex = @(t) y0*exp(lambda*(t - t0));
% definisco la derivata parziale di f rispetto a y per poter applicare il
% metodo di newton all'interno di EI
df_dy = @(t, y) lambda;
% creo il vettore di query su istanti di tempo arbitrari (1000)
t_dis = linspace(0, t_max, 1000);

%% 1

% plotto la soluzione esatta
plot(t_dis, y_ex(t_dis))

%% 2

% con passo 0.05 calcolo la soluzione numerica del problema con EA
h1 = 0.05;
[t_h1_EA, u_h1_EA] = eulero_avanti(f, t_max, y0, h1);

% plotto confrontando con la soluzione esatta
subplot(1, 2, 1)
plot(t_dis, y_ex(t_dis), t_h1_EA, u_h1_EA)
legend('Soluzione esatta', 'Eulero Avanti')
title('Eulero Avanti con h = 0.05')

% notiamo che la soluzione non converge a 0
% sappiamo che EA è condizionatamente assolutamente stabile, in particolare
% solo se:
%       |lambda*h + 1| < 1
% ovvero in caso reale solo se:
%       -1 < lambda*h < 0
% e in questo caso abbiamo che lambda*h = -2.1, quindi non entra nella
% regione di assoluta stabilità

% con passo 0.05 calcolo la soluzione numerica del problema con EI usando
% come risolutore dell'equazione non lineare il metodo di newton
[t_h1_EI, u_h1_EI, vett_it_newton] = eulero_indietro_newton(f, df_dy, t_max, y0, h1);

% plotto confrontando con la soluzione esatta
subplot(1, 2, 2)
plot(t_dis, y_ex(t_dis), t_h1_EI, u_h1_EI)
legend('Soluzione esatta', 'Eulero Indietro')
title('Eulero Indietro (Newton) con h = 0.05')

% notiamo che la soluzione converge a 0
% sappiamo che EI è incondizionatamente assolutamente stabile, in
% particolare
%       |lambda*h - 1| > 1
% ovvero in caso reale solo se:
%       lambda*h < 0
% e in questo caso abbiamo che lambda*h = -2.1, quindi entra nella
% regione di assoluta stabilità

%% 3

% modifico il passo a 0.01 ed eseguo di nuovo per i due metodi
h2 = 0.01;
[t_h2_EA, u_h2_EA] = eulero_avanti(f, t_max, y0, h2);
[t_h2,u_h2,vett_it_newton] = eulero_indietro_newton (f, df_dy, t_max, y0, h2);

% plotto le soluzioni
figure
plot(t_dis, y_ex(t_dis), t_h2_EA, u_h2_EA, t_h2, u_h2)
legend('Soluzione esatta', 'Eulero Avanti', 'Eulero Indietro')
title('Eulero Avanti e Indietro (Newton) con h = 0.01')

% questa volta h * lambda = -0.42, quindi risulta assolutamente stabile
% anche il metodo di EA

% notiamo ancora che EA ha convergenza "da sotto" (sottostima) mentre EI ha
% convergenza "da sopra" (sovrastima)

%% 4.1

% chiamando i metodi di EA ed EI (punto fisso) per h1 = 0.05
h1 = 0.05;
[t_h1_EA, u_h1_EA] = eulero_avanti(f, t_max, y0, h1);
[t_h1_EI, u_h1_EI, iter_pf_h1] = eulero_indietro_pto_fisso(f, t_max, y0, h1);

% plotto le soluzioni trovate
figure
subplot(1, 2, 1)
plot(t_dis, y_ex(t_dis), t_h1_EA, u_h1_EA)
legend('Soluzione esatta', 'Eulero Avanti')
subplot(1, 2, 2)
plot(t_dis, y_ex(t_dis), t_h1_EI, u_h1_EI)
legend('Soluzione esatta', 'Eulero Indietro')
sgtitle('Eulero Avanti e Indietro (Punto Fisso) con h = 0.05')

% come osservato prima, EA non è assolutamente stabile

% osserviamo un errore di overflow nella soluzione di EI con punto fisso,
% bisogna prestare attenzione anche alla regione di convergenza del metodo
% utilizzato per risolvere le equazioni non lineari!
%   - la scelta del metodo utilizzato per la risoluzione dell equazione
%     (non lineare in generale) associata ad ogni passo di un metodo
%     implicito può modificare le proprietò di convergenza generali

% applicando punto fisso, ovvero con la funzione di iterazione
%       phi(y) = u_n + h*f(t_(n+1), y)
% che nel caso del problema modello diventa
%       phi(y) = u_n + h*lamba*y
% ricordiamo che deve valere la condizione
%       |phi'(u_(n+1))| < 1
% affinché il metodo delle iterazioni di punto fisso converga
% che nel caso del problema modello diventa
%       |h * lambda| < 1
% che nel piamo complesso è una circonferenza unitaria centrata in (0,0)

% essendo per h1 = 0.05
%       lambda*h1 = -2.1
% fuori dalla circonferenza descritta, il metodo non converge non per una
% proprietà del metodo in sé, quanto del risolutore non lineare

%% 4.2

% chiamando i metodi di EA ed EI (punto fisso) per h2 = 0.03
h2 = 0.03;
[t_h2_EA, u_h2_EA] = eulero_avanti(f, t_max, y0, h2);
[t_h2_EI, u_h2_EI, iter_pf_h2] = eulero_indietro_pto_fisso(f, t_max, y0, h2);

% plotto le soluzioni trovate
figure
subplot(1, 2, 1)
plot(t_dis, y_ex(t_dis), t_h2_EA, u_h2_EA)
legend('Soluzione esatta', 'Eulero Avanti')
subplot(1, 2, 2)
plot(t_dis, y_ex(t_dis), t_h2_EI, u_h2_EI)
legend('Soluzione esatta', 'Eulero Indietro')
sgtitle('Eulero Avanti e Indietro (Punto Fisso) con h = 0.03')

% per quanto oscillatoria e con sovraelongazione abbastanza elevata, EA
% converge a 0 per t crescenti

% anche in questo caso EI non converge a 0 a causa delle condizioni del
% metodo delle iterazioni di punto fisso (lambda*h2 = -1.26)

%% 4.3

% chiamando i metodi di EA ed EI (punto fisso) per h3 = 0.01
h3 = 0.01;
[t_h3_EA, u_h3_EA] = eulero_avanti(f, t_max, y0, h3);
[t_h3_EI, u_h3_EI, iter_pf_h3] = eulero_indietro_pto_fisso(f, t_max, y0, h3);

% plotto le soluzioni trovate
figure
subplot(1, 2, 1)
plot(t_dis, y_ex(t_dis), t_h3_EA, u_h3_EA)
legend('Soluzione esatta', 'Eulero Avanti')
subplot(1, 2, 2)
plot(t_dis, y_ex(t_dis), t_h3_EI, u_h3_EI)
legend('Soluzione esatta', 'Eulero Indietro')
sgtitle('Eulero Avanti e Indietro (Punto Fisso) con h = 0.01')

% in questo caso il passo è abbastanza fino da evitare instabilità sia per
% EI che per EA

% in tutti questi casi EA sottostima ed EI sovrastima