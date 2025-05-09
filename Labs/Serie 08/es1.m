clear
clc
close all

% definisco la funzione come da testo e ne calcolo derivate prima e
% seconda in modo analtico
f = @(x) exp(-x.^2) .* sin(2*x + 1);
df = @(x) -2*x.*exp(-x.^2) .* sin(2*x + 1) + 2*exp(-x.^2) .* cos(2*x + 1);
d2f = @(x) 4*x.^2.*exp(-x.^2) .* sin(2*x + 1) + ...
           -8*x.*exp(-x.^2) .* cos(2*x + 1) + ...
           -6*exp(-x.^2) .* sin(2*x + 1);

% definisco il punto in cui valutare la derivata
xb = 0;
% definisco i vari passi
hv = 0.4 ./ (2.^(0:5));

%% 1

% inizializzo i vettori degli errori
err_dfa = [];
err_dfi = [];
err_dfc = [];
% valuto la derivata esatta
df_ex = df(xb);

% ciclo su ogni passo calcolo l'approssimazione della derivata e salvo il
% valore dell'errore nel vettore corrispondente
for h = hv
    % differenze finite in avanti
    dfa = (f(xb + h) - f(xb)) / h;
    err_dfa = [err_dfa, abs(dfa - df_ex)];
    % differenze finite all'indietro
    dfi = (f(xb) - f(xb - h)) / h;
    err_dfi = [err_dfi, abs(dfi - df_ex)];
    % differenze finite centrate
    dfc = (f(xb + h) - f(xb - h)) / (2*h);
    err_dfc = [err_dfc, abs(dfc - df_ex)];
end

% plotto gli errori in scala logaritmica
figure
loglog(hv, err_dfa, hv, err_dfi, hv, err_dfc)
hold on
loglog(hv, hv, hv, hv.^2, LineStyle="--")
legend('Diff Finite Avanti', 'Diff Finite Indietro', 'Diff Finite Centrate', 'h', 'h^2', 'Location', 'NorthWest')

% troviamo il risultato teorico, DFA e DFI hanno ordine 1 e DFC ha ordine 2

%% 2

% inizializzo il vettore degli errori
err_d2fc = [];
% valuto la derivata seconda esatta
d2f_ex = d2f(xb);

% per ogni passo calcolo l'approssimazione della derivata seconda
for h = hv
    % differenze finite centrate (derivata seconda)
    df2c = (f(xb + h) - 2*f(xb) + f(xb - h)) / h^2;
    err_d2fc = [err_d2fc, abs(df2c - d2f_ex)];
end

% plotto in scala logaritmica
figure
loglog(hv, err_d2fc)
hold on
loglog(hv, hv, hv, hv.^2, LineStyle="--")
legend('Diff Finite Centrate (derivata seconda)', 'h', 'h^2', 'Location', 'NorthWest')

% troviamo il risultato teorico, DFC per derivata seconda ha oordine 2