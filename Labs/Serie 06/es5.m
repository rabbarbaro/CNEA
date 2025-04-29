clear
clc
close all

% definisco la funzione e l'intervallo (vettore query)

f = @(x) exp(-x.^2) .*sin(x);

a = -2;
b = 3;
x_q = linspace(a, b, 1000);
f_q = f(x_q);

%% 1

% definisco numero di sottointervalli e loro ampiezza
n = 3;
h = (b - a)/n;
% definnisco nodi equispaziati e valuto la funzione nei nodi
x_nod = a:h:b;
f_nod = f(x_nod);
% faccio l'interpolazione lineare sui 3 sottointervalli per ogni puto del
% vettore query
v_q = interp1(x_nod, f_nod, x_q);

% plotto la funzione
plot(x_q, f_q)
hold on
% plotto l'interpolazione e i nodi
plot(x_q, v_q)
plot(x_nod, f_nod, '*')
legend('f(x)', 'Interpolante', 'Nodi')
hold off

%% 2

% per calcolare la norma infinito dell'errore cerco il suo massimo
err = max(abs(f_q - v_q));
fprintf('Errore in norma infinito (massimo): %f\n', err)

%% 3

% inizializzo vettori contenenti la lunghezza dei sottointervalli e degli
% errori
H = [];
err = [];

% itero sul numero di sottointervalli 4,8,16,32,64,128
for n = 2 .^ (2:7)
    % calcolo la lunghezza dei sottointervalli e salvo nel vettore
    h = (b - a)/n;
    H = [H, h];
    % definnisco nodi equispaziati e valuto la funzione nei nodi
    x_nod = a:h:b;
    f_nod = f(x_nod);
    % faccio l'interpolazione lineare sui 3 sottointervalli per ogni puto
    % del vettore query
    v_q = interp1(x_nod, f_nod, x_q);
    % calcolo il massimo dell'errore e lo salvo nel vettore
    err = [err, max(abs(f_q - v_q))];
end

% plotto in scala logaritmica errore VS intervallo calcolato
loglog(H, err, 'o-', H, H, '--', H, H.^2, '--')
legend('Errore', 'H', 'H^2')

% ritroviamo il risultao teorico, ordine di convergenza 2