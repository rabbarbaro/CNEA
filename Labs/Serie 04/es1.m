clear
clc
close all

% definiamo la funzione (e l'intervallo in cui valutarla) come da testo
f = @(x) x.^3 - (2 + exp(1)) * x.^2 + (2*exp(1) + 1) * x + (1 - exp(1)) - cosh(x - 1);
x = 0.5:0.01:6.5;

%% 1

% disegnamo il grafico della funzione ed evidenzio gli zeri della funzione
plot(x, f(x))
grid on
hold on
plot(x, 0*x, '--')

%% 2

% osserviamoche primo zero (~1) non verrà trovato, mentre il secondo (3~4)
% e il terzo (6~6.5) sì

% verifichiamo
if f(0.5) * f(1.5) < 0
    disp("Bisezione applicabile per alpha1")
else
    disp("Bisezione NON applicabile per alpha1")
end
if f(3) * f(4) < 0
    disp("Bisezione applicabile per alpha2")
else
    disp("Bisezione NON applicabile per alpha1")
end
if f(6) * f(6.5) < 0
    disp("Bisezione applicabile per alpha3")
else
    disp("Bisezione NON applicabile per alpha1")
end

%% 3

% creiamo la funzione che implemennta il metodo di bisezione
% vedi file bisez.m

%% 4

% definiamo la tolleranza
toll = 1e-12;

% per il primo zero non è possibile applicarlo, darebbe errore
% a1 = 0.5;
% b1 = 1.5;
% [xvect1, it1] = bisez(a1, b1, toll, f);

% applichiamo per il secondo zero
a2 = 3;
b2 = 4;
[xvect2, it2] = bisez(a2, b2, toll, f);

% applichiamo per il terzo zero
a3 = 6;
b3 = 6.5;
[xvect3, it3] = bisez(a3, b3, toll, f);