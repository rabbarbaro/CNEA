clear
clc
close all

% definisco la funzione e il suo intervallo
f = @(x) atan(7*(x - pi/2)) + sin((x - pi/2).^3);
a = -1;
b = 6;
x = linspace(a, b, 100);

%% 1

% plotto nell'intervallo
plot(x, f(x), x, 0*x)
% visivamente la derivata non si annulla nello zero

%% 2

alpha = pi/2;

% derivo la funzione e ne valuto la derivata nello zero
df = @(x) 7/(49/4 * (pi - 2*x).^2 + 1) + cos((x - pi/2).^3) * 3 * (x - pi/2).^2;
df(alpha)
% conferma analitica che la derivata prima è diversa da 1, zero semplice

% definisco tolleranza e un numero di iterazioni massimo
toll = 1e-3;
nmax = 1000;

% applico newton con guess iniziale x0 = 1.5
x0 = 1.5;
[xvectN1, itN1] = newton(x0, nmax, toll, f, df, 1);
[p1, c1] = stimap(xvectN1);

% applico newton con guess iniziale x0 = 4
x0 = 4;
[xvectN2, itN2] = newton(x0, nmax, toll, f, df, 1);
[p4, c4] = stimap(xvectN2);

err_ass_N1 = alpha - xvectN1(end);
err_ass_N2 = alpha - xvectN2(end);
% con il secondo valore di x0 siamo troppo lontani e newton non converge! 

%% 3

% definisco la nuova tolleranza e calcolo con bisezione sull'intervallo
toll = (b - a)/2^31;
[xvectB, itB] = bisez(a, b, toll, f);

err_ass_B = alpha - xvectB(end);

%% 4

% implementiamo la funzione che esegue qualche iterazione di bisezione per
% avvicinarsi a sufficienza allo zero per poi usare il risultato ottenuto
% come stima iniziale per newton
% vedi file biseznewton.m

%% 5

% calcolo con il metodo bisezione + newton
iterBisez = 5;
[xvectBN, itBN] = biseznewton(a, b, iterBisez, nmax, toll, f, df);