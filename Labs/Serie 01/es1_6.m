clear
clc
close all

% (MIA GUESS) ci aspettiamo una parabola, essendo il log(x) riportato a
% una retta (poi da elevare al quadrato)

% (SOL):  Il grafico x-logaritmico di f(x) = (log(x)).^2 è formato dai
% punti del piano (log(x),(log(x))2) (e non dai punti (x,f(log(x)))!!).
% Quindi il grafico che si ottiene è una parabola con la concavità
% verso l’alto (tipo y = x.^2).

f = @(x) (log(x)).^2;

x = 0.1:0.01:10;

semilogx(x, f(x));