clear
clc
close all

n = 10;

T = checkBoard(n);

function T = checkBoard(n)
    % inizializzo la matrice a tutti zeri e poi sostituisco a schacchi
    T = zeros(n);
    T(1:2:n, 1:2:n) = 1;
    T(2:2:n, 2:2:n) = 1;
end