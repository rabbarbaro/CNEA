clear
clc
close all

a = 3;
b = 5;
c = 4;

isRightTriangle(a, b, c);

function isRightTriangle(a, b, c)
    % metto i lati in un vettore e li ordino (l'ipotenusa è l'ultimo)
    sides = [a, b, c];
    sides = sort(sides);
    %verifico il teorema di pitagora
    if sides(1)^2 + sides(2)^2 == sides(3)^2
        disp("Il triangolo è rettangolo")
    else
        disp("Il triangolo non è rettangolo")
    end
end