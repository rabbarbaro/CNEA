clear
clc
close all

n = 25;
isPerfectSquare(n);

% la funzione non ha outargs, come una 'void' in C
function isPerfectSquare(n)
    % se la differenza fra la radice quadrata e l'intero minore più vicino
    % è non nulla allora non è un quadrato perfetto
    if sqrt(n)-floor(sqrt(n))
        disp("n non è un quadrato perfetto")
    else
        disp("n è un quadrato perfetto")
    end
end