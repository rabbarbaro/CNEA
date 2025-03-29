clear
clc
close all

k = 15;
a = 1:10+k;

mediaGeometrica = geometricMean(a);
% il formato short ha 4 cifre decimali fisse
fprintf("%.4f\n", mediaGeometrica)

function gm = geometricMean(a)
    % estraggo la lunghezza e applico la definizione data nel testo
    n = length(a);
    gm = prod(a)^(1/n);
end