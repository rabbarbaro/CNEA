clear
clc
close all

k = 15;
a = 1:10+k;

mediaArmonica = harmonicMean(a);
% il formato short ha 4 cifre decimali fisse
fprintf("%.4f\n", mediaArmonica)

function hm = harmonicMean(a)
    % estraggo la lunghezza e applico la definizione data nel testo
    n = length(a);
    hm = n / sum(a.^-1);
end 
