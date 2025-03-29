clear
clc
close all

% inizializzo la sequenza a 1, poi faccio il calcolo (mi fermo a < 88 e non
% <= perché aggiungendo n+1 potrei sforare
n = 1;
while sum(1:n) < 88
    n = n + 1;
end