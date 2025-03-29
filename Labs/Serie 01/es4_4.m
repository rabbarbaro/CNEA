clear
clc
close all

% numero di punti
n = 1e6;

%% MIA: più veloce, impiega una valanga di memoria (ma tipo tantissima)

% ma tipo con n = 1e9 matlab letteralmente esplode

tic

% riempio una matrice con n coppie di numeri casuali fra 0 e 1
pnt = rand(2, n);
% calcolo il quadrato della distanza dal centro dei punti
pnt_dist = pnt(1, :).^2 + pnt (2, :).^2;
% verifico quali coppie sono nel cerchio e le conto
idx = pnt_dist <= 1;
pi_est = 4 * sum(idx) / n

toc

%% SOL: più lenta (~ 3 volte tanto) ma occupa pochissima memoria

tic

% inizializzo il counter dei punti e l'approssimazione di pi a 0
points_counter = 0;
approximated_pi = 0;

for i = 1:n
    % coordinata x: numero random in [0, 1]
    x = rand;

    % coordinata y: numero random in [0, 1]
    y = rand;

    % controllo che il punto estratto sia all'interno del cerchio
    if ( (x^2 + y^2) <= 1 )
        points_counter = points_counter + 1;
    end
end

% calcolo approssimazione di pi greco
approximated_pi = 4*points_counter/n

toc