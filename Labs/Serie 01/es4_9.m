clear
clc
close all

n = 3;
sqrtn = sqrtErone(n);

function sqrtn = sqrtErone(n)
    % definisco tolleranza e la differenza alla prima iterazione (di fatto
    % è necessario solo che sia > toll
    toll = 1e-6;
    delta = n;
    % inizializzo la successione
    r_k(1) = n;
    % l'indice 0 non esiste per Matlab, devo partire da 1, il secondo è 2
    k = 2;

    % finché rimango fuori dalla tolleranza itero
    while delta >= toll
        r_k(end+1) = 0.5 * (r_k(k-1) + n / r_k(k-1));
        % uso la differenza dall'iterazione precedente per poi verificare
        % se ho raggiunto la tolleranza richiesta
        delta = abs(r_k(k) - r_k(k-1));
        k = k + 1;
    end

    plot(r_k);
    sqrtn = r_k(end);
end