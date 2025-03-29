clear
clc
close all

n = 10;
T = n2triu(n);

function T = n2triu(n)
    T = zeros(n);
    % sommo ogni volta la (sopra)diagonale in questione
    for ii = 1:n
        T = T + diag(ii.*ones(n+1-ii, 1), ii-1);
    end
    % la soluzione invece costruisce elemento per elemento (MIGLIORE)
    % for ii = 1:n
    %     T(ii, ii:n) = 1:n-ii+1;
    % end
end