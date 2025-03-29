clear
clc
close all

K = 8;
v_k = zeros(1, K+1);

% con un ciclo for riempio gli elementi
for k = 0:K
    v_k(k + 1) = (2 * k + 1)^2;
end