clear
clc
close all

A = hilb(3);

K_A = max(eig(A)) / min(eig(A));
d = (K_A - 1) / (K_A + 1);

% fatt. abb. = d^k
% d^k = 1/200

k = log(1/200) / log(d)