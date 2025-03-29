clear
clc
close all

% uso questo script come verifia

gamma = 1;

A = [2*gamma 2 -8
    gamma 1 8
    2 0 1];

[L, U, P] = lu(A);