clear
clc
close all

a = -1;
b = 1;
f = @(x) cos(pi * x);
N = 1;

I = simpcomp(a, b, N, f)