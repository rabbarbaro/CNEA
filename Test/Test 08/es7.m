clear
clc
close all

M1 = 10;
err_interp_max1 = 1e-2;
err_interp_max2 = 1e-5;

% H1 = (b-a) / M1

% max|f(x) - s3(x)| <= C * H^4 * max|f''''(x)|

% err_interp_max1 = C_dot_max_d4f_dot_b-a4 * 1/M1^4
Const = err_interp_max1 * M1^4;

% err_interp_max1 = C_dot_max_d4f_dot_b-a4 * 1/M2^4
M2 = (Const / err_interp_max2)^(1/4)
M2 = ceil(M2)