clear
clc
close all

M1 = 10;
err_interp_max1 = 1e-2;
err_interp_max2 = 1e-5;

% max|f(x) - s3(x)| <= C * H^4 * max|f''''(x)|

% err_interp_max1 = C_dot_max_d4f * 1/M1
C_dot_max_d4f = err_interp_max1 * M1;

% err_interp_max1 = C_dot_max_d4f * 1/M2
M2 = C_dot_max_d4f / err_interp_max2