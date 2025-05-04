clear
clc
close all

f = @(x) 1 + x;

n = 2;
a = -1;
b = 1;

% I_cc_n/2(f) = a0 + sum_k=1^n/2 ( a_2k / (1 - (2k)^2) )
% a_k = 2/pi * int_0^pi f(cos(theta)) * cos(k * theta) * dtheta
% k = 0, 1, 2, ..., n;

% a0 = 2/pi * int_0^pi f(cos(theta)) * cos(0) * dtheta
% a0 = 2/pi * int_0^pi (1 + cos(theta)) * 1 * dtheta
% a0 = 2/pi * int_0^pi (1 + cos(theta)) * dtheta
% a0 = 2/pi * [theta + sin(theta)]_0^pi
% a0 = 2/pi * [pi + sin(pi) - 0 - sin(0)]
% a0 = 2/pi * pi
a0 = 2;

% a2 = 2/pi * int_0^pi f(cos(theta)) * cos(2 * theta) * dtheta
% a2 = 2/pi * int_0^pi (2 + cos(theta)) * cos(2 * theta) * dtheta
%   int 2cos(2 * theta) da 0 a pi si annulla
%   int cos(theta) * cos(2 * theta) da 0 a pi si annulla
%   dato che cos(theta) è dispari e cos(2 * theta) è dispari se centrate
%   in pi/2
a2 = 0;

% I_cc_n/2(f) = a0 + sum_k=1^n/2 ( a_2k / (1 - (2k)^2) )

I_cc = a0 + a2 / (1 - (2)^2)