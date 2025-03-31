clear
clc
close all

F = @(x1, x2) [6*x1 - 2*exp(-2*x1.*x2).*x2 - 1/4; ...
               2*x2 - 2*exp(-2*x1.*x2).*x1 - 1/6];

JF = @(x1, x2) [6 + 4*exp(-2*x1.*x2).*x2^2, exp(-2*x1.*x2).*(4*x1.*x2-2); ...
               exp(-2*x1.*x2).*(4*x1.*x2-2), 2 + 4*exp(-2*x1.*x2).*x1^2,];

alpha = [0.099299729019640, 0.179161952163217]';

x0 = [-0.14, 0.14]';

kmax = 5;

x = x0;
xvec = [];

for k = 1:kmax
    Fvec = F(x(1), x(2));
    Jmat = JF(x(1), x(2));
    x = x - Jmat\Fvec;
    xvec = [xvec x];
end

det(JF(alpha(1), alpha(2)))         % deve essere != 0

xvec(:, 5)

% CONTROLLARE
% normxvec = ((xvec(1, :) + xvec(2, :))/2)';
% [p, c] = stimap(normxvec);

errvec = [];
for k = 1:kmax-1
    errvec = [errvec, norm(xvec(:, k) - alpha)];
end

conv_ord1 = errvec(2:end) ./ errvec(1:end-1)
% la convergenza di ordine 1 tende a 0 --> ordine superiore
conv_ord2 = errvec(2:end) ./ errvec(1:end-1).^2
% non tende a zero: ordine 2
% lim ||x^(k+1)-alpha||/||x^(k)-alpha||^p = mu
