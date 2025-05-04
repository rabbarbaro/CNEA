clear
clc
close all

% Define the function
f = @(x) (x <= pi) .* x/pi + (x > pi) .* (2 - x/pi).^2;

% Query points
xq = linspace(0, 2*pi, 1000);

% Plot the original function
figure;
plot(xq, f(xq), 'LineWidth', 1.5);
hold on;

% Number of points for interpolation
n = 4;
M = floor(n/2);

% Grid points
h = 2 * pi / (n + 1);
x_j = (0:n) * h;

% Preallocate Fourier coefficients
c_k = zeros(1, 2*M + 1);

% Compute Fourier coefficients
for k = -M:M
    S = 0;
    for j = 0:n
        S = S + f(x_j(j+1)) * exp(-1i * k * j * h);
    end
    c_k(k + M + 1) = 1/(n+1) * S;
end

% Trigonometric interpolation function
I = @(x) real(sum(c_k * exp(1i * (-M:M)' * x), 1));

% Plot the interpolated function
plot(xq, I(xq), '--', 'LineWidth', 1.5);
legend('Original Function', 'Interpolated Function');
title('Trigonometric Interpolation');
xlabel('x');
ylabel('f(x)');
grid on;