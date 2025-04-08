function [x_vec, it] = es5_metodo_iterativo(A, b, x0, mu, tol)

x_vec = x0;
it = 0;
n = size(A, 1);
x = x0;

while norm(b - A * x0)/norm(b) > tol
    for ii = 1:n
        S1 = 0;
        for jj = 1:ii-1
            S1 = S1 + A(ii, jj) * x(jj);
        end
        S2 = 0;
        for jj = ii+1:n
            S2 = S2 + A(ii, jj) * x0(jj);
        end
        x(ii) = (1 - mu) * x0(ii) + mu/A(ii, ii) * (b(ii) - S1 - S2);
    end
    x0 = x;
    it = it + 1;
    x_vec = [x_vec x];
end

end

