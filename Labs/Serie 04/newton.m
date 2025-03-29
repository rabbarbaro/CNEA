function [xvect, it] = newton(x0, nmax, toll, fun, dfun, mol)

if nargin == 5
    mol = 1;
end

xvect = x0;
it = 0;
err = toll + 1;

while err > toll && it < nmax
    if dfun(x0) == 0
        error(' Arresto per azzeramento di dfun');
    else
        x = x0 - mol * fun(x0) / dfun(x0);
        err = abs(x - x0);
        x0 = x;
        it = it + 1;
        xvect = [xvect; x];
    end
end



end

