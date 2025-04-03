function K = K_SDP(A,tol,nmax,x0)

[n,m] = size(A);
if n ~= m
    error('Solo per matrici quadrate');
end

if nargin == 1
   tol = 1.e-06;   x0 = ones(n,1);   nmax = 100;
end

% iterazione zero fuori dal ciclo while
iter = 0;
y_max = x0/norm(x0);
err = tol + 1;
y_min = y_max;
K = 1;

% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(A);

while err > tol && abs(y_min'*A*y_min) ~= 0 && iter<nmax
   iter = iter + 1;
   x_max = A*y_max;
   y_max= x_max/norm(x_max);

   % risolvo Ax^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y_min);
   x_min=bksub(U,z);
   y_min= x_min/norm(x_min);

   Knew = (y_max'*A*y_max) / (y_min'*A*y_min);
   err = abs(Knew - K);
   K = Knew;
end

if (err <= tol)
     fprintf('eigpower converge in %d iterazioni \n', iter);
else
     fprintf('eigpower non converge in %d iterazioni. \n', iter)
end

return
